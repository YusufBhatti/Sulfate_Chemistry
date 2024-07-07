! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_sl_helmholtz_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SL_HELMHOLTZ_MOD'

CONTAINS
SUBROUTINE eg_sl_helmholtz(                                              &
           exner_np1, p_star_np1, u_np1, v_np1, w_np1,                   &
           etadot_np1, rho_np1, thetav_np1, exner_prime_term,            &
           m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1,    &
           cycleno, Ih,  InnIts, pre_type,l_rel_tol, GCR_Diagnostics,    &
           R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,L_eliminate_rho,      &
           R_m_v_d, R_m_cl_d, R_m_cf_d,R_m_r_d, R_m_gr_d, R_m_cf2_d,     &
           tol, tol_sc_fact,n_rows, row_length, rows, model_levels,      &
           S_u,S_v,S_w,R_theta_a,S_m_v,S_m_cl,S_m_cf,S_m_cf2,            &
           S_m_r,S_m_gr,psi_w_surf, psi_w_lid)

USE eg_alpha_mod,      ONLY: alpha_w

USE um_parvars,        ONLY: offx, offy, at_extremity
USE um_parparams

USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,          &
                              xi3_at_theta=>r_theta_levels,              &
                              xi3_at_rho=>r_rho_levels, xi3_at_u=>r_at_u,&
                              xi3_at_v=>r_at_v

USE timestep_mod,      ONLY: timestep
USE um_parcore,        ONLY: mype

USE ereport_mod, ONLY: ereport
USE Field_Types

USE eg_helm_var_rhs_mod
USE eg_helm_fixd_rhs_mod
USE eg_add_pressure_mod
USE eg_helm_rhs_star_mod
USE eg_bicgstab_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE metric_terms_mod

USE gravity_mod
USE helmholtz_const_matrix_mod
USE coriolis_mod

USE eg_parameters_mod, ONLY: total_conv_inner

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_cyclic_lam

USE atm_step_local, ONLY: L_print_L2norms
USE eg_helm_norm_mod
USE moist_norm_mod
USE lam_config_inputs_mod, ONLY: n_rims_to_do
USE turb_diff_mod, ONLY: norm_lev_start, norm_lev_end
USE umPrintMgr, ONLY: umMessage, umPrint

IMPLICIT NONE

!
! Description: Solves the nonlinear Helmholtz problem
!
! Method: ENDGame formulation version 3.01
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

! Array dimensions


REAL, PARAMETER :: sc_err_min = 1.0e-4

INTEGER, INTENT(IN) :: row_length
INTEGER, INTENT(IN) :: rows
INTEGER, INTENT(IN) :: model_levels
INTEGER, INTENT(IN) :: n_rows
INTEGER             :: InnIts, pre_type
INTEGER, INTENT(IN) :: cycleno, GCR_Diagnostics
                           ! Switch controlling diagnostic output.
                           ! 0 = none
                           ! 1 = initial and final residuals
                           ! 2 = all
                           ! 3 = iteration count processing

REAL, INTENT(IN) :: Ih
REAL :: tol, tol_sc_fact
REAL :: psi_w_surf(row_length,rows), psi_w_lid (row_length,rows)

! Logical flags

LOGICAL, INTENT(IN) ::  l_rel_tol

! Input arrays from the departure point calculation

REAL, INTENT(IN) ::                                                     &
 R_u_d(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),          &
 R_w_d(row_length,rows,0:model_levels),                                 &
 R_theta_d(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),     &
 R_rho_d(1-offx:row_length+offx,1-offy:rows+offy,model_levels),         &
 R_m_v_d(row_length,rows,0:model_levels),                               &
 R_m_cl_d(row_length,rows,0:model_levels),                              &
 R_m_cf_d(row_length,rows,0:model_levels),                              &
 R_m_r_d(row_length,rows,0:model_levels),                               &
 R_m_gr_d(row_length,rows,0:model_levels),                              &
 R_m_cf2_d(row_length,rows,0:model_levels)

! Modified if running a cyclic LAM
REAL, INTENT(INOUT) ::                                                  &
 R_v_d(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels)

LOGICAL, PARAMETER  :: L_accel_convergence = .FALSE.


LOGICAL, INTENT(IN) :: L_eliminate_rho

! Pressure perturbation

REAL, INTENT(INOUT)  ::                                                  &
 exner_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels+1)

REAL ::                                                                  &
 exner0(1:row_length,1:rows,model_levels+1)

REAL, INTENT(INOUT)   ::                                                 &
      p_star_np1(1-offx:row_length+offx,1-offy:rows+offy)

! Fields at next time level

REAL, INTENT(INOUT)   ::                                                 &
      u_np1(-offx:row_length-1+offx,1-offy:rows+offy,model_levels),      &
      v_np1(1-offx:row_length+offx,-offy:n_rows-1+offy,model_levels),    &
      rho_np1(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(INOUT)   ::                                                 &
      m_v_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
     m_cl_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
     m_cf_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
      m_r_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
     m_gr_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
    m_cf2_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
   thetav_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
        w_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),   &
   etadot_np1(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

! Output array for plotting

REAL, INTENT(OUT)   ::                                                   &
      exner_prime_term(1-offx:row_length+offx,1-offy:rows+offy,          &
                       model_levels)

! Local arrays used for departure points

REAL                ::                                                   &
     R_u_a(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),       &
     R_v_a(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels),     &
     R_w_a(row_length,rows,0:model_levels),                              &
   R_rho_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
  R_etadot(row_length,rows,0:model_levels),                              &
     R_p_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL               ::                                                   &
      rhs(1-offx:row_length+offx,1-offy:rows+offy,model_levels),        &
       Rn(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

INTEGER :: k, iter
REAL    :: del_rho
REAL    :: ex_extrema(2)
REAL    :: init_err, fin_err
INTEGER :: no_its
INTEGER :: istat1, istat2

!     Fast Physics source terms ( R_theta_a = S_theta )

REAL, INTENT(INOUT) ::                                                   &
 R_theta_a(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels)

REAL, INTENT(IN) ::                                                      &
       S_u( -offx:row_length-1+offx,1-offy:  rows  +offy,model_levels),  &
       S_v(1-offx:row_length  +offx, -offy:n_rows-1+offy,model_levels),  &
       S_w(row_length,rows,0:model_levels),                              &
     S_m_v(row_length,rows,0:model_levels),                              &
    S_m_cl(row_length,rows,0:model_levels),                              &
    S_m_cf(row_length,rows,0:model_levels),                              &
     S_m_r(row_length,rows,0:model_levels),                              &
    S_m_gr(row_length,rows,0:model_levels),                              &
   S_m_cf2(row_length,rows,0:model_levels)

CHARACTER(LEN=errormessagelength)            :: Cmessage
CHARACTER(LEN=*), PARAMETER             :: RoutineName = 'EG_SL_HELMHOLTZ'
INTEGER                       :: icode
LOGICAL, PARAMETER  :: l_inc= .FALSE.   ! is this a increment routine

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0
IF ( pre_type < 0 .OR. pre_type > 4 ) THEN
  icode = 5
  cmessage = 'Invalid preconditioner type'
  CALL Ereport(RoutineName,icode,cmessage)
END IF

!     del_rho = 0.0
!     IF ( .NOT. L_SLICE ) del_rho = 1.0
del_rho = 1.0

! Adjust R_v_d at North and South boundaries for use in a cyclic LAM

IF (model_type == mt_cyclic_lam) THEN
  IF ( at_extremity(Psouth) ) THEN
    R_v_d(:,vdims%j_start,:) = 0.0
  END IF
  IF ( at_extremity(Pnorth) ) THEN
    R_v_d(:,vdims%j_end,:)   = 0.0
  END IF
END IF

! Calculate fixed part of Right-Hand-Side
CALL eg_helm_fixd_rhs(Rn,eta_rho_levels,row_length, rows, model_levels,&
                      offx, offy)

! Calculate pressure deviation from reference profile and initialise rhs
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                         &
!$OMP          PRIVATE(k) SHARED(rhs, model_levels, exner_prime_term,    &
!$OMP                            exner_np1, exner_ref_pro)
DO k = 1, model_levels
  rhs(:,:,k) = 0.0         ! mainly to ensure the halos are set to zero
  exner_prime_term(:,:,k) = ( exner_np1(:,:,k) - exner_ref_pro(:,:,k) )
END DO
!$OMP END PARALLEL DO

IF ( mype == 0 .AND. GCR_Diagnostics > 0 ) THEN
  WRITE(umMessage,"(A44)") "********************************************"
  CALL umPrint(umMessage,src='eg_sl_helmholtz')
  WRITE(umMessage,"(A44)") "*    Linear solve for Helmholtz problem    *"
  CALL umPrint(umMessage,src='eg_sl_helmholtz')
  IF ( GCR_Diagnostics == 1 ) THEN
    WRITE(umMessage,"(A44)")                                           &
          "* Outer Inner Iterations InitialError      *"
  ELSE
    WRITE(umMessage,"(A44)")                                           &
          "*   ====================================   *"
  END IF
  CALL umPrint(umMessage,src='eg_sl_helmholtz')
END IF


IF ( total_conv_inner) InnIts=20

DO iter = 1, InnIts

  ! Calculate arrival point quantities

  CALL eg_helm_rhs_star(R_u_a, R_v_a, R_w_a, R_etadot,                  &
                     R_p_a, R_theta_a, R_rho_a,                         &
                     R_m_v_d, R_m_cl_d, R_m_cf_d,                       &
                     R_m_r_d, R_m_gr_d, R_m_cf2_d,                      &
                     u_np1, v_np1, w_np1, exner_np1,                    &
                     p_star_np1, exner_prime_term, rho_np1, thetav_np1, &
                     etadot_np1, m_v_np1, m_cl_np1, m_cf_np1,           &
                     m_r_np1, m_gr_np1, m_cf2_np1, R_w_d, del_rho, Ih,  &
                     timestep, alpha_w,row_length, rows, n_rows,        &
                     model_levels,                                      &
                     iter, S_u,S_v,S_w,S_m_v,S_m_cl,S_m_cf,             &
                     S_m_cf2,S_m_r,S_m_gr,psi_w_surf, psi_w_lid)

  ! Update RHS of Helmholtz operator to include "star" values

  CALL eg_helm_var_rhs(rhs,Rn,Ih, eta_rho_levels,              &
      exner_prime_term,                                        &
      R_u_a, R_v_a, R_w_a, R_theta_a, R_rho_a, R_p_a,          &
      R_etadot,row_length, rows, n_rows, model_levels,         &
      offx, offy)

  ! Solve the linear system

#if defined(C_DP_HLM)
           ! The default, double precision solver
  IF ( mype == 0 .AND. GCR_Diagnostics > 1 ) THEN
    WRITE(umMessage,'(A44)')                                          &
          "*    Double precision solver               *"
    CALL umPrint(umMessage,src='eg_sl_helmholtz')
  END IF

  CALL eg_bicgstab(exner_prime_term,rhs,tol,pre_type,          &
       l_rel_tol, row_length, rows, n_rows, model_levels, &
       sc_err_min, init_err, fin_err, no_its, l_inc)
#else
  IF (tol <= 1e-5) THEN
     ! Safety first, we're too close to the bone
    cmessage = &
         'mixed precision solver not sufficient for accuracy of solver'
    icode = 3
    CALL Ereport(RoutineName,icode,cmessage)
  END IF

  IF ( mype == 0 .AND. GCR_Diagnostics > 1 ) THEN
    WRITE(umMessage,'(A44)')                                          &
          "*   Mixed precision solver                 *"
    CALL umPrint(umMessage,src='eg_sl_helmholtz')
  END IF
  ! call the mixed precision solver
  CALL eg_bicgstab_mixed_prec(exner_prime_term,rhs,tol,pre_type,   &
       l_rel_tol, row_length, rows, n_rows, model_levels, &
       sc_err_min, init_err, fin_err, no_its, l_inc)
#endif

  ! Calculate time level n+1 quantities

  CALL eg_add_pressure(exner_prime_term, exner_np1, p_star_np1,      &
          u_np1, v_np1, w_np1, thetav_np1, etadot_np1, rho_np1,      &
          m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1, &
          R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,                   &
          R_u_a, R_v_a, R_w_a, R_theta_a, R_rho_a,                   &
          R_etadot, R_p_a,                                           &
          alpha_w, timestep,  l_eliminate_rho,                       &
          Ih,  offx, offy,n_rows, row_length, rows, model_levels)


  IF ( mype == 0 .AND. GCR_Diagnostics > 0 ) THEN
    IF ( GCR_Diagnostics == 1 ) THEN
      WRITE(umMessage,"(A3,I3,3X,I3,5X,I4,4X,E13.6,5X,A1)")           &
            "*  ", cycleno, iter, no_its, init_err, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz')
    ELSE
      WRITE(umMessage,"(A24,I3,A1,I3,12X,A1)")                        &
            "* Outer/Inner iteration ", cycleno, "/", iter, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz')
      WRITE(umMessage,"(A34,I4,5X,A1)")                               &
            "* No. Of linear solver iterations ", no_its, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz')
      WRITE(umMessage,"(A16,E13.6,14X,A1)")                           &
            "* Initial error ", init_err, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz')
      WRITE(umMessage,"(A16,E13.6,14X,A1)")                           &
            "*   Final error ", fin_err, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz')
    END IF
  END IF

  tol = tol/tol_sc_fact
END DO

IF ( mype == 0 .AND. GCR_Diagnostics > 0) THEN
  WRITE(umMessage,'(A44)') "********************************************"
  CALL umPrint(umMessage,src='eg_sl_helmholtz')
  CALL umPrint(" ",src='eg_sl_helmholtz')
END IF

IF ( L_print_L2norms ) THEN
  WRITE(umMessage,'(A, I2, A)') '  **   cycleno = ', cycleno            &
                       , '  **  L2 norms after eg_sl_helmholtz   **'
  CALL umPrint(umMessage,src='EG_SL_HELMHOLTZ')
  CALL eg_helm_norm(                                                    &
                    norm_lev_start, norm_lev_end, n_rims_to_do,         &
                    exner_np1, p_star_np1, u_np1, v_np1, w_np1,         &
                    etadot_np1, rho_np1, thetav_np1, exner_prime_term,  &
                    .TRUE., .FALSE., .FALSE. )
  CALL moist_norm(                                                      &
                  norm_lev_start, norm_lev_end, n_rims_to_do,           &
                  m_v_np1, m_cl_np1, m_cf_np1,                          &
                  m_r_np1, m_gr_np1, m_cf2_np1,                         &
                  .TRUE., .FALSE., .FALSE. )
END IF !  L_print_L2norms

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_SL_Helmholtz
END MODULE eg_sl_helmholtz_mod
