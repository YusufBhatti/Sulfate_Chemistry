! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_sl_helmholtz_inc_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SL_HELMHOLTZ_INC_MOD'

CONTAINS
SUBROUTINE eg_SL_Helmholtz_inc(                                 &
           exner_np1, p_star_np1, u_np1, v_np1, w_np1,          &
           etadot_np1, rho_np1, thetav_np1,                     &
           exner_prime_term, m_v_np1, m_cl_np1, m_cf_np1,       &
           m_r_np1,m_gr_np1, m_cf2_np1, cycleno, Ih, InnIts,    &
           pre_type, l_rel_tol, GCR_Diagnostics,                &
           R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,             &
           R_m_v_d, R_m_cl_d, R_m_cf_d,                         &
           R_m_r_d, R_m_gr_d, R_m_cf2_d, tol, tol_sc_fact,      &
           n_rows, row_length,                                  &
           rows, model_levels, S_u, S_v, S_w, S_thetav,         &
           S_m_v, S_m_cl, S_m_cf, S_m_cf2, S_m_r, S_m_gr        &
           ,psi_w_surf, psi_w_lid)

USE um_parvars,        ONLY: offx, offy, halo_i, halo_j,        &
                             datastart

USE yomhook,           ONLY: lhook, dr_hook
USE parkind1,          ONLY: jprb, jpim

USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels, &
                             xi3_at_theta=>r_theta_levels,     &
                             xi3_at_rho=>r_rho_levels,         &
                             xi3_at_u=>r_at_u,                 &
                             xi3_at_v=>r_at_v

USE timestep_mod,      ONLY: timestep
USE um_parcore,        ONLY: mype, nproc

USE ereport_mod,       ONLY: ereport
USE umPrintMgr
USE Field_Types

USE eg_add_pressure_inc_mod
USE eg_helm_rhs_inc_mod
USE eg_bicgstab_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE metric_terms_mod

USE gravity_mod
USE helmholtz_const_matrix_mod
USE coriolis_mod

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: model_type, mt_cyclic_lam

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

REAL,    PARAMETER                :: sc_err_min = 1.0e-20

INTEGER,            INTENT(IN)    :: row_length
INTEGER,            INTENT(IN)    :: rows
INTEGER,            INTENT(IN)    :: model_levels
INTEGER,            INTENT(IN)    :: n_rows
INTEGER,            INTENT(IN)    :: cycleno, InnIts, pre_type
REAL,               INTENT(IN)    :: Ih
REAL                              :: tol, tol_sc_fact

REAL ::                                                                  &
psi_w_surf(row_length,rows),                                             &
psi_w_lid (row_length,rows)

INTEGER, INTENT(IN) :: GCR_Diagnostics
                           ! Switch controlling diagnostic output.
                           ! 0 = none
                           ! 1 = initial and final residuals
                           ! 2 = all
                           ! 3 = iteration count processing

! Logical flags

LOGICAL, INTENT(IN) ::  l_rel_tol


! Input arrays from the departure point calculation

REAL, INTENT(IN) ::    R_u_d(udims_s%i_start:udims_s%i_end,     &
                             udims_s%j_start:udims_s%j_end,     &
                             udims_s%k_start:udims_s%k_end)

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


REAL, INTENT(INOUT) :: R_v_d(vdims_s%i_start:vdims_s%i_end,     &
                             vdims_s%j_start:vdims_s%j_end,     &
                             vdims_s%k_start:vdims_s%k_end)

! Pressure perturbation

REAL, INTENT(INOUT) :: exner_np1(pdims_s%i_start:pdims_s%i_end, &
                                 pdims_s%j_start:pdims_s%j_end, &
                                 pdims_s%k_start:pdims_s%k_end+1)


REAL, INTENT(INOUT) :: p_star_np1(pdims_s%i_start:pdims_s%i_end,&
                                  pdims_s%j_start:pdims_s%j_end)

! Fields at next time level

REAL, INTENT(INOUT) :: u_np1(udims_s%i_start:udims_s%i_end,     &
                             udims_s%j_start:udims_s%j_end,     &
                             udims_s%k_start:udims_s%k_end)

REAL, INTENT(INOUT) :: v_np1(vdims_s%i_start:vdims_s%i_end,     &
                             vdims_s%j_start:vdims_s%j_end,     &
                             vdims_s%k_start:vdims_s%k_end)

REAL, INTENT(INOUT) :: rho_np1(pdims_s%i_start:wdims_s%i_end,   &
                               pdims_s%j_start:pdims_s%j_end,   &
                               pdims_s%k_start:pdims_s%k_end)

REAL, INTENT(INOUT) ::  m_v_np1(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: m_cl_np1(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: m_cf_np1(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::  m_r_np1(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: m_gr_np1(tdims_s%i_start:wdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: m_cf2_np1(tdims_s%i_start:wdims_s%i_end, &
                                 tdims_s%j_start:tdims_s%j_end, &
                                 tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: thetav_np1(tdims_s%i_start:wdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) :: w_np1(wdims_s%i_start:wdims_s%i_end,     &
                             wdims_s%j_start:wdims_s%j_end,     &
                             wdims_s%k_start:wdims_s%k_end)

REAL, INTENT(INOUT) :: etadot_np1(wdims_s%i_start:wdims_s%i_end,&
                                  wdims_s%j_start:wdims_s%j_end,&
                                  wdims_s%k_start:wdims_s%k_end)

! Output array for plotting

REAL, INTENT(OUT)  ::                                           &
   exner_prime_term(pdims_s%i_start:wdims_s%i_end,              &
                    pdims_s%j_start:pdims_s%j_end,              &
                    pdims_s%k_start:pdims_s%k_end)

! Local arrays used for departure points

REAL                :: R_u_a(udims_s%i_start:udims_s%i_end,     &
                             udims_s%j_start:udims_s%j_end,     &
                             udims_s%k_start:udims_s%k_end)

REAL                :: R_v_a(vdims_s%i_start:vdims_s%i_end,     &
                             vdims_s%j_start:vdims_s%j_end,     &
                             vdims_s%k_start:vdims_s%k_end)

REAL                ::   R_w_a(wdims%i_start:wdims%i_end,       &
                               wdims%j_start:wdims%j_end,       &
                               wdims%k_start:wdims%k_end)

REAL                ::  R_theta_a(tdims_s%i_start:wdims_s%i_end,&
                                  tdims_s%j_start:tdims_s%j_end,&
                                  tdims_s%k_start:tdims_s%k_end)

REAL                :: R_rho_a(pdims_s%i_start:wdims_s%i_end,   &
                               pdims_s%j_start:pdims_s%j_end,   &
                               pdims_s%k_start:pdims_s%k_end)


REAL                ::  R_p_a(pdims_s%i_start:wdims_s%i_end,    &
                              pdims_s%j_start:pdims_s%j_end,    &
                              pdims_s%k_start:pdims_s%k_end)

REAL                :: R_etadot_a(wdims%i_start:wdims%i_end,    &
                                  wdims%j_start:wdims%j_end,    &
                                  wdims%k_start:wdims%k_end)
REAL ::                                                         &
                 rho_divu_out(pdims_s%i_start:wdims_s%i_end,    &
                              pdims_s%j_start:pdims_s%j_end,    &
                              pdims_s%k_start:pdims_s%k_end)

REAL ::                                                         &
                          rhs(pdims_s%i_start:wdims_s%i_end,    &
                              pdims_s%j_start:pdims_s%j_end,    &
                              pdims_s%k_start:pdims_s%k_end)

REAL ::                                                         &
       R_v_North_fixd(vdims_s%i_start:vdims_s%i_end,            &
                      vdims_s%k_start:vdims_s%k_end),           &
       R_v_South_fixd(vdims_s%i_start:vdims_s%i_end,            &
                      vdims_s%k_start:vdims_s%k_end),           &
        R_v_North_var(vdims_s%i_start:vdims_s%i_end,            &
                      vdims_s%k_start:vdims_s%k_end),           &
        R_v_South_var(vdims_s%i_start:vdims_s%i_end,            &
                      vdims_s%k_start:vdims_s%k_end)

INTEGER                              :: i, j, k, iter
REAL                                 :: del_rho
REAL                                 :: ex_err(2)
REAL                                 :: init_err, fin_err
INTEGER                              :: no_its
INTEGER                              :: istat1

REAL :: RHS_err, tol_save

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

LOGICAL, PARAMETER  :: l_inc= .TRUE.   ! is this an increment routine?

CHARACTER(LEN=errormessagelength)            :: Cmessage
CHARACTER(LEN=15)             :: Routine = 'eg_sl_helmholtz'
INTEGER                       :: icode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SL_HELMHOLTZ_INC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

icode = 0

IF ( pre_type < 0 .OR. pre_type > 4 ) THEN
  icode = 5
  cmessage = 'Invalid preconditioner type'
  CALL Ereport(Routine,icode,cmessage)
END IF

tol_save = tol

del_rho = 1.0

! Calculate pressure deviation from reference profile

IF (model_type == mt_cyclic_lam) THEN
  icode    = 6
  cmessage = 'LAM not available.'
  CALL Ereport(Routine,icode,cmessage)
END IF

IF ( mype == 0 .AND. GCR_Diagnostics >0 ) THEN
  WRITE(umMessage,'(A44)') "********************************************"
  CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
  WRITE(umMessage,'(A44)') "*   iLinear solve for Helmholtz problem    *"
  CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
  WRITE(umMessage,'(A44)') "*   ====================================   *"
  CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
END IF

inner:DO iter = 1, InnIts

  DO k = 1, model_levels
    exner_prime_term(:,:,k) = 0.0
  END DO

  ! Calculate arrival point quantities1

  CALL eg_helm_rhs_inc(R_u_a, R_v_a, R_w_a,                     &
                      R_p_a, R_theta_a, R_rho_a,                &
                      R_etadot_a, rho_divu_out,                 &
                      R_m_v_d, R_m_cl_d, R_m_cf_d,              &
                      R_m_r_d, R_m_gr_d, R_m_cf2_d,             &
                      u_np1, v_np1, w_np1,                      &
                      exner_np1, p_star_np1,                    &
                      rho_np1, thetav_np1, etadot_np1,          &
                      m_v_np1, m_cl_np1, m_cf_np1,              &
                      m_r_np1, m_gr_np1, m_cf2_np1,             &
                      R_u_d, R_v_d, R_w_d, R_theta_d, R_rho_d,  &
                      Ih, row_length, rows, n_rows,             &
                      model_levels, iter,                       &
                      S_u,S_v,S_w,S_thetav,S_m_v,S_m_cl,S_m_cf, &
                      S_m_cf2,S_m_r,S_m_gr,                     &
                      rhs,psi_w_surf, psi_w_lid)


  IF ( GCR_Diagnostics > 0 ) THEN
    IF ( PrintStatus >= PrStatus_Normal) THEN
      RHS_err = -1.0
      DO k = 1, model_levels
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            ! Compute error
            RHS_err = MAX(RHS_err,ABS(rhs(i,j,k)))
          END DO
        END DO
      END DO
    CALL gc_rmax(1,nproc,istat1,RHS_err)
    END IF
  END IF

  ! Solve the linear system
  CALL eg_bicgstab(exner_prime_term,rhs,tol,pre_type,l_rel_tol, &
         row_length, rows, n_rows, model_levels,sc_err_min,     &
         init_err, fin_err, no_its,l_inc)

  IF ( GCR_Diagnostics > 0 ) THEN
    IF ( PrintStatus >= PrStatus_Normal) THEN

      IF ( mype == 0 ) THEN 
        WRITE(umMessage,"(A20,E13.6,10X,A1)")                           &
              "*   Error in RHS    ", RHS_err, "*"
        CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
      END IF
    END IF

    ex_err(1) = - MINVAL(exner_prime_term)
    ex_err(2) =   MAXVAL(exner_prime_term)
    CALL gc_rmax(2,nproc,istat1,ex_err)

    IF ( mype == 0 ) THEN
      WRITE(umMessage,"(A24,I3,A1,I3,12X,A1)")                          &
            "* Outer/Inner iteration ", cycleno, "/", iter, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
      WRITE(umMessage,"(A34,I4,5X,A1)")                                 &
            "* No. Of linear solver iterations ", no_its, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
      WRITE(umMessage,"(A16,ES13.6,15X,A1)")                             &
            "* Initial error ", init_err, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
      WRITE(umMessage,"(A16,ES13.6,15X,A1)")                             &
            "*   Final error ", fin_err, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
      WRITE(umMessage,"(A20,E13.6,10X,A1)")                             &
            "*   Error in RHS    ", RHS_err, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
      WRITE(umMessage,"(A20,E13.6,10X,A1)")                             &
            "*   Solver tol      ", tol, "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')

      WRITE(umMessage,"(A20,2E13.6,10X,A1)")                            &
            "*   Exner err      ",-ex_err(1),ex_err(2), "*"
      CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
    END IF
  END IF

  ! Calculate time level n+1 quantities
  CALL eg_add_pressure_inc(                                            &
            exner_prime_term,exner_np1,p_star_np1,                     &
            R_v_South_var, R_v_North_var,                              &
            u_np1, v_np1, w_np1, thetav_np1, etadot_np1, rho_np1,      &
            m_v_np1, m_cl_np1, m_cf_np1, m_r_np1, m_gr_np1, m_cf2_np1, &
            R_w_d,R_u_a, R_v_a, R_w_a, R_theta_a, R_etadot_a, R_p_a,   &
            Ih, n_rows, row_length, rows, model_levels,                &
            psi_w_surf, psi_w_lid)

  tol = tol/tol_sc_fact

END DO inner

IF ( mype == 0 .AND. GCR_Diagnostics > 0 ) THEN
  WRITE(umMessage,'(A44)') "********************************************"
  CALL umPrint(umMessage,src='eg_sl_helmholtz_inc')
  CALL umPrint(" ",src='eg_sl_helmholtz_inc')
END IF
tol = tol_save

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_SL_Helmholtz_inc
END MODULE eg_sl_helmholtz_inc_mod
