! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_helm_var_rhs_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_HELM_VAR_RHS_MOD'

CONTAINS
SUBROUTINE eg_helm_var_rhs(rhs,Rn,Ih, eta_rho_levels,           &
       exner_prime,                                             &
       R_u_a, R_v_a, R_w_a, R_theta_a, R_rho_a, R_p_a,          &
       R_etadot,row_length, rows, n_rows, model_levels,         &
       offx, offy)

USE eg_vert_damp_mod, ONLY: mu_w
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE eg_helmholtz_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE um_parvars, ONLY: at_extremity
USE UM_ParParams
USE Field_Types
USE helmholtz_const_matrix_mod
USE coriolis_mod

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE

!
! Description: Code to calculate the Fixed RHS terms
!              in the Helmholtz problem
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

INTEGER,  INTENT(IN)    :: offx, offy
INTEGER,  INTENT(IN)    :: row_length, rows, n_rows, model_levels

! Hydrostatic switch

REAL,      INTENT(IN)    :: Ih

REAL, INTENT(IN) ::   eta_rho_levels(model_levels)

REAL, INTENT(IN) :: exner_prime(1-offx:row_length+offx,                  &
                                1-offy:rows+offy,model_levels)

REAL, INTENT(OUT)   ::                                                   &
    rhs(1-offx:row_length+offx,1-offy:rows+offy,model_levels)

REAL, INTENT(IN)    ::                                                   &
     Rn(1-offx:row_length+offx,1-offy:rows+offy,model_levels)


REAL, INTENT(IN) ::                                                      &
 R_u_a(-offx:row_length+offx-1,1-offy:rows+offy,model_levels),           &
 R_v_a(1-offx:row_length+offx,-offy:n_rows+offy-1,model_levels),         &
 R_w_a(row_length,rows,0:model_levels),                                  &
 R_theta_a(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels),      &
 R_rho_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels),          &
 R_p_a(1-offx:row_length+offx,1-offy:rows+offy,model_levels),            &
 R_etadot(row_length,rows,0:model_levels)


! Local temporary variables

REAL    :: Rk(1-offx:row_length+offx,1-offy:rows+offy,model_levels)
REAL    :: rdxi1(pdims%i_start:pdims%i_end)
REAL    :: rdxi2(pdims%j_start:pdims%j_end)
REAL    :: rdxi3
INTEGER :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_HELM_VAR_RHS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = pdims%i_start, pdims%i_end
  rdxi1(i) = 1.0/(xi1_u(i) - xi1_u(i-1))
END DO

DO j = pdims%j_start, pdims%j_end
  rdxi2(j) = 1.0/(xi2_v(j) - xi2_v(j-1))
END DO

! Calculate total RHS = RHS^n + RHS^star, eqn(9.44) of EG3.01

!$OMP  PARALLEL DEFAULT(NONE)                                     &
!$OMP& PRIVATE(k,j,i,rdxi3)                                       &
!$OMP& SHARED(model_levels,pdims,rhs,Rn,rho_ref_pro,              &
!$OMP& HM_p,R_p_a,intw_w2rho,R_theta_a,thetav_ref_pro,R_rho_a,    &
!$OMP& HM_vol,HM_rhox,R_u_a,HM_rhoy,R_v_a,Rk,Hlm_Ck,R_w_a,Ih,     &
!$OMP& mu_w,R_etadot,HM_w,exner_ref_pro, eta_rho_levels,Hlm_Lp,   &
!$OMP& Hlm_Ek,Hlm_Fk,rdxi1,rdxi2)
!$OMP  DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    Rk(i,j,1) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
  rdxi3 = 1.0/(eta_rho_levels(k+1)-eta_rho_levels(k))
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rhs(i,j,k) = Rn(i,j,k) - (rho_ref_pro(i,j,k)/HM_p(i,j,k)) *  &
                   ( R_p_a(i,j,k) +                                &
                   ( intw_w2rho(k,1)*R_theta_a(i,j,k) /            &
                                    thetav_ref_pro(i,j,k) +        &
                     intw_w2rho(k,2)*R_theta_a(i,j,k-1) /          &
                                    thetav_ref_pro(i,j,k-1) ))     &
                   - R_rho_a(i,j,k)
      rhs(i,j,k) = rhs(i,j,k) + HM_vol(i,j,k)*(                    &
                           ( HM_rhox(i,j,k)*R_u_a(i,j,k)           &
                           - HM_rhox(i-1,j,k)*R_u_a(i-1,j,k)       &
                           )*rdxi1(i)                              &
                          +( HM_rhoy(i,j,k)*R_v_a(i,j,k)           &
                           - HM_rhoy(i,j-1,k)*R_v_a(i,j-1,k)       &
                           )*rdxi2(j) )

      ! Use same definintion of D_1(R) as in LHS :
      !            D_1(X) = E_k.R_k + F_k.R_(k-1)

      Rk(i,j,k+1) = Hlm_Ck(i,j,k) *  ( R_w_a(i,j,k) +              &
                    (Ih+mu_w(i,j,k))*R_etadot(i,j,k) -             &
                    HM_w(i,j,k)*( exner_ref_pro(i,j,k+1)           &
                                 -exner_ref_pro(i,j,k  )) *        &
                   R_theta_a(i,j,k)*rdxi3/thetav_ref_pro(i,j,k) )

    END DO
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end


      rhs(i,j,k) = rhs(i,j,k) + Hlm_Ek(i,j,k)*Rk(i,j,k+1)        &
                                + Hlm_Fk(i,j,k)*Rk(i,j,k)
! Now rescale by the diagonal of linear operator: Hlm_Lp contains the
! reciprocal of the diagonal, see eg_set_helm_lhs.
      rhs(i,j,k) = rhs(i,j,k)*Hlm_Lp(i,j,k)

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

k = model_levels
!$OMP  DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    rhs(i,j,k) = Rn(i,j,k) - (rho_ref_pro(i,j,k)/HM_p(i,j,k)) *  &
                 ( R_p_a(i,j,k) +                                &
                 ( intw_w2rho(k,1)*R_theta_a(i,j,k) /            &
                                  thetav_ref_pro(i,j,k) +        &
                   intw_w2rho(k,2)*R_theta_a(i,j,k-1) /          &
                                  thetav_ref_pro(i,j,k-1) ))     &
                 - R_rho_a(i,j,k)

    rhs(i,j,k) = rhs(i,j,k) + HM_vol(i,j,k)*(                    &
                         ( HM_rhox(i,j,k)*R_u_a(i,j,k)           &
                         - HM_rhox(i-1,j,k)*R_u_a(i-1,j,k)       &
                         )*rdxi1(i)                              &
                        +( HM_rhoy(i,j,k)*R_v_a(i,j,k)           &
                         - HM_rhoy(i,j-1,k)*R_v_a(i,j-1,k)       &
                         )*rdxi2(j) )
! Fix up top
    rhs(i,j,k) = rhs(i,j,k) + Hlm_Fk(i,j,k)*Rk(i,j,k)

! Now rescale by the diagonal of linear operator: Hlm_Lp contains the
! reciprocal of the diagonal, see eg_set_helm_lhs.
    rhs(i,j,k) = rhs(i,j,k)*Hlm_Lp(i,j,k)

  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! Exner fixed on boundary (exner prime =0)

IF (model_type == mt_lam) THEN
  i = pdims%i_end - 1         ! one less p-point in LAM's
  IF ( at_extremity(PWest) ) THEN
    rhs(1,:,:) = exner_prime(1,:,:)
  END IF
  IF ( at_extremity(PEast) ) THEN
    rhs(i,:,:)   = exner_prime(i,:,:)
    rhs(i+1,:,:) = exner_prime(i+1,:,:)
  END IF

  j = pdims%j_end
  IF ( at_extremity(PSouth) ) THEN
    rhs(:,1,:) = exner_prime(:,1,:)
  END IF
  IF ( at_extremity(PNorth) ) THEN
    rhs(:,j,:) = exner_prime(:,j,:)
  END IF
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_helm_var_rhs
END MODULE eg_helm_var_rhs_mod
