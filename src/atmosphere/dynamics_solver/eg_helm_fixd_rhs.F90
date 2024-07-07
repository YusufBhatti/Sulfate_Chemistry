! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_helm_fixd_rhs_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_HELM_FIXD_RHS_MOD'

CONTAINS
SUBROUTINE eg_helm_fixd_rhs(rhs,eta_rho_levels,row_length, rows,&
                            model_levels,offx, offy)

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE eg_helmholtz_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE Field_Types
USE metric_terms_mod
USE fields_rhs_mod

USE helmholtz_const_matrix_mod
USE coriolis_mod

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
INTEGER,  INTENT(IN)    :: row_length, rows,  model_levels

REAL, INTENT(IN) ::   eta_rho_levels(model_levels)

REAL, INTENT(OUT)        ::                                              &
  rhs(1-offx:row_length+offx,1-offy:rows+offy,model_levels)


! Local temporary variables

REAL    :: ajm(1-offx:row_length+offx,1-offy:rows+offy,0:model_levels-1)
REAL    :: rdxi1, rdxi2, rdxi3
INTEGER :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_HELM_FIXD_RHS'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate RHS^n eqn(9.45) of EG3.01

!$OMP PARALLEL  DEFAULT(NONE) PRIVATE(i,j,k,rdxi1,rdxi2,rdxi3)    &
!$OMP& SHARED(model_levels,pdims,xi2_v,xi1_u,Hlm_Ek,Hlm_Fk,ajm,   &
!$OMP& rhs,rho_ref_pro,intw_w2rho,R_theta_d,thetav_ref_pro,HM_p,  &
!$OMP& R_rho_d,HM_vol,HM_rhox,R_u_d,HM_rhoy,R_v_d,eta_rho_levels, &
!$OMP& Hlm_Ck,R_w_d,HM_w,exner_ref_pro)

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
    DO i = pdims%i_start, pdims%i_end
      rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

      rhs(i,j,k) = -rho_ref_pro(i,j,k)*(                           &
                    intw_w2rho(k,1)*R_theta_d(i,j,k)               &
                                   /thetav_ref_pro(i,j,k)          &
                   +intw_w2rho(k,2)*R_theta_d(i,j,k-1)             &
                                   /thetav_ref_pro(i,j,k-1) )      &
                   /HM_p(i,j,k) - R_rho_d(i,j,k)


      rhs(i,j,k) = rhs(i,j,k) + HM_vol(i,j,k)*(                    &
                              ( HM_rhox(i,j,k)*R_u_d(i,j,k)        &
                              - HM_rhox(i-1,j,k)*R_u_d(i-1,j,k)    &
                              )*rdxi1                              &
                             +( HM_rhoy(i,j,k)*R_v_d(i,j,k)        &
                              - HM_rhoy(i,j-1,k)*R_v_d(i,j-1,k)    &
                              )*rdxi2   )
    END DO
  END DO
END DO
!$OMP END DO

! These are the terms arising from D_1 in (9.47)

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ajm(i,j,0) = 0.0
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
  rdxi3 = 1.0/( eta_rho_levels(k+1) - eta_rho_levels(k))
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      ! Use same definintion of D_1(X) as in LHS :
      !            D_1(X) = E_k.X_k + F_k.X_(k-1)

      ! NB. HM_vol is subsumed into Hlm_Ek and Hlm_Fk, see eg_set_helm_lhs.

      ajm(i,j,k) = Hlm_Ck(i,j,k) *                                            &
           ( R_w_d(i,j,k) - HM_w(i,j,k) * R_theta_d(i,j,k) *          &
             ( exner_ref_pro(i,j,k+1) - exner_ref_pro(i,j,k) )        &
             *rdxi3/thetav_ref_pro(i,j,k) )

    END DO
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rhs(i,j,k) = rhs(i,j,k) + Hlm_Ek(i,j,k)*ajm(i,j,k)              &
                              + Hlm_Fk(i,j,k)*ajm(i,j,k-1)
    END DO
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

! Fix up top

k = model_levels
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    rhs(i,j,k) = rhs(i,j,k) + Hlm_Fk(i,j,k)*ajm(i,j,k-1)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_helm_fixd_rhs
END MODULE eg_helm_fixd_rhs_mod
