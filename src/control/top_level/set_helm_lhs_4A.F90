! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE  eg_set_helm_lhs_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_SET_HELM_LHS_MOD'

CONTAINS
SUBROUTINE eg_set_helm_lhs(                                              &
       row_length, rows, n_rows, model_levels, ih)

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE eg_vert_damp_mod,      ONLY: mu_w
USE level_heights_mod,     ONLY: eta_theta_levels, eta_rho_levels
USE dynamics_input_mod,    ONLY: l_inc_solver
USE eg_helmholtz_mod
USE horiz_grid_mod
USE ref_pro_mod
USE atm_fields_bounds_mod
USE metric_terms_mod
USE Field_Types
USE helmholtz_const_matrix_mod
USE UM_ParVars
USE UM_ParParams, ONLY: nodomain, peast, pwest, pnorth, psouth

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE
!
! Description:calculates the constant
!              part of the Helmholtz coefficient matrix.
!
!
! Method: Appendix F, ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments



INTEGER, INTENT(IN) :: row_length, rows, n_rows, model_levels

REAL, INTENT(IN) :: ih

! Local variables

REAL    :: c, t, cp, tp, rdxi1, rdxi2, rdxi3
REAL    :: rd_theta(model_levels), rd_rho(2:model_levels)
INTEGER :: i, j, k, k_end, mid_point
INTEGER :: l, n, m            ! more looping/indexing

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_SET_HELM_LHS'


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

k_end = model_levels
DO k = 2, model_levels
  rd_theta(k) = 1.0/(eta_theta_levels(k)-eta_theta_levels(k-1))
  rd_rho(k)   = 1.0/(eta_rho_levels(k)-eta_rho_levels(k-1))
END DO
rd_theta(1) = 1.0/(eta_theta_levels(1)-eta_theta_levels(0))
t = (eta_rho_levels(1)-eta_theta_levels(0))
c = (eta_theta_levels(k_end)-eta_rho_levels(k_end))
mid_point = pdims%i_start + (pdims%i_end - pdims%i_start)/2 + 1


! Precalculate some useful constants for reuse in the
! RHS of the Helmholtz equation
! In (9.40) we  D1_k = E_k*P_k + F_k*P_(k-1)
! NOTE : (a) We've absorbed the common HM_vol factor into the definitions
!            of Ek and Fk.
!        (b) Version of model using (1/rho)*grad(p) had A_k and B_k
!            coefficients for a more comlicated D2 operator.


!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,j,k,l,n,cp,tp,rdxi1,rdxi2)   
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      cp = HM_etadot(i,j,k)*(Ih+mu_w(i,j,k))                         &
            -hm_w(i,j,k)*hm_theta(i,j,k)*hm_b(i,j,k)

      hlm_ck(i,j,k) = 1.0/cp

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    hlm_ck(i,j,0)            = 0.0
    hlm_ck(i,j,model_levels) = 0.0
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 2, model_levels-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      cp            = rho_ref_pro(i,j,k)/hm_p(i,j,k)
      tp            = hm_vol(i,j,k)*rd_theta(k)

      hlm_ek(i,j,k) = hm_rhoz(i,j,k)*tp                          &
                     +cp*intw_w2rho(k,1)*hm_theta(i,j,k)         &
                       /thetav_ref_pro(i,j,k)

      hlm_fk(i,j,k) = -hm_rhoz(i,j,k-1)*tp                       &
                      +cp*intw_w2rho(k,2)*hm_theta(i,j,k-1)      &
                        /thetav_ref_pro(i,j,k-1)

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

! Top and bottom values for Ek and Fk

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    hlm_ek(i,j,1)  = ( hm_rhoz(i,j,1)*hm_vol(i,j,1) +            &
                       rho_ref_pro(i,j,1) * hm_theta(i,j,1) *    &
                       t/(hm_p(i,j,1)*thetav_ref_pro(i,j,1)))*   &
                       rd_theta(1)
    hlm_fk(i,j,1) = 0.0

    hlm_ek(i,j,k_end) = 0.0
    hlm_fk(i,j,k_end) = (-hm_rhoz(i,j,k_end-1)*hm_vol(i,j,k_end) +     &
                     rho_ref_pro(i,j,k_end) *hm_theta(i,j,k_end-1) *   &
                     c/(hm_p(i,j,k_end)*thetav_ref_pro(i,j,k_end-1)))* &
                     rd_theta(k_end)

  END DO
END DO
!$OMP END DO

! Calculate the constant part of the Helmholtz matrix
! Simplest way is to build it in pieces!

!$OMP DO SCHEDULE(STATIC)
DO k=1, model_levels
  DO j=pdims%j_start, pdims%j_end
    rdxi2 = 1.0/(xi2_v(j) - xi2_v(j-1))
    DO i=pdims%i_start, pdims%i_end
      rdxi1 = 1.0/(xi1_u(i) - xi1_u(i-1))

      ! Main diagonal entry

      hlm_lp(i,j,k) = -rho_ref_pro(i,j,k)*hm_pp/                 &
                   (hm_p(i,j,k)*exner_ref_pro(i,j,k))

      ! Now do "e" and "w" points with update at "p"

      hlm_le(i,j,k) = hm_rhox(i,j,k)*hm_u(i,j,k)/                &
                      (xi1_p(i+1) - xi1_p(i))
      hlm_lw(i,j,k) = hm_rhox(i-1,j,k)*hm_u(i-1,j,k)/            &
                      (xi1_p(i) - xi1_p(i-1))

      hlm_le(i,j,k) = hm_vol(i,j,k)*hlm_le(i,j,k)*rdxi1
      hlm_lw(i,j,k) = hm_vol(i,j,k)*hlm_lw(i,j,k)*rdxi1

      ! Now do "n" and "s" points with update at "p"

      hlm_ln(i,j,k) = hm_rhoy(i,j,k)*hm_v(i,j,k)/                &
                      (xi2_p(j+1) - xi2_p(j))
      hlm_ls(i,j,k) = hm_rhoy(i,j-1,k)*hm_v(i,j-1,k)/            &
                      (xi2_p(j) - xi2_p(j-1))

      hlm_ln(i,j,k) = hm_vol(i,j,k)*hlm_ln(i,j,k)*rdxi2
      hlm_ls(i,j,k) = hm_vol(i,j,k)*hlm_ls(i,j,k)*rdxi2

      hlm_lp(i,j,k)=hlm_lp(i,j,k)-((hlm_ln(i,j,k)+hlm_ls(i,j,k))+      &
                                   (hlm_le(i,j,k)+hlm_lw(i,j,k)))      &
                                  *deta_xi3(i,j,k)

      hlm_le(i,j,k) = hlm_le(i,j,k)*deta_xi3(i+1,j,k)
      hlm_lw(i,j,k) = hlm_lw(i,j,k)*deta_xi3(i-1,j,k)
      hlm_ln(i,j,k) = hlm_ln(i,j,k)*deta_xi3(i,j+1,k)
      hlm_ls(i,j,k) = hlm_ls(i,j,k)*deta_xi3(i,j-1,k)
    END DO
  END DO
END DO
!$OMP END DO NOWAIT


! Now do "u" and "d" points with final update of "p"
! The top and bottom levels are done separately

!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels

  IF ( k == 1 ) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        hlm_lu(i,j,k) = hlm_ek(i,j,k)*hlm_ck(i,j,k)*hm_w(i,j,k)*       &
                        rd_rho(2)

        hlm_lp(i,j,k) = hlm_lp(i,j,k) - hlm_lu(i,j,k)
        hlm_ld(i,j,k) = 0.0

      END DO
    END DO
  ELSE IF ( k == model_levels ) THEN
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        hlm_ld(i,j,k) = -hlm_fk(i,j,k)*hlm_ck(i,j,k-1)*hm_w(i,j,k-1)*  &
                        rd_rho(k)

        hlm_lp(i,j,k) = hlm_lp(i,j,k) - hlm_ld(i,j,k)
        hlm_lu(i,j,k) = 0.0

      END DO
    END DO
  ELSE
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        hlm_lu(i,j,k) = hlm_ek(i,j,k)*hlm_ck(i,j,k)*hm_w(i,j,k)*       &
                        rd_rho(k+1)
        hlm_ld(i,j,k) =-hlm_fk(i,j,k)*hlm_ck(i,j,k-1)*hm_w(i,j,k-1)*   &
                        rd_rho(k)

        hlm_lp(i,j,k) = hlm_lp(i,j,k) - (hlm_lu(i,j,k) + hlm_ld(i,j,k))

      END DO
    END DO
  END IF
  ! Now rescale to make the diagonal unity
  ! NB. Rescaling is applied to RHS in eg_helm_var_rhs.
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      hlm_lp(i,j,k) = 1.0/hlm_lp(i,j,k)
      hlm_lu(i,j,k) = hlm_lu(i,j,k)*hlm_lp(i,j,k)
      hlm_ld(i,j,k) = hlm_ld(i,j,k)*hlm_lp(i,j,k)
      hlm_le(i,j,k) = hlm_le(i,j,k)*hlm_lp(i,j,k)
      hlm_lw(i,j,k) = hlm_lw(i,j,k)*hlm_lp(i,j,k)
      hlm_ln(i,j,k) = hlm_ln(i,j,k)*hlm_lp(i,j,k)
      hlm_ls(i,j,k) = hlm_ls(i,j,k)*hlm_lp(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

! Dirichelt BCs for a LAM

IF (model_type == mt_lam) THEN
  IF (neighbour(pwest)  ==  nodomain) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        HLM_Lp(1,j,k) = 1.0
        HLM_Ln(1,j,k) = 0.0
        HLM_Ls(1,j,k) = 0.0
        HLM_Le(1,j,k) = 0.0
        HLM_Lw(1,j,k) = 0.0
        HLM_Lu(1,j,k) = 0.0
        HLM_Ld(1,j,k) = 0.0
      END DO
    END DO
!$OMP END DO
  END IF
  IF (neighbour(peast)  ==  nodomain) THEN
    i = pdims%i_end - 1     ! one less p point in LAM's
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        HLM_Lp(i,j,k) = 1.0
        HLM_Ln(i,j,k) = 0.0
        HLM_Ls(i,j,k) = 0.0
        HLM_Le(i,j,k) = 0.0
        HLM_Lw(i,j,k) = 0.0
        HLM_Lu(i,j,k) = 0.0
        HLM_Ld(i,j,k) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT

    i = i + 1               ! Force last column to input
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
        HLM_Lp(i,j,k) = 1.0
        HLM_Ln(i,j,k) = 0.0
        HLM_Ls(i,j,k) = 0.0
        HLM_Le(i,j,k) = 0.0
        HLM_Lw(i,j,k) = 0.0
        HLM_Lu(i,j,k) = 0.0
        HLM_Ld(i,j,k) = 0.0
      END DO
    END DO
!$OMP END DO
  END IF
  IF (neighbour(psouth)  ==  nodomain) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO i = pdims%i_start, pdims%i_end
        HLM_Lp(i,1,k) = 1.0
        HLM_Ln(i,1,k) = 0.0
        HLM_Ls(i,1,k) = 0.0
        HLM_Le(i,1,k) = 0.0
        HLM_Lw(i,1,k) = 0.0
        HLM_Lu(i,1,k) = 0.0
        HLM_Ld(i,1,k) = 0.0
      END DO
    END DO
!$OMP END DO
  END IF
  IF (neighbour(pnorth)  ==  nodomain) THEN
    j = pdims%j_end
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO i = pdims%i_start, pdims%i_end
        HLM_Lp(i,j,k) = 1.0
        HLM_Ln(i,j,k) = 0.0
        HLM_Ls(i,j,k) = 0.0
        HLM_Le(i,j,k) = 0.0
        HLM_Lw(i,j,k) = 0.0
        HLM_Lu(i,j,k) = 0.0
        HLM_Ld(i,j,k) = 0.0
      END DO
    END DO
!$OMP END DO
  END IF
END IF

! Build tridiagonal factorization for SOR routine

k = 1
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    Hu_k(i,k,j) =  hlm_lu(i,j,k)
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO k = 2, model_levels
    DO i = pdims%i_start, pdims%i_end
      Hd_k(i,k,j) = 1.0/(1.0 - Hlm_Ld(i,j,k)*Hu_k(i,k-1,j))
      Hu_k(i,k,j) = Hlm_Lu(i,j,k)*Hd_k(i,k,j)
    END DO
  END DO
END DO
!$OMP END DO

! Set the red-black ordered versions of the Helmholtz coefficients
!$OMP DO SCHEDULE(STATIC)
DO k = 1, model_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end, 2
      l = pdims%i_start + (i - pdims%i_start)/2
      thlm_ln(l,j,k) = hlm_ln(i,j,k)
      thlm_ls(l,j,k) = hlm_ls(i,j,k)
      thlm_le(l,j,k) = hlm_le(i,j,k)
      thlm_lw(l,j,k) = hlm_lw(i,j,k)
      thlm_lu(l,j,k) = hlm_lu(i,j,k)
      thlm_ld(l,j,k) = hlm_ld(i,j,k)
      thu_k  (l,k,j) = hu_k  (i,k,j)
      thd_k  (l,k,j) = hd_k  (i,k,j)
    END DO
    DO i = pdims%i_start+1, pdims%i_end, 2
      n = mid_point + (i - 1 - pdims%i_start)/2
      thlm_ln(n,j,k) = hlm_ln(i,j,k)
      thlm_ls(n,j,k) = hlm_ls(i,j,k)
      thlm_le(n,j,k) = hlm_le(i,j,k)
      thlm_lw(n,j,k) = hlm_lw(i,j,k)
      thlm_lu(n,j,k) = hlm_lu(i,j,k)
      thlm_ld(n,j,k) = hlm_ld(i,j,k)
      thu_k  (n,k,j) = hu_k  (i,k,j)
      thd_k  (n,k,j) = hd_k  (i,k,j)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL


! Make adjustments for incremental solver

IF( l_inc_solver ) THEN
   DO k = 1, model_levels
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end
            HM_b(i,j,k)   = HM_b(i,j,k)*HM_w(i,j,k)

             HM_theta(i,j,k) = HM_theta(i,j,k)/deta_xi3_theta(i,j,k)

             c = (Ih+mu_w(i,j,k))-hm_theta(i,j,k)*hm_b(i,j,k)

             hlm_ck(i,j,k) = 1.0/c

             HM_w(i,j,k) = HM_w(i,j,k)*hlm_ck(i,j,k)
          END DO
       END DO
    END DO
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_set_helm_lhs
END MODULE eg_set_helm_lhs_mod
