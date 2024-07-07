! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE pressure_grad_mod
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PRESSURE_GRAD_MOD'

CONTAINS
SUBROUTINE pressure_grad(p,Pu, Pv, Pw, Pn)
USE atm_fields_bounds_mod
USE metric_terms_mod
USE eg_helmholtz_mod
USE horiz_grid_mod
USE ref_pro_mod
USE level_heights_mod, ONLY: eta_theta_levels, eta_rho_levels,           &
                             xi3_at_theta=>r_theta_levels
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod, ONLY: swap_field_is_vector
USE helmholtz_const_matrix_mod
USE Field_Types
USE um_parvars,    ONLY: at_extremity
USE UM_ParParams

USE eg_v_at_poles_mod

USE model_domain_mod, ONLY: model_type, mt_global


USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!
! Description: Calculate the three components of the 
!              unapproximated pressure gradient.
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

REAL,           INTENT(IN)  :: p(pdims_s%i_start:pdims_s%i_end,                &
                                 pdims_s%j_start:pdims_s%j_end,                &
                                 pdims_s%k_start:pdims_s%k_end)
REAL,           INTENT(OUT) :: Pu(udims_s%i_start:udims_s%i_end,               &
                                  udims_s%j_start:udims_s%j_end,               &
                                  udims_s%k_start:udims_s%k_end)
REAL,           INTENT(OUT) :: Pv(vdims_s%i_start:vdims_s%i_end,               &
                                  vdims_s%j_start:vdims_s%j_end,               &
                                  vdims_s%k_start:vdims_s%k_end)
REAL,           INTENT(OUT) :: Pw(wdims%i_start:wdims%i_end,                   &
                                  wdims%j_start:wdims%j_end,                   &
                                  wdims%k_start:wdims%k_end)

REAL, OPTIONAL, INTENT(OUT) :: Pn(wdims%i_start:wdims%i_end,                   &
                                  wdims%j_start:wdims%j_end,                   &
                                  wdims%k_start:wdims%k_end)

REAL                        :: Pu_km1(udims%i_start:udims%i_end,               &
                                      udims%j_start:udims%j_end)
REAL                        :: Pv_km1(vdims%i_start:vdims%i_end,               &
                                      vdims%j_start:vdims%j_end)

REAL                        :: r_dxi1, r_dxi2, r_deta, P_k
REAL                        :: u_at_w, v_at_w

INTEGER                     :: i, j, k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRESSURE_GRAD'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! First find exner at xi1 pressure gradient

! Surface terms

Pu_km1(:,:) = 0.0
Pv_km1(:,:) = 0.0

! u-pressure gradient

DO k = pdims%k_start, pdims%k_end-1
   r_deta = 1.0/(eta_theta_levels(k) - eta_theta_levels(k-1))
   IF( k == 1 ) r_deta=0.0

   DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
         r_dxi1 = 1.0/(xi1_p(i+1) - xi1_p(i))

         P_k = (intw_rho2w(k,1)*(intw_p2u(i,1)*p(i,j,k+1)                      &
                                +intw_p2u(i,2)*p(i+1,j,k+1))                   &
              + intw_rho2w(k,2)*(intw_p2u(i,1)*p(i,j,k)                        &
                                +intw_p2u(i,2)*p(i+1,j,k)))                    &
              *(xi3_at_theta(i+1,j,k)-xi3_at_theta(i,j,k))

         Pu(i,j,k) = p(i+1,j,k)*deta_xi3(i+1,j,k) - p(i,j,k)*deta_xi3(i,j,k)
         Pu(i,j,k) = HM_u(i,j,k)*(Pu(i,j,k)                                    &
                        - r_deta*(P_k - Pu_km1(i,j)))*r_dxi1

         Pu_km1(i,j) = P_k
      END DO
   END DO

! v-pressure gradient

   DO j = vdims%j_start, vdims%j_end
      r_dxi2 = 1.0/(xi2_p(j+1) - xi2_p(j))
      DO i = vdims%i_start, vdims%i_end

         P_k = (intw_rho2w(k,1)*(intw_p2v(j,1)*p(i,j,k+1)                      &
                                +intw_p2v(j,2)*p(i,j+1,k+1))                   &
              + intw_rho2w(k,2)*(intw_p2v(j,1)*p(i,j,k)                        &
                                +intw_p2v(j,2)*p(i,j+1,k)))                    &
              *(xi3_at_theta(i,j+1,k)-xi3_at_theta(i,j,k))

         Pv(i,j,k) = p(i,j+1,k)*deta_xi3(i,j+1,k) - p(i,j,k)*deta_xi3(i,j,k)
         Pv(i,j,k) = HM_v(i,j,k)*(Pv(i,j,k)                                    &
                        - r_deta*(P_k - Pv_km1(i,j)))*r_dxi2

         Pv_km1(i,j) = P_k
      END DO
   END DO

! w-points

   r_deta = 1.0/(eta_rho_levels(k+1) - eta_rho_levels(k))

   DO j = wdims%j_start, wdims%j_end
      DO i = wdims%i_start, wdims%i_end

         Pw(i,j,k) = p(i,j,k+1) - p(i,j,k)
         Pw(i,j,k) = HM_w(i,j,k)*Pw(i,j,k)*r_deta

      END DO
   END DO

END DO

! lid

k = udims%k_end
r_deta = 1.0/(eta_theta_levels(k) - eta_theta_levels(k-1))
DO j = udims%j_start, udims%j_end
  DO i = udims%i_start, udims%i_end
    r_dxi1 = 1.0/(xi1_p(i+1) - xi1_p(i))

    Pu(i,j,k) = p(i+1,j,k)*deta_xi3(i+1,j,k) - p(i,j,k)*deta_xi3(i,j,k)
    Pu(i,j,k) = HM_u(i,j,k)*Pu(i,j,k)*r_dxi1

  END DO
END DO

DO j = vdims%j_start, vdims%j_end
  r_dxi2 = 1.0/(xi2_p(j+1) - xi2_p(j))
  DO i = vdims%i_start, vdims%i_end
    Pv(i,j,k) = p(i,j+1,k)*deta_xi3(i,j+1,k) - p(i,j,k)*deta_xi3(i,j,k)
    Pv(i,j,k) = HM_v(i,j,k)*Pv(i,j,k)*r_dxi2
  END DO
END DO

Pw(:,:,0) = 0.0
Pw(:,:,pdims%k_end) = 0.0

! Convert to etatot if required

IF( PRESENT(Pn) ) THEN
   IF ( model_type == mt_global ) THEN
      IF ( at_extremity(PSouth) ) THEN

         CALL eg_v_at_poles(Pu,Pv, 1.0, udims%j_start, vdims%j_start,          &
                                        udims_s,vdims_s)

      END IF

      IF ( at_extremity(PNorth) ) THEN

         CALL eg_v_at_poles(Pu,Pv,-1.0, udims%j_end, vdims%j_end,              &
                                        udims_s,vdims_s)
      END IF
   END IF

   CALL swap_bounds(Pu, udims%i_len, udims%j_len, udims%k_len,                 &
                        udims_s%halo_i, udims_s%halo_j,                        &
                        fld_type_u,swap_field_is_vector)
   CALL swap_bounds(Pv, vdims%i_len, vdims%j_len, vdims%k_len,                 &
                        vdims_s%halo_i, vdims_s%halo_j,                        &
                        fld_type_v,swap_field_is_vector)


   DO k = pdims%k_start, pdims%k_end-1
      DO j = pdims%j_start, pdims%j_end
         DO i = pdims%i_start, pdims%i_end

            u_at_w = intw_rho2w(k,1)*( intw_u2p(i,1)*Pu(i-1,j,k+1)             &
                                      +intw_u2p(i,2)*Pu(i,j,k+1) )             &
                    +intw_rho2w(k,2)*( intw_u2p(i,1)*Pu(i-1,j,k)               &
                                      +intw_u2p(i,2)*Pu(i,j,k) )

            v_at_w = intw_rho2w(k,1)*( intw_v2p(j,1)*Pv(i,j-1,k+1)             &
                                      +intw_v2p(j,2)*Pv(i,j,k+1) )             &
                    +intw_rho2w(k,2)*( intw_v2p(j,1)*Pv(i,j-1,k)               &
                                      +intw_v2p(j,2)*Pv(i,j,k) )

            Pn(i,j,k) = ( Pw(i,j,k)/h3_p_eta(i,j,k)                            &
                         -u_at_w*dxi1_xi3(i,j,k)/h1_p_eta(i,j,k)               &
                         -v_at_w*dxi2_xi3(i,j,k)/h2_p_eta(i,j,k)               &
                        )/deta_xi3_theta(i,j,k)

         END DO
      END DO
   END DO

   Pn(:,:,0) = 0.0
   Pn(:,:,pdims%k_end) = 0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE pressure_grad
END MODULE pressure_grad_mod
