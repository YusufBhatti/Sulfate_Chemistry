! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE init_psiw_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_PSIW_MOD'

CONTAINS
SUBROUTINE init_psiw( psi_w_surf, psi_w_lid, l_shallow )

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE atm_fields_mod,        ONLY: u,v
USE level_heights_mod,     ONLY: r_theta_levels
USE timestep_mod,          ONLY: timestep
USE horiz_grid_mod
USE metric_terms_mod
USE atm_fields_bounds_mod
USE missing_data_mod,      ONLY: rmdi
USE dynamics_input_mod,    ONLY: ih
USE eg_alpha_mod,          ONLY: alpha_w

USE model_domain_mod, ONLY: model_type, mt_lam

IMPLICIT NONE
!
! Description:
!
!
! Method:
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.


REAL ::                                                                 &
  psi_w_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  psi_w_lid (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

LOGICAL :: l_shallow

! Local
INTEGER :: i, j, k
REAL :: u_at_w, v_at_w


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_PSIW'

! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


IF (       (psi_w_surf(pdims%i_start,pdims%j_end) == rmdi)               &
      .OR. (model_type == mt_lam)) THEN

  ! the reconfiguration cannot do this due to the timestep dependency
  ! (and the shallow flag - but we have ignored that elsewhere)
  ! therefore the reconfiguration fills it with RMDI, which we pick up here
  ! and fill it if required.

  k = 0
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      psi_w_surf(i,j)   = Ih*(                                           &
                          (h3_p_eta(i,j,k)/h1_p_eta(i,j,k))*             &
                          dxi1_xi3(i,j,k)*(intw_u2p(i,1)*u(i-1,j,k+1)    &
                                           +intw_u2p(i,2)*u(i,j,k+1)) +  &
                          (h3_p_eta(i,j,k)/h2_p_eta(i,j,k))*             &
                          dxi2_xi3(i,j,k)*(intw_v2p(j,1)*v(i,j-1,k+1)    &
                                           +intw_v2p(j,2)*v(i,j,k+1)))

      psi_w_surf(i,j) = psi_w_surf(i,j)/(alpha_w*timestep)
      psi_w_lid(i,j) = 0.0
    END DO
  END DO

  IF ( .NOT. l_shallow ) THEN
    k =  pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        psi_w_surf(i,j) = psi_w_surf(i,j)                               &
                          -( (intw_u2p(i,1)*u(i-1,j,1)                  &
                             +intw_u2p(i,2)*u(i,j,1) )**2               &
                            +(intw_v2p(j,1)*v(i,j-1,1)                  &
                             +intw_v2p(j,2)*v(i,j,1) )**2               &
                           )/r_theta_levels(i,j,0)

        psi_w_lid(i,j) = psi_w_lid(i,j)                                 &
                            -( (intw_u2p(i,1)*u(i-1,j,k)                &
                               +intw_u2p(i,2)*u(i,j,k) )**2             &
                              +(intw_v2p(j,1)*v(i,j-1,k)                &
                               +intw_v2p(j,2)*v(i,j,k) )**2             &
                             )/r_theta_levels(i,j,k)

      END DO
    END DO
  END IF
END IF



IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE init_psiw
END MODULE init_psiw_mod
