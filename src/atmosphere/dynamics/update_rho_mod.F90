! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE update_rho_mod
IMPLICIT NONE

INTEGER, PARAMETER :: tl_n     = 1
INTEGER, PARAMETER :: tl_np1   = 2
INTEGER, PARAMETER :: tl_np1_2 = 3

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_RHO_MOD'

CONTAINS

SUBROUTINE update_dryrho(rho_in)

USE eg_q_to_mix_mod
USE atm_fields_bounds_mod
USE atm_fields_mod
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE field_types
USE level_heights_mod
USE atm_step_local
USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE horiz_grid_mod, ONLY: intw_w2rho
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
!
! Description:
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Solver
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

REAL :: rho_in      (pdims_s%i_start:pdims_s%i_end,                        &
                  pdims_s%j_start:pdims_s%j_end,                        &
                  pdims_s%k_start:pdims_s%k_end)

REAL :: mixing_r

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_DRYRHO'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

CALL eg_q_to_mix                                                &
          (tdims_l,tdims_s,                                     &
           q, qcl, qcf,                                         &
           qcf2, qrain, qgraup,                                 &
           l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
           m_v_np1, m_cl_np1, m_cf_np1                          &
          ,m_cf2_np1, m_r_np1, m_gr_np1)

DO k=pdims_s%k_start, pdims_s%k_end
  DO j=pdims_s%j_start, pdims_s%j_end
    DO i=pdims_s%i_start, pdims_s%i_end

      mixing_r=(intw_w2rho(k,1)*(m_v_np1  (i,j,k)      +        &
                                 m_r_np1  (i,j,k)      +        &
                                 m_gr_np1 (i,j,k)      +        &
                                 m_cl_np1 (i,j,k)      +        &
                                 m_cf_np1 (i,j,k)      +        &
                                 m_cf2_np1(i,j,k))     +        &
                intw_w2rho(k,2)*(m_v_np1  (i,j,k-1)    +        &
                                 m_r_np1  (i,j,k-1)    +        &
                                 m_gr_np1 (i,j,k-1)    +        &
                                 m_cl_np1 (i,j,k-1)    +        &
                                 m_cf_np1 (i,j,k-1)    +        &
                                 m_cf2_np1(i,j,k-1)))

      !     convert to dry rho
      dryrho(i,j,k) = rho_in(i,j,k)/(1.0+mixing_r)

      !     dividing by r**2
      dryrho(i,j,k) = dryrho(i,j,k)/( r_rho_levels(i,j,k)       &
                                       *r_rho_levels(i,j,k))
    END DO
  END DO
END DO

CALL swap_bounds(dryrho,                                       &
         pdims_s%i_len - 2*pdims_s%halo_i,                     &
         pdims_s%j_len - 2*pdims_s%halo_j,                     &
         pdims_s%k_len,                                        &
         pdims_s%halo_i, pdims_s%halo_j, fld_type_p, swap_field_is_scalar)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE


SUBROUTINE update_dryrho_mix_tl_n(rho_in)

USE eg_q_to_mix_mod
USE atm_fields_bounds_mod
USE atm_fields_mod
USE field_types
USE level_heights_mod
USE atm_step_local
USE mphys_inputs_mod,  ONLY:  &
         l_mcr_qcf2,              &
         l_mcr_qrain,             &
         l_mcr_qgraup
USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar
USE model_domain_mod, ONLY: model_type, mt_lam
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
!
! Description:
!
!
!
! Method:
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

REAL :: rho_in      (pdims_s%i_start:pdims_s%i_end,                              &
                  pdims_s%j_start:pdims_s%j_end,                              &
                  pdims_s%k_start:pdims_s%k_end)

REAL :: mixing_r, intw_w2rho(2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_DRYRHO_MIX_TL_N'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP             PRIVATE(k,j,i,mixing_r,intw_w2rho) SHARED(pdims,     &
!$OMP                     eta_rho_levels, eta_theta_levels, m_v, m_r,  &
!$OMP                     m_gr, m_cl, m_cf, m_cf2, dryrho, rho_in,     &
!$OMP                     r_rho_levels)
DO k=pdims%k_start, pdims%k_end

  intw_w2rho(1) = ( eta_rho_levels(k)  -eta_theta_levels(k-1) ) / &
                      ( eta_theta_levels(k)-eta_theta_levels(k-1) )
  intw_w2rho(2) = 1.0 - intw_w2rho(1)

  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end


      mixing_r=(intw_w2rho(1)*(m_v  (i,j,k)        +                  &
                                 m_r  (i,j,k)      +                  &
                                 m_gr (i,j,k)      +                  &
                                 m_cl (i,j,k)      +                  &
                                 m_cf (i,j,k)      +                  &
                                 m_cf2(i,j,k))     +                  &
                intw_w2rho(2)*(m_v  (i,j,k-1)      +                  &
                                 m_r  (i,j,k-1)    +                  &
                                 m_gr (i,j,k-1)    +                  &
                                 m_cl (i,j,k-1)    +                  &
                                 m_cf (i,j,k-1)    +                  &
                                 m_cf2(i,j,k-1)))

      !     convert to dry rho
      dryrho(i,j,k) = rho_in(i,j,k)/(1.0+mixing_r)

      !     dividing by r**2
      dryrho(i,j,k) = dryrho(i,j,k)/( r_rho_levels(i,j,k)             &
                                     *r_rho_levels(i,j,k))
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

CALL swap_bounds(dryrho,                                       &
         pdims_s%i_len - 2*pdims_s%halo_i,                     &
         pdims_s%j_len - 2*pdims_s%halo_j,                     &
         pdims_s%k_len,                                        &
         pdims_s%halo_i, pdims_s%halo_j, fld_type_p,           &
         swap_field_is_scalar)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE


SUBROUTINE update_wetrho_r_sq(time_level)

USE atm_fields_bounds_mod
USE atm_fields_mod
USE level_heights_mod
USE atm_step_local
USE horiz_grid_mod, ONLY: intw_w2rho
USE ereport_mod
USE model_domain_mod, ONLY: model_type,mt_lam
USE field_types
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER :: time_level,errorstatus

REAL :: mixing_r
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_WETRHO_R_SQ'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

SELECT CASE(time_level)

CASE (tl_n)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,mixing_r) &
!$OMP& SHARED(pdims_s,intw_w2rho,m_v,m_r,m_gr,m_cl,m_cf, &
!$OMP& m_cf2,wetrho_r_sq_n,dryrho,r_rho_levels)
  DO k=pdims_s%k_start, pdims_s%k_end
    DO j=pdims_s%j_start, pdims_s%j_end
      DO i=pdims_s%i_start, pdims_s%i_end

        mixing_r = (intw_w2rho(k,1)*(m_v  (i,j,k)    +              &
                               m_r  (i,j,k)      +                  &
                               m_gr (i,j,k)      +                  &
                               m_cl (i,j,k)      +                  &
                               m_cf (i,j,k)      +                  &
                               m_cf2(i,j,k))     +                  &
              intw_w2rho(k,2)*(m_v  (i,j,k-1)    +                  &
                               m_r  (i,j,k-1)    +                  &
                               m_gr (i,j,k-1)    +                  &
                               m_cl (i,j,k-1)    +                  &
                               m_cf (i,j,k-1)    +                  &
                               m_cf2(i,j,k-1)))

        !           convert to wet rho
        wetrho_r_sq_n(i,j,k) = dryrho(i,j,k)*(1.0+mixing_r)
        !           multiply by r**2
        wetrho_r_sq_n(i,j,k) = wetrho_r_sq_n(i,j,k)*                &
                               ( r_rho_levels(i,j,k) &
                                *r_rho_levels(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

CASE (tl_np1)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,mixing_r) &
!$OMP& SHARED(pdims_s,intw_w2rho,m_v_np1,m_r_np1,m_gr_np1,m_cl_np1,m_cf_np1, &
!$OMP& m_cf2_np1,wetrho_r_sq_np1,dryrho_np1,r_rho_levels)
  DO k=pdims_s%k_start, pdims_s%k_end
    DO j=pdims_s%j_start, pdims_s%j_end
      DO i=pdims_s%i_start, pdims_s%i_end

        mixing_r = (intw_w2rho(k,1)*(m_v_np1  (i,j,k)    +              &
                               m_r_np1  (i,j,k)      +                  &
                               m_gr_np1 (i,j,k)      +                  &
                               m_cl_np1 (i,j,k)      +                  &
                               m_cf_np1 (i,j,k)      +                  &
                               m_cf2_np1(i,j,k))     +                  &
              intw_w2rho(k,2)*(m_v_np1  (i,j,k-1)    +                  &
                               m_r_np1  (i,j,k-1)    +                  &
                               m_gr_np1 (i,j,k-1)    +                  &
                               m_cl_np1 (i,j,k-1)    +                  &
                               m_cf_np1(i,j,k-1)     +                  &
                               m_cf2_np1(i,j,k-1)))

        !           convert to wet rho
        wetrho_r_sq_np1(i,j,k) = dryrho_np1(i,j,k)*(1.0+mixing_r)

        !           multiply by r**2
        wetrho_r_sq_np1(i,j,k) = wetrho_r_sq_np1(i,j,k)*&
                                  ( r_rho_levels(i,j,k) &
                                   *r_rho_levels(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

IF (model_type == mt_lam) THEN
!DEPENDS ON: fill_external_halos
  CALL fill_external_halos(wetrho_r_sq_np1,                        &
             pdims_s%i_len - 2*pdims_s%halo_i,                     &
             pdims_s%j_len - 2*pdims_s%halo_j,                     &
             pdims_s%k_len,                                        &
             pdims_s%halo_i, pdims_s%halo_j)
END IF


CASE (tl_np1_2)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k,mixing_r) &
!$OMP& SHARED(pdims_s,intw_w2rho,m_v,m_r,m_gr,m_cl,m_cf, &
!$OMP& m_cf2,wetrho_r_sq_n,dryrho,r_rho_levels)
  DO k=pdims_s%k_start, pdims_s%k_end
    DO j=pdims_s%j_start, pdims_s%j_end
      DO i=pdims_s%i_start, pdims_s%i_end

        mixing_r = (intw_w2rho(k,1)*(m_v  (i,j,k)    +              &
                               m_r  (i,j,k)      +                  &
                               m_gr (i,j,k)      +                  &
                               m_cl (i,j,k)      +                  &
                               m_cf (i,j,k)      +                  &
                               m_cf2(i,j,k))     +                  &
              intw_w2rho(k,2)*(m_v  (i,j,k-1)    +                  &
                               m_r  (i,j,k-1)    +                  &
                               m_gr (i,j,k-1)    +                  &
                               m_cl (i,j,k-1)    +                  &
                               m_cf(i,j,k-1)    +                   &
                               m_cf2(i,j,k-1)))

        !           convert to wet rho
        wetrho_r_sq_n(i,j,k) = dryrho(i,j,k)*(1.0+mixing_r)
                ! it should use dryrho_np1, to be consistent, but
                ! then throughout atm_step, once the outer loop is
                ! done, the np1 variable is no longer used and only
                ! dryrho etc updated. We have to go through carefully
                ! and change this at some point!

                ! also  it should use  wetrho_r_sq_np1, however, this is not
                ! implemented this way due to some (arguable) systems
                ! constraints (see comments in atm_step_4a).
        !           multiply by r**2
        wetrho_r_sq_n(i,j,k) = wetrho_r_sq_n(i,j,k)*                &
                                 ( r_rho_levels(i,j,k)              &
                                  *r_rho_levels(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

IF (model_type == mt_lam) THEN
!DEPENDS ON: fill_external_halos
  CALL fill_external_halos(wetrho_r_sq_n,                          &
             pdims_s%i_len - 2*pdims_s%halo_i,                     &
             pdims_s%j_len - 2*pdims_s%halo_j,                     &
             pdims_s%k_len,                                        &
             pdims_s%halo_i, pdims_s%halo_j)
END IF


CASE DEFAULT

  ErrorStatus = 1

  CALL ereport('update_dryrho',ErrorStatus,'unknown time level selected')


END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

END MODULE
