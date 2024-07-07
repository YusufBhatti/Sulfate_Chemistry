! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE update_moisture_fields_mod
USE atm_fields_bounds_mod, ONLY : tdims, tdims_s, tdims_l
USE atm_fields_mod
USE atm_step_local
USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain
USE eg_mix_to_q_mod
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
IMPLICIT NONE

PUBLIC :: update_nd_moisture,update_eg_moisture,update_m_star,update_q_star
PRIVATE


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_MOISTURE_FIELDS_MOD'

CONTAINS

SUBROUTINE update_nd_moisture(l_mr_pc2,l_mr_qtbalcld)

USE eg_q_to_mix_mod

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

LOGICAL :: l_mr_pc2,l_mr_qtbalcld
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_ND_MOISTURE'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (l_mr_pc2 .OR. l_mr_qtbalcld) THEN

  q  (tdims_s%i_start:tdims_s%i_end,                                   &
      tdims_s%j_start:tdims_s%j_end,                                   &
      tdims_s%k_start:tdims_s%k_end)  =                                &
                m_v (:,:,:)
  qcl(tdims_s%i_start:tdims_s%i_end,                                   &
      tdims_s%j_start:tdims_s%j_end,                                   &
      tdims_s%k_start:tdims_s%k_end)  =                                &
                m_cl(:,:,:)
  qcf(tdims_s%i_start:tdims_s%i_end,                                   &
      tdims_s%j_start:tdims_s%j_end,                                   &
      tdims_s%k_start:tdims_s%k_end)  =                                &
                m_cf(:,:,:)

  IF (l_mcr_qcf2  )                                                    &
     qcf2  (tdims_s%i_start:tdims_s%i_end,                             &
            tdims_s%j_start:tdims_s%j_end,:) =                         &
     m_cf2(:,:,:)

  IF (l_mcr_qrain )                                                    &
     qrain (tdims_s%i_start:tdims_s%i_end,                             &
            tdims_s%j_start:tdims_s%j_end, :) =                        &
     m_r(:,:,:)

  IF (l_mcr_qgraup)                                                    &
     qgraup(tdims_s%i_start:tdims_s%i_end,                             &
            tdims_s%j_start:tdims_s%j_end,:) =                         &
     m_gr(:,:,:)

ELSE

  CALL eg_mix_to_q                                                     &
                  (tdims_l,tdims_s,                                    &
                   m_v, m_cl, m_cf,                                    &
                   m_cf2, m_r, m_gr,                                   &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,              &
                   q, qcl, qcf,                                        &
                   qcf2, qrain, qgraup)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

SUBROUTINE update_eg_moisture(l_mr_pc2,l_mr_qtbalcld)

USE eg_q_to_mix_mod

IMPLICIT NONE

LOGICAL :: l_mr_pc2,l_mr_qtbalcld
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_EG_MOISTURE'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (l_mr_pc2 .OR. l_mr_qtbalcld) THEN

  m_v (:,:,:)  =                                                       &
     q  (tdims_s%i_start:tdims_s%i_end,                                &
         tdims_s%j_start:tdims_s%j_end,                                &
         tdims_s%k_start:tdims_s%k_end)

  m_cl(:,:,:)  =                                                       &
     qcl(tdims_s%i_start:tdims_s%i_end,                                &
         tdims_s%j_start:tdims_s%j_end,                                &
         tdims_s%k_start:tdims_s%k_end)
  m_cf(:,:,:)  =                                                       &
     qcf(tdims_s%i_start:tdims_s%i_end,                                &
         tdims_s%j_start:tdims_s%j_end,                                &
         tdims_s%k_start:tdims_s%k_end)

  IF (l_mcr_qcf2  )                                                    &
  m_cf2(:,:,:)  =                                                      &
     qcf2(tdims_s%i_start:tdims_s%i_end,                               &
          tdims_s%j_start:tdims_s%j_end,                               &
          tdims_s%k_start:tdims_s%k_end)

  IF (l_mcr_qrain )                                                    &
  m_r (:,:,:)  =                                                       &
     qrain(tdims_s%i_start:tdims_s%i_end,                              &
           tdims_s%j_start:tdims_s%j_end,                              &
           tdims_s%k_start:tdims_s%k_end)

  IF (l_mcr_qgraup)                                                    &
     m_gr(:,:,:)  =                                                    &
     qgraup(tdims_s%i_start:tdims_s%i_end,                             &
            tdims_s%j_start:tdims_s%j_end,                             &
            tdims_s%k_start:tdims_s%k_end)

  CALL eg_mix_to_q                                                   &
                (tdims_l,tdims_s,                                    &
                 m_v, m_cl, m_cf,                                    &
                 m_cf2, m_r, m_gr,                                   &
                 l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,              &
                 q, qcl, qcf,                                        &
                 qcf2, qrain, qgraup)

ELSE

  CALL eg_q_to_mix                                                   &
                (tdims_l,tdims_s,                                    &
                 q, qcl, qcf,                                        &
                 qcf2, qrain, qgraup,                                &
                 l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,              &
                 m_v, m_cl, m_cf                                     &
                ,m_cf2, m_r, m_gr)

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE


SUBROUTINE update_m_star()
USE eg_q_to_mix_mod
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE eg_star_mod

IMPLICIT NONE

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_M_STAR'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (l_mr_physics) THEN

  m_star  (:,:,:) = q_star  (:,:,:)
  mcl_star(:,:,:) = qcl_star(:,:,:)
  mcf_star(:,:,:) = qcf_star(:,:,:)

  IF (l_mcr_qcf2  )    mcf2_star  (:,:,:) = qcf2_star  (:,:,:)
  IF (l_mcr_qrain )    mrain_star (:,:,:) = qrain_star (:,:,:)
  IF (l_mcr_qgraup)    mgraup_star(:,:,:) = qgraup_star(:,:,:)

  ! NOTE: not using store flavours because they were not used to start with.

ELSE

  ! temporarily convert _star fields from increments to full fields
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)       &
!$OMP& SHARED(tdims,q_star,q,qcl_star,qcl,qcf_star,qcf,l_mcr_qcf2,    &
!$OMP& qcf2_star,qcf2,qrain_star,qrain,qgraup_star,qgraup,            &
!$OMP& l_mcr_qgraup,l_mcr_qrain)
  DO k = tdims%k_start,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end             
        q_star(i,j,k)   = q(i,j,k) + q_star(i,j,k)
        qcl_star(i,j,k) = qcl(i,j,k) + qcl_star(i,j,k)
        qcf_star(i,j,k) = qcf(i,j,k) + qcf_star(i,j,k)
        IF (l_mcr_qcf2) qcf2_star(i,j,k)=qcf2(i,j,k) + qcf2_star(i,j,k)
        IF (l_mcr_qrain) qrain_star(i,j,k)=qrain(i,j,k) +            &
                                           qrain_star(i,j,k)
        IF (l_mcr_qgraup) qgraup_star(i,j,k)=qgraup(i,j,k) +         &
                                             qgraup_star(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  CALL eg_q_to_mix                                                   &
                (tdims,tdims,                                        &
                 q_star, qcl_star, qcf_star,                         &
                 qcf2_star, qrain_star, qgraup_star,                 &
                 l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,              &
                 m_star, mcl_star, mcf_star                          &
                ,mcf2_star, mrain_star, mgraup_star,swap_in=.FALSE.)

  ! now convert back to increments
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)         &
!$OMP& SHARED(tdims,q_star,q,qcl_star,qcl,qcf_star,qcf,l_mcr_qcf2,      &
!$OMP& qcf2_star,qcf2,qrain_star,qrain,qgraup_star,qgraup,              &
!$OMP& l_mcr_qgraup,l_mcr_qrain,m_star,m_v,mcl_star,m_cl,mcf_star,m_cf, &
!$OMP& mrain_star,mcf2_star,mgraup_star,m_cf2,m_r,m_gr)
  DO k = tdims%k_start,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end             
        q_star(i,j,k)   = q_star(i,j,k) - q(i,j,k)
        m_star(i,j,k)   = m_star(i,j,k) - m_v(i,j,k)
        qcl_star(i,j,k) = qcl_star(i,j,k) -qcl(i,j,k)
        mcl_star(i,j,k) = mcl_star(i,j,k) -m_cl(i,j,k)
        qcf_star(i,j,k) = qcf_star(i,j,k) - qcf(i,j,k)
        mcf_star(i,j,k) = mcf_star(i,j,k) - m_cf(i,j,k)
        IF (l_mcr_qcf2) qcf2_star(i,j,k)= qcf2_star(i,j,k) - qcf2(i,j,k)
        IF (l_mcr_qcf2) mcf2_star(i,j,k)= mcf2_star(i,j,k) - m_cf2(i,j,k)
        IF (l_mcr_qrain) qrain_star(i,j,k)= qrain_star(i,j,k) -  &
                                            qrain(i,j,k)
        IF (l_mcr_qrain) mrain_star(i,j,k)= mrain_star(i,j,k) -  &
                                            m_r(i,j,k)
        IF (l_mcr_qgraup) qgraup_star(i,j,k)=qgraup_star(i,j,k) - &
                                             qgraup(i,j,k)
        IF (l_mcr_qgraup) mgraup_star(i,j,k)=mgraup_star(i,j,k) - &
                                             m_gr(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

SUBROUTINE update_q_star()
USE atm_fields_bounds_mod, ONLY : tdims
USE gen_phys_inputs_mod,  ONLY: l_mr_physics
USE eg_star_mod

IMPLICIT NONE

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_Q_STAR'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF (l_mr_physics) THEN
  q_star  (:,:,:) = m_star  (:,:,:)
  qcl_star(:,:,:) = mcl_star(:,:,:)
  qcf_star(:,:,:) = mcf_star(:,:,:)

  IF (l_mcr_qcf2  )    qcf2_star  (:,:,:) = mcf2_star  (:,:,:)
  IF (l_mcr_qrain )    qrain_star (:,:,:) = mrain_star (:,:,:)
  IF (l_mcr_qgraup)    qgraup_star(:,:,:) = mgraup_star(:,:,:)

ELSE
  CALL eg_mix_to_q                                                      &
                  (tdims,tdims,                                         &
                   m_star, mcl_star, mcf_star,                          &
                   mcf2_star, mrain_star, mgraup_star,                  &
                   l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,               &
                   q_star, qcl_star, qcf_star,                          &
                   qcf2_star,qrain_star,qgraup_star                     &
                   )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

END MODULE
