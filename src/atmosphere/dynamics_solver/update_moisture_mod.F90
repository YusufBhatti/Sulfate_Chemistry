! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE update_moisture_mod
IMPLICIT NONE

! Description:         updates moisture to new time level
!
! Method:
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics Solver
!
! Code description:
! Language: Fortran 95.
! This code is written to UMDP3 standards.
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_MOISTURE_MOD'

CONTAINS
SUBROUTINE update_moisture(m_v_np1, m_cl_np1, m_cf_np1, m_r_np1,&
                        m_gr_np1, m_cf2_np1,                    &
                        R_m_v_d, R_m_cl_d, R_m_cf_d,            &
                        R_m_r_d, R_m_gr_d, R_m_cf2_d,           &
                        R_m_v_a, R_m_cl_a, R_m_cf_a,            &
                        R_m_r_a, R_m_gr_a, R_m_cf2_a)


USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim
USE atm_fields_bounds_mod
USE mpp_conf_mod,        ONLY: swap_field_is_scalar
USE UM_ParParams
USE mphys_inputs_mod,    ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

IMPLICIT NONE


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


REAL, INTENT(IN) ::   R_m_v_a(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::  R_m_cl_a(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::  R_m_cf_a(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::   R_m_r_a(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) ::  R_m_gr_a(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(IN) :: R_m_cf2_a(tdims%i_start:tdims%i_end,        &
                              tdims%j_start:tdims%j_end,        &
                              tdims%k_start:tdims%k_end)

REAL, INTENT(OUT) ::    m_v_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) ::   m_cl_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) ::   m_cf_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) ::    m_r_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) ::   m_gr_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) ::  m_cf2_np1(tdims_s%i_start:tdims_s%i_end,  &
                                tdims_s%j_start:tdims_s%j_end,  &
                                tdims_s%k_start:tdims_s%k_end)

INTEGER :: i, j, k
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_MOISTURE'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP&    PRIVATE(i,j,k)                                               &
!$OMP&    SHARED(tdims, l_mcr_qrain,l_mcr_qgraup,l_mcr_qcf2,           &
!$OMP&           m_v_np1, m_cl_np1, m_cf_np1,                          &
!$OMP&           m_r_np1, m_gr_np1, m_cf2_np1,                         &
!$OMP&           R_m_v_d, R_m_cl_d, R_m_cf_d,                          &
!$OMP&           R_m_v_a, R_m_cl_a, R_m_cf_a,                          &
!$OMP&           R_m_r_d, R_m_gr_d, R_m_cf2_d,                         &
!$OMP&           R_m_r_a, R_m_gr_a, R_m_cf2_a)
DO k = tdims%k_start,tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      m_v_np1(i,j,k)   = R_m_v_d(i,j,k)   + R_m_v_a(i,j,k)
      m_cl_np1(i,j,k)  = R_m_cl_d(i,j,k)  + R_m_cl_a(i,j,k)
      m_cf_np1(i,j,k)  = R_m_cf_d(i,j,k)  + R_m_cf_a(i,j,k)
    END DO
  END DO

  IF ( l_mcr_qrain ) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        m_r_np1(i,j,k)   = R_m_r_d(i,j,k)   + R_m_r_a(i,j,k)
      END DO
    END DO
  END IF

  IF ( l_mcr_qgraup ) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        m_gr_np1(i,j,k)  = R_m_gr_d(i,j,k)  + R_m_gr_a(i,j,k)
      END DO
    END DO
  END IF

  IF ( l_mcr_qcf2 ) THEN
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        m_cf2_np1(i,j,k) = R_m_cf2_d(i,j,k) + R_m_cf2_a(i,j,k)
      END DO
    END DO
  END IF

END DO
!$OMP END PARALLEL DO


! DEPENDS ON: swap_bounds
CALL swap_bounds(m_v_np1,tdims%i_len,                            &
                         tdims%j_len,                            &
                         tdims%k_len,                            &
              tdims_s%halo_i, tdims_s%halo_j,fld_type_p, swap_field_is_scalar)

! DEPENDS ON: swap_bounds
CALL swap_bounds(m_cl_np1,tdims%i_len,                           &
                          tdims%j_len,                           &
                          tdims%k_len,                           &
              tdims_s%halo_i, tdims_s%halo_j,fld_type_p, swap_field_is_scalar)

! DEPENDS ON: swap_bounds
CALL swap_bounds(m_cf_np1,tdims%i_len,                           &
                          tdims%j_len,                           &
                          tdims%k_len,                           &
              tdims_s%halo_i, tdims_s%halo_j,fld_type_p, swap_field_is_scalar)

! DEPENDS ON: swap_bounds
IF ( l_mcr_qrain )                                               &
  CALL swap_bounds(m_r_np1,tdims%i_len,                          &
                           tdims%j_len,                          &
                           tdims%k_len,                          &
              tdims_s%halo_i, tdims_s%halo_j,fld_type_p, swap_field_is_scalar)

! DEPENDS ON: swap_bounds
IF ( l_mcr_qgraup )                                              &
  CALL swap_bounds(m_gr_np1,tdims%i_len,                         &
                            tdims%j_len,                         &
                            tdims%k_len,                         &
              tdims_s%halo_i, tdims_s%halo_j,fld_type_p, swap_field_is_scalar)

! DEPENDS ON: swap_bounds
IF ( l_mcr_qcf2 )                                                &
  CALL swap_bounds(m_cf2_np1,tdims%i_len,                        &
                             tdims%j_len,                        &
                             tdims%k_len,                        &
              tdims_s%halo_i, tdims_s%halo_j,fld_type_p, swap_field_is_scalar)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE update_moisture

END MODULE update_moisture_mod
