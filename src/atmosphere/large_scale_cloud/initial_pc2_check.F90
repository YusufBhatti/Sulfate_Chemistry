! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Routine to call pc2_checks with the inputs available in the different
! versions of initial
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Cloud
!
! Code Description:
!   Language: FORTRAN 90
!

MODULE initial_pc2_check_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'INITIAL_PC2_CHECK_MOD'
CONTAINS

SUBROUTINE initial_pc2_check(exner_theta_levels, p_theta_levels, theta,    &
                             cf_bulk, cf_liquid, cf_frozen, q, qcl, qcf,   &
                             l_mr_iau)

USE atm_fields_bounds_mod, ONLY:   &
    pdims, pdims_s,                &
    tdims, tdims_s, tdims_l

USE parkind1, ONLY:                 &
    jprb,                           &
    jpim
USE yomhook, ONLY:                  &
    lhook,                          &
    dr_hook

USE pc2_checks_mod, ONLY: pc2_checks
IMPLICIT NONE

REAL,    INTENT(IN) ::                                     &
  exner_theta_levels    ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(IN) ::                                     &
  p_theta_levels        ( tdims_s%i_start : tdims_s%i_end, &
                          tdims_s%j_start : tdims_s%j_end, &
                          tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT)::                    &
  theta  ( tdims_s%i_start : tdims_s%i_end, &
           tdims_s%j_start : tdims_s%j_end, &
           tdims_s%k_start : tdims_s%k_end )

REAL,    INTENT(INOUT) ::                   &
  q      ( tdims_l%i_start : tdims_l%i_end, &
           tdims_l%j_start : tdims_l%j_end, &
           tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  qcl    ( tdims_l%i_start : tdims_l%i_end, &
           tdims_l%j_start : tdims_l%j_end, &
           tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                   &
  qcf    ( tdims_l%i_start : tdims_l%i_end, &
           tdims_l%j_start : tdims_l%j_end, &
           tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                      &
  cf_bulk ( tdims_l%i_start : tdims_l%i_end,   &
            tdims_l%j_start : tdims_l%j_end,   &
            tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                        &
  cf_liquid ( tdims_l%i_start : tdims_l%i_end,   &
              tdims_l%j_start : tdims_l%j_end,   &
              tdims_l%k_start : tdims_l%k_end )

REAL,    INTENT(INOUT) ::                        &
  cf_frozen ( tdims_l%i_start : tdims_l%i_end,   &
              tdims_l%j_start : tdims_l%j_end,   &
              tdims_l%k_start : tdims_l%k_end )

! are we using mixing ratio code? This is kept consistent with the 
! setting hardwired in inital for the iau:
LOGICAL, INTENT(IN) :: l_mr_iau

! working temperature:
REAL :: t_work(tdims%i_start : tdims%i_end, &
               tdims%j_start : tdims%j_end, &
               1 : tdims%k_end )

! loop counters:
INTEGER :: k,i,j

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INITIAL_PC2_CHECK'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! convert theta to temperature
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i= tdims%i_start, tdims%i_end
      t_work  (i,j,k) = exner_theta_levels(i,j,k) * theta(i,j,k)
    END DO
  END DO
END DO

CALL pc2_checks(p_theta_levels(tdims%i_start:tdims%i_end,             &
                               tdims%j_start:tdims%j_end,             &
                               1:tdims%k_end),                        &
       t_work, cf_bulk(tdims%i_start:tdims%i_end,                     &
                       tdims%j_start:tdims%j_end,                     &
                                   1:tdims%k_end),                    &
       cf_liquid(tdims%i_start:tdims%i_end,                           &
                 tdims%j_start:tdims%j_end,                           &
                             1:tdims%k_end),                          &
       cf_frozen(tdims%i_start:tdims%i_end,                           &
                 tdims%j_start:tdims%j_end,                           &
                             1:tdims%k_end),                          &
       q(tdims%i_start:tdims%i_end,                                   &
         tdims%j_start:tdims%j_end,                                   &
                     1:tdims%k_end),                                  &
       qcl(tdims%i_start:tdims%i_end,                                 &
           tdims%j_start:tdims%j_end,                                 &
                       1:tdims%k_end),                                &
       qcf(tdims%i_start:tdims%i_end,                                 &
           tdims%j_start:tdims%j_end,                                 &
                       1:tdims%k_end),                                &
       l_mr_iau)

! convert back to theta:
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i= tdims%i_start, tdims%i_end
      theta(i,j,k) = t_work  (i,j,k) / exner_theta_levels(i,j,k)
    END DO
  END DO
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE initial_pc2_check

END MODULE initial_pc2_check_mod
