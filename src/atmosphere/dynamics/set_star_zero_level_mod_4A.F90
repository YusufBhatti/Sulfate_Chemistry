! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! set the zeroth level for physics _star variables
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics

MODULE set_star_zero_level_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_STAR_ZERO_LEVEL_MOD'

CONTAINS
SUBROUTINE set_star_zero_level(                                       &
             theta_star,                                              &
             q_star,                                                  &
             qcl_star,                                                &
             qcf_star,                                                &
             cf_star,                                                 &
             cfl_star,                                                &
             cff_star,                                                &
             qcf2_star,                                               &
             qrain_star,                                              &
             qgraup_star,                                             &
             L_mcr_qgraup,                                            &
             L_mcr_qrain,                                             &
             L_mcr_qcf2)

USE atm_fields_bounds_mod, ONLY : tdims, tdims_s
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

REAL, INTENT (INOUT) :: theta_star                                    &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: q_star                                        &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qcl_star                                      &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qcf_star                                      &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qcf2_star                                     &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qrain_star                                    &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: qgraup_star                                   &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: cf_star                                       &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: cfl_star                                      &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)
REAL, INTENT (INOUT) :: cff_star                                      &
                                 (tdims%i_start:tdims%i_end,          &
                                  tdims%j_start:tdims%j_end,          &
                                  tdims%k_start:tdims%k_end)

LOGICAL, INTENT (IN) ::  l_mcr_qgraup
LOGICAL, INTENT (IN) ::  l_mcr_qrain
LOGICAL, INTENT (IN) ::  l_mcr_qcf2

INTEGER :: i,j
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_STAR_ZERO_LEVEL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j)                                &
!$OMP          SHARED(tdims,theta_star,q_star,qcl_star,qcf_star,cf_star, &
!$OMP                 cfl_star,cff_star,qcf2_star,qrain_star,qgraup_star,&
!$OMP                 l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup)
!$OMP DO SCHEDULE(STATIC)
DO j=tdims%j_start, tdims%j_end
  DO i=tdims%i_start, tdims%i_end
    theta_star(i,j,0) = theta_star(i,j,1)
    q_star(i,j,0)     = q_star(i,j,1)
    qcl_star(i,j,0)   = qcl_star(i,j,1)
    qcf_star(i,j,0)   = qcf_star(i,j,1)
    cf_star(i,j,0)    = cf_star(i,j,1)
    cfl_star(i,j,0)   = cfl_star(i,j,1)
    cff_star(i,j,0)   = cff_star(i,j,1)
  END DO
END DO
!$OMP END DO NOWAIT


IF (l_mcr_qcf2) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      qcf2_star  (i,j,0) = qcf2_star(i,j,1)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (l_mcr_qrain) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      qrain_star(i,j,0) = qrain_star(i,j,1)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF 

IF (l_mcr_qgraup) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j=tdims%j_start, tdims%j_end
    DO i=tdims%i_start, tdims%i_end
      qgraup_star(i,j,0) = qgraup_star(i,j,1)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE set_star_zero_level
END MODULE set_star_zero_level_mod
