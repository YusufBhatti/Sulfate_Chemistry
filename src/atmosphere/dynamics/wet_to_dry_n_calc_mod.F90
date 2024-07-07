! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE wet_to_dry_n_calc_mod
IMPLICIT NONE
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics

REAL, SAVE, ALLOCATABLE :: wet_to_dry_n(:,:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='WET_TO_DRY_N_CALC_MOD'

CONTAINS

SUBROUTINE destroy_wet_to_dry_n()
IMPLICIT NONE
IF (ALLOCATED(wet_to_dry_n)) DEALLOCATE (wet_to_dry_n)
END SUBROUTINE

SUBROUTINE wet_to_dry_n_calc(q,qcl,qcf,qcf2,qrain,qgraup)

USE atm_fields_bounds_mod, ONLY: tdims_s,pdims,tdims,tdims_l
USE level_heights_mod,     ONLY: r_theta_levels,r_rho_levels
USE ereport_mod,           ONLY: ereport
USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE mphys_inputs_mod,      ONLY: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup

IMPLICIT NONE
REAL, INTENT (IN) ::                                                           &
            q(tdims_l%i_start:tdims_l%i_end,                                   &
              tdims_l%j_start:tdims_l%j_end,                                   &
              tdims_l%k_start:tdims_l%k_end),                                  &
          qcl(tdims_l%i_start:tdims_l%i_end,                                   &
              tdims_l%j_start:tdims_l%j_end,                                   &
              tdims_l%k_start:tdims_l%k_end),                                  &
          qcf(tdims_l%i_start:tdims_l%i_end,                                   &
              tdims_l%j_start:tdims_l%j_end,                                   &
              tdims_l%k_start:tdims_l%k_end)
REAL, INTENT (IN), TARGET ::                                                   &
         qcf2(tdims_l%i_start:tdims_l%i_end,                                   &
              tdims_l%j_start:tdims_l%j_end,                                   &
              tdims_l%k_start:tdims_l%k_end),                                  &
        qrain(tdims_l%i_start:tdims_l%i_end,                                   &
              tdims_l%j_start:tdims_l%j_end,                                   &
              tdims_l%k_start:tdims_l%k_end),                                  &
       qgraup(tdims_l%i_start:tdims_l%i_end,                                   &
              tdims_l%j_start:tdims_l%j_end,                                   &
              tdims_l%k_start:tdims_l%k_end)

REAL ::                                                                        &
         sumq(tdims_l%i_start:tdims_l%i_end,                                   &
              tdims_l%j_start:tdims_l%j_end,                                   &
              tdims_l%k_start:tdims_l%k_end)


REAL :: weight1,weight2,weight3

INTEGER :: i,j,k,ierr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='WET_TO_DRY_N_CALC'


INTEGER :: used_ptr

TYPE pp
  REAL, POINTER :: ptr(:,:,:)
END TYPE pp

TYPE(pp), ALLOCATABLE :: parray(:)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE(parray(3))

used_ptr    = 0
sumq(:,:,:) = 0.0

IF (l_mcr_qcf2) THEN
  used_ptr = used_ptr+1
  parray(used_ptr)%ptr => qcf2
END IF

IF (l_mcr_qrain) THEN
  used_ptr = used_ptr+1
  parray(used_ptr)%ptr => qrain
END IF

IF (l_mcr_qgraup) THEN
  used_ptr = used_ptr+1
  parray(used_ptr)%ptr => qgraup
END IF

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC) PRIVATE(i,j,k)   &
!$OMP& SHARED(tdims_s,tdims,sumq,q,qcl,qcf,parray,used_ptr)
  DO k = tdims_s%k_start, MIN(tdims%k_end+1,tdims_s%k_end)
    DO j =  tdims%j_start,tdims%j_end
      DO i =  tdims%i_start,tdims%i_end
        sumq(i,j,k) = q(i,j,k)  + qcl(i,j,k) + qcf(i,j,k) 
        IF (used_ptr >= 1) sumq(i,j,k)=sumq(i,j,k)+   &
                                       parray(1)%ptr(i,j,k) 
        IF (used_ptr >= 2) sumq(i,j,k)=sumq(i,j,k)+   &
                                       parray(2)%ptr(i,j,k) 
        IF (used_ptr >= 3) sumq(i,j,k)=sumq(i,j,k)+   &
                                       parray(3)%ptr(i,j,k) 
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO


parray(1)%ptr => NULL()
parray(2)%ptr => NULL()
parray(3)%ptr => NULL()

DEALLOCATE(parray)


IF (.NOT. ALLOCATED(wet_to_dry_n)) THEN
  ALLOCATE (wet_to_dry_n(tdims_s%i_start:tdims_s%i_end,                        &
                         tdims_s%j_start:tdims_s%j_end,                        &
                         tdims_s%k_start:tdims_s%k_end),stat=ierr)

  IF (ierr /= 0) THEN
    CALL ereport('wet_to_dry_n_calc_mod',ierr,'failure when allocating')
  END IF
END IF

! Section to calculate wet to dry conversion array.
! Note the code here is only appropriate for moisture variables q, qcl
! and qcf. If qrain, qcf2 etc are present the code needs altering.
!

DO j = pdims%j_start,pdims%j_end
  DO i = pdims%i_start,pdims%i_end

    wet_to_dry_n(i,j,tdims_s%k_start:1)= 1.0 - sumq(i,j,1)

  END DO
END DO

IF (tdims%k_end <  tdims_s%k_end) THEN
  k=tdims%k_end+1
  DO j =  tdims%j_start,tdims%j_end
    DO i =  tdims%i_start,tdims%i_end

      weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
      weight2 = r_rho_levels(i,j,k)   - r_theta_levels(i,j,k-1)
      weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)

      wet_to_dry_n(i,j,k)= ( weight2 + weight1*(1.0-sumq(i,j,k-1)))/weight3

    END DO
  END DO
END IF

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)     &
!$OMP& PRIVATE(i,j,k,weight1,weight2,weight3)        &
!$OMP& SHARED(tdims,pdims,r_theta_levels,r_rho_levels,sumq,wet_to_dry_n)
DO k = 2, tdims%k_end ! this is tdims, intentionally.
  DO j = pdims%j_start,pdims%j_end
    DO i =  pdims%i_start,pdims%i_end

      weight1 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
      weight2 = r_rho_levels(i,j,k)   - r_theta_levels(i,j,k-1)
      weight3 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)

      wet_to_dry_n(i,j,k)= ( weight2 * (1.0 - sumq(i,j,k))               &
                            +weight1 * (1.0 - sumq(i,j,k-1)) ) / weight3

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE
END MODULE
