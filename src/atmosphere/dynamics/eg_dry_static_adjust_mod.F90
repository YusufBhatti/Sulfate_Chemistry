! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE eg_dry_static_adj_mod
IMPLICIT NONE
!
! Description:
!          Enforces static stability on input theta field.
!          Uses very simple sort algorithm resulting in no mixing of
!          temperature just a simple re-arrangement.
!
!
! Method:
!          Adopted from ND version
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_DRY_STATIC_ADJ_MOD'

CONTAINS

SUBROUTINE eg_dry_static_adj(theta,rho,exner,intw_w2rho)

USE atm_fields_bounds_mod, ONLY: pdims,pdims_s,tdims_s
USE planet_constants_mod,  ONLY: p_zero,r,kappa

USE um_parcore,            ONLY: mype
USE ereport_mod,           ONLY: ereport

USE parkind1,              ONLY: jpim, jprb       !DrHook
USE yomhook,               ONLY: lhook, dr_hook   !DrHook
USE umPrintMgr

IMPLICIT NONE

! Subroutine arguments


REAL ::                                                               &
  theta (tdims_s%i_start:tdims_s%i_end,                               &
         tdims_s%j_start:tdims_s%j_end,                               &
         tdims_s%k_start:tdims_s%k_end)

REAL ::                                                               &
  rho  (pdims_s%i_start:pdims_s%i_end,                                &
        pdims_s%j_start:pdims_s%j_end,                                &
        pdims_s%k_start:pdims_s%k_end)

REAL ::                                                               &
  exner  (pdims_s%i_start:pdims_s%i_end,                              &
          pdims_s%j_start:pdims_s%j_end,                              &
          pdims_s%k_start:pdims_s%k_end)

REAL, INTENT (IN) :: intw_w2rho(pdims_s%k_start:pdims_s%k_end,2)


! Local variables

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_DRY_STATIC_ADJ'

INTEGER ::                                                            &
  i,                                                                  &
  j,                                                                  &
  k,                                                                  &
  kk,                                                                 &
  list,                                                               &
                    ! Loop indices
  cnt,                                                                &
  count_old,                                                          &
  iterations,                                                         &
  n_search,                                                           &
  i_search,                                                           &
  n_match

REAL ::                                                               &
  temp
! Local arrays

INTEGER ::                                                            &
  unstable    (pdims%i_end*pdims%j_end*pdims%k_end),                  &
  unstable_old(pdims%i_end*pdims%j_end*pdims%k_end),                  &
  i_prev      (pdims%i_end*pdims%j_end),                              &
  j_prev      (pdims%i_end*pdims%j_end)

INTEGER :: maxit                     ! Maximum of iteration for trying
                                  ! to remove instabillity

INTEGER :: errorstatus

! No External Routines


! 1.0 Start of subroutine code: perform the calculation.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Section 1. Check for inertial instability and remove.
!            Iterative sort procedure is used.
! ----------------------------------------------------------------------

maxit = 2*pdims%k_end

iterations = 0
cnt = 0

DO k = pdims%k_start,pdims%k_end
  DO j = pdims%j_start,pdims%j_end
    DO i = pdims%i_start,pdims%i_end

      IF (theta(i,j,k)  <   theta(i,j,k-1)) THEN
        temp = theta(i,j,k)
        theta(i,j,k) = theta(i,j,k-1)
        theta(i,j,k-1) = temp

        DO kk=MAX(1,k-1),k+1
          !         recompute density to satisfy equation of state
          rho (i,j,kk) = p_zero/(r*                                   &
                          (intw_w2rho(kk,1)*theta(i,j,kk)             &
                         + intw_w2rho(kk,2)*theta(i,j,kk-1) ))        &
                         *(exner(i,j,kk))**((1.0-kappa)/kappa)
        END DO
        cnt = cnt + 1
        unstable(cnt) = k*1e6+j*1000 + i
      END IF
    END DO
  END DO
END DO

IF ( PrintStatus == PrStatus_Diag .AND. mype == 0  .AND. cnt  /=  0) THEN
  CALL umPrint( '================================================',&
      src='eg_dry_static_adjust_mod')
  CALL umPrint( 'static instabillity at (i,j,k,count):',&
      src='eg_dry_static_adjust_mod')
  CALL umPrint( '------------------------------------------------',&
      src='eg_dry_static_adjust_mod')

  DO list = 1, MIN(10,cnt)

    k = unstable(list) * 1e-6
    j = (unstable(list) - k * 1e6) *1e-3
    i = unstable(list) - k * 1e6 - j*1e3

    WRITE(umMessage,'(4I4)') i,j,k,list
    CALL umPrint(umMessage,src='eg_dry_static_adjust_mod')

  END DO

  IF (cnt >  10) THEN
    WRITE(umMessage,'(A,I7)')' further printing suppressed. total count:',cnt
    CALL umPrint(umMessage,src='eg_dry_static_adjust_mod')
  END IF

  CALL umPrint( '================================================',&
      src='eg_dry_static_adjust_mod')

END IF

DO WHILE (cnt  >   0 .AND. iterations  <   maxit)
  count_old = cnt
  DO list = 1, cnt
    unstable_old(list) = unstable(list)
  END DO
  cnt = 0
  iterations = iterations + 1

  n_search = 0
  n_match = 0
  DO list = 1, count_old
    k = unstable_old(list) * 1e-6
    j = (unstable_old(list) - k * 1e6) *1e-3
    i = unstable_old(list) - k * 1e6 - j*1e3
    DO i_search = 1, n_search
      IF (i  ==  i_prev(i_search) .AND.                               &
          j  ==  j_prev(i_search)) THEN
        n_match = 1
      END IF
    END DO
    IF (n_match  ==  0) THEN
      n_search = n_search + 1
      i_prev(n_search) = i
      j_prev(n_search) = j
    END IF
    n_match = 0
  END DO

  DO list = 1, n_search
    i = i_prev(list)
    j = j_prev(list)
    DO k = 1, pdims%k_end
      IF (theta(i,j,k)  <   theta(i,j,k-1)) THEN

        temp = theta(i,j,k)
        theta(i,j,k) = theta(i,j,k-1)
        theta(i,j,k-1) = temp

        DO kk=MAX(1,k-1),k+1
          !   recompute density to satisfy equation of state
          rho (i,j,kk) = p_zero/(r*                                   &
                          (intw_w2rho(kk,1)*theta(i,j,kk)             &
                         + intw_w2rho(kk,2)*theta(i,j,kk-1) ))        &
                         *(exner(i,j,kk))**((1.0-kappa)/kappa)
        END DO

        cnt = cnt + 1
        unstable(cnt) = k*1e6+j*1000 + i
      END IF
    END DO
  END DO

  IF ( PrintStatus == PrStatus_Diag .AND. mype == 0  .AND. cnt  /=  0) THEN
    CALL umPrint( '================================================',&
        src='eg_dry_static_adjust_mod')
    CALL umPrint( 'static instabillity at (i,j,k,count):',&
        src='eg_dry_static_adjust_mod')
    CALL umPrint( '------------------------------------------------',&
        src='eg_dry_static_adjust_mod')

    DO list = 1, MIN(10,cnt)

      k = unstable(list) * 1e-6
      j = (unstable(list) - k * 1e6) *1e-3
      i = unstable(list) - k * 1e6 - j*1e3

      WRITE(umMessage,'(4I4)') i,j,k,list
      CALL umPrint(umMessage,src='eg_dry_static_adjust_mod')

    END DO

    IF (cnt >  10) THEN
      WRITE(umMessage,'(A,I7)')' further printing suppressed. total count:',cnt
      CALL umPrint(umMessage,src='eg_dry_static_adjust_mod')
    END IF

    CALL umPrint( '================================================',&
        src='eg_dry_static_adjust_mod')

  END IF

END DO

IF (iterations  ==  maxit .AND. cnt  /=  0) THEN
  errorstatus = 1
  CALL ereport("eg_dry_static_adj",errorstatus ,                   &
               " unable to remove static instability" )
END IF

IF (iterations > 0) THEN

  errorstatus = -iterations
  CALL ereport("eg_dry_static_adj",errorstatus ,                      &
                  " static instability removed" )

END IF

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_dry_static_adj
END MODULE eg_dry_static_adj_mod
