MODULE gravity_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame gravity fields
!
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics Advection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

IMPLICIT NONE

! Horizontal coordinates


REAL, SAVE, ALLOCATABLE ::                                             &
           g_theta(:,:,:) ,                                            &
           g_rho  (:,:,:)

REAL, SAVE, ALLOCATABLE :: g_ref_theta(:)
REAL, SAVE, ALLOCATABLE :: g_ref_rho(:)


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GRAVITY_MOD'

CONTAINS

SUBROUTINE init_gravity()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

USE ereport_mod
USE atm_fields_bounds_mod
USE dynamics_testing_mod, ONLY: problem_number
USE problem_mod,          ONLY: standard

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_GRAVITY'

INTEGER :: ierr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE (    g_theta (tdims_s%i_start:tdims_s%i_end,                  &
                       tdims_s%j_start:tdims_s%j_end,                  &
                       tdims_s%k_start:tdims_s%k_end),                 &
              g_rho   (pdims_s%i_start:pdims_s%i_end,                  &
                       pdims_s%j_start:pdims_s%j_end,                  &
                       pdims_s%k_start:pdims_s%k_end),                 &
                STAT=ierr)

IF (ierr/=0) CALL Ereport(RoutineName, ierr, "Unable to allocate 3D.")

IF (problem_number /= standard) THEN
  ALLOCATE(g_ref_theta(tdims%k_start:tdims%k_end),                     &
           g_ref_rho(pdims%k_start:pdims%k_end),                       &
           STAT=ierr)

  IF (ierr/=0) CALL Ereport(RoutineName,ierr, "Unable to allocate 1D.")
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_gravity
END MODULE gravity_mod
