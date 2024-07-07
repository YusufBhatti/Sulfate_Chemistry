! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE eg_star_mod

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


REAL, ALLOCATABLE ::                            &
      q_star(:,:,:),                            &
    qcl_star(:,:,:),                            &
    qcf_star(:,:,:),                            &
     cf_star(:,:,:),                            &
    cfl_star(:,:,:),                            &
    cff_star(:,:,:),                            &
      m_star(:,:,:),                            &
    mcl_star(:,:,:),                            &
    mcf_star(:,:,:),                            &
   mcf2_star(:,:,:),                            &
  mrain_star(:,:,:),                            &
 mgraup_star(:,:,:)

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_STAR_MOD'

CONTAINS

SUBROUTINE init_star()

USE atm_fields_bounds_mod, ONLY : tdims, tdims_s
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_STAR'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

ALLOCATE(                                                               &
      q_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
    qcl_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
    qcf_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
     cf_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
    cfl_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
    cff_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
      m_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
    mcl_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
    mcf_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
   mcf2_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
  mrain_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end),                                &
 mgraup_star(tdims%i_start:tdims%i_end,                                 &
             tdims%j_start:tdims%j_end,                                 &
             tdims%k_start:tdims%k_end))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)

END SUBROUTINE

END MODULE
