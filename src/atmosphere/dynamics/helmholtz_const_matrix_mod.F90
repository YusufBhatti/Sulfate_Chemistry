MODULE helmholtz_const_matrix_mod
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Description: Contains the ENDGame departure points
!
!
! Method:
!
! Documentation: ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

USE um_types, ONLY: real64,real32
IMPLICIT NONE

! Horizontal coordinates
#if defined(C_DP_HLM)
INTEGER, PARAMETER :: real_eg_hlm_kind=real64
#else
INTEGER, PARAMETER :: real_eg_hlm_kind=real32
#endif

REAL (KIND=real_eg_hlm_kind) , SAVE, ALLOCATABLE ::                    &
         Hlm_Lp(:,:,:),                                                &
         Hlm_Ln(:,:,:),                                                &
         Hlm_Ls(:,:,:),                                                &
         Hlm_Le(:,:,:),                                                &
         Hlm_Lw(:,:,:),                                                &
         Hlm_Lu(:,:,:),                                                &
         Hlm_Ld(:,:,:),                                                &
         Hlm_Ek(:,:,:),                                                &
         Hlm_Fk(:,:,:),                                                &
         Hlm_Ck(:,:,:),                                                &
         Hu_k(:,:,:), Hd_k(:,:,:)


! Red-black ordered versions
REAL (KIND=real_eg_hlm_kind) , SAVE, ALLOCATABLE ::                    &
         tHlm_Ln(:,:,:),                                               &
         tHlm_Ls(:,:,:),                                               &
         tHlm_Le(:,:,:),                                               &
         tHlm_Lw(:,:,:),                                               &
         tHlm_Lu(:,:,:),                                               &
         tHlm_Ld(:,:,:),                                               &
         tHu_k(:,:,:), tHd_k(:,:,:)


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='HELMHOLTZ_CONST_MATRIX_MOD'

CONTAINS

SUBROUTINE init_helmholtz_const_matrix()

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod
USE atm_fields_bounds_mod

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_HELMHOLTZ_CONST_MATRIX'

INTEGER :: ierr

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE (   Hlm_Lp   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Ln   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Ls   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Le   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Lw   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Lu   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Ld   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Ek   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Fk   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
             Hlm_Ck   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       0:pdims%k_end),                             &
             Hu_k     (pdims%i_start:pdims%i_end,                  &
                       pdims%k_start:pdims%k_end,                  &
                       pdims%j_start:pdims%j_end),                 &
             Hd_k     (pdims%i_start:pdims%i_end,                  &
                       pdims%k_start:pdims%k_end,                  &
                       pdims%j_start:pdims%j_end),                 &
                STAT=ierr)

IF (ierr/=0) CALL Ereport("init_helmholtz_const_matrix",ierr,      &
                            "Unable to allocate.")

ALLOCATE (  tHlm_Ln   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
            tHlm_Ls   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
            tHlm_Le   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
            tHlm_Lw   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
            tHlm_Lu   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
            tHlm_Ld   (pdims%i_start:pdims%i_end,                  &
                       pdims%j_start:pdims%j_end,                  &
                       pdims%k_start:pdims%k_end),                 &
            tHu_k     (pdims%i_start:pdims%i_end,                  &
                       pdims%k_start:pdims%k_end,                  &
                       pdims%j_start:pdims%j_end),                 &
            tHd_k     (pdims%i_start:pdims%i_end,                  &
                       pdims%k_start:pdims%k_end,                  &
                       pdims%j_start:pdims%j_end),                 &
                STAT=ierr)

IF (ierr/=0) CALL Ereport("init_helmholtz_const_matrix",ierr,      &
                            "Unable to allocate.")

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE init_helmholtz_const_matrix
END MODULE helmholtz_const_matrix_mod
