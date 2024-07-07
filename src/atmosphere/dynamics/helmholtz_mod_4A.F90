! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Subroutine Interface:
MODULE eg_helmholtz_mod

IMPLICIT NONE

REAL, SAVE, ALLOCATABLE ::  ec_vol(:,:,:)
REAL, SAVE, ALLOCATABLE ::  ec_area(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_u(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_v(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_w(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_theta(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_etadot(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_vol(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_rhox(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_rhoy(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_rhoz(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_p(:,:,:)
REAL, SAVE, ALLOCATABLE ::  hm_b(:,:,:)
REAL, SAVE              ::  hm_pp

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_HELMHOLTZ_MOD'

CONTAINS


SUBROUTINE eg_init_helmholtz(row_length,rows,n_rows,model_levels,offx,offy)
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE ereport_mod, ONLY: ereport

IMPLICIT NONE
!
! Description: allocates constants required for the
!              solution to the  Helmholtz problem.
!
!
! Method: Appendix F, ENDGame formulation version 1.01
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Control
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


INTEGER :: row_length,rows,n_rows,model_levels,offx,offy

INTEGER :: alloc_stat

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_INIT_HELMHOLTZ'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


ALLOCATE ( ec_vol(row_length,rows,model_levels), STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  ec_area(row_length,rows,model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_u(-offx:row_length-1+offx,1-offy:rows+offy,            &
                  model_levels), STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_v(1-offx:row_length+offx,-offy:n_rows-1+offy,          &
                 model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_w(row_length,rows,0:model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_b(row_length,rows, 0:model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")


ALLOCATE (  hm_theta(row_length,rows,0:model_levels),STAT=alloc_stat )
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_etadot(row_length,rows,0:model_levels), STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_vol(row_length,rows, model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_rhox(-offx:row_length-1+offx,1-offy:rows+offy,         &
                         model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_rhoy(1-offx:row_length+offx,-offy:n_rows-1+offy,       &
          model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_rhoz(row_length,rows,0:model_levels), STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

ALLOCATE (  hm_p(row_length,rows,0:model_levels),STAT=alloc_stat)
IF (alloc_stat/=0) CALL Ereport("ATM_STEP_ALLOC",                   &
                       alloc_stat, "Unable to allocate.")

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_init_helmholtz

SUBROUTINE eg_destroy_helmholtz()
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE
!
! Description: frees constants required for the
!              solution to the  Helmholtz problem.
!
!
! Method: Appendix F, ENDGame formulation version 1.01
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Subroutine arguments


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_DESTROY_HELMHOLTZ'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DEALLOCATE ( hm_p )
DEALLOCATE ( hm_rhoz )
DEALLOCATE ( hm_rhoy )
DEALLOCATE ( hm_rhox )
DEALLOCATE ( hm_vol )
DEALLOCATE ( hm_etadot )
DEALLOCATE ( hm_theta )
DEALLOCATE ( hm_w )
DEALLOCATE ( hm_v )
DEALLOCATE ( hm_u )
DEALLOCATE ( ec_area )
DEALLOCATE ( ec_vol )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE eg_destroy_helmholtz

END MODULE eg_helmholtz_mod
