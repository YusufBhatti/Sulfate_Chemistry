! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to store increments from the Leonard terms.

MODULE leonard_incs_mod

IMPLICIT NONE

!
! Description:
!   Module stores increments from the Leonard terms, so that they can
!   just be calculated on the first EG outer loop call to atmos_physics2
!   and then just re-used on subsequent calls.
!
! Method:
!   Module contains the allocatable increment arrays.  If Leonard terms
!   are active, the arrays are allocated by the contained alloc subroutine
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Diffusion and filtering
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName =                    &
                                        'LEONARD_INCS_MOD'

! Arrays to store Leonard term increments of:
REAL, ALLOCATABLE :: u_inc_leonard(:,:,:)      ! Zonal wind
REAL, ALLOCATABLE :: v_inc_leonard(:,:,:)      ! Meridional wind
REAL, ALLOCATABLE :: w_inc_leonard(:,:,:)      ! Vertical wind
REAL, ALLOCATABLE :: thetal_inc_leonard(:,:,:) ! Water potential temp
REAL, ALLOCATABLE :: qw_inc_leonard(:,:,:)     ! Total water content


CONTAINS


! Subroutine to allocate the above arrays; called in the first
! call to atmos_physics2, by atmo_physics2_alloc
SUBROUTINE leonard_incs_alloc()

USE atm_fields_bounds_mod, ONLY: udims, vdims, wdims, tdims
USE nlsizes_namelist_mod,  ONLY: bl_levels

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LEONARD_INCS_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ALLOCATE( u_inc_leonard      ( udims%i_start:udims%i_end,               &
                               udims%j_start:udims%j_end,               &
                               bl_levels ) )

ALLOCATE( v_inc_leonard      ( vdims%i_start:vdims%i_end,               &
                               vdims%j_start:vdims%j_end,               &
                               bl_levels ) )

ALLOCATE( w_inc_leonard      ( wdims%i_start:wdims%i_end,               &
                               wdims%j_start:wdims%j_end,               &
                               bl_levels ) )

ALLOCATE( thetal_inc_leonard ( tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                               bl_levels ) )

ALLOCATE( qw_inc_leonard     ( tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                               bl_levels ) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE leonard_incs_alloc

END MODULE leonard_incs_mod
