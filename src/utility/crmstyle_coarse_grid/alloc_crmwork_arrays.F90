! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! ALLOCATE full grid arrays
MODULE alloc_crmwork_arrays_mod

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ALLOC_CRMWORK_ARRAYS_MOD'

CONTAINS

SUBROUTINE alloc_crmwork_arrays(in_cols, in_rows, in_lev)

USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels

USE crmwork_arrays_mod, ONLY: orog_full, xcoslat_full,        &
  r_theta_levels_local, r_rho_levels_local,                   &
  th_km_level, uv_km_level, th_km_level_full,  mask,          &
  th_weight, uv_weight, th_weight_full,                       &
  index_col,index_row, prec_full

USE crmstyle_cntl_mod,  ONLY: num_x, num_y, new_res, mlevs,   &
                              l_class_col

USE crmstyle_grid_info_mod,  ONLY: local_new_x, local_new_y, &
           local_row_len, local_rows

USE UM_ParCore, ONLY: mype

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   ALLOCATE and initialise work arrays required by crmstyle_coarse_grid
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Utilty - crmstyle_coarse_grid
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------
INTEGER, INTENT(IN) ::  &
  in_cols               & ! Number of columns
 ,in_rows               & ! Number of local_rows
 ,in_lev                  ! Number of levels

!-------------------------------------------------------------------------
! Local variables
INTEGER ::               &
  i, j, k,ii                  ! loop counters

INTEGER, PARAMETER ::  &
  offx  = 1            &  ! one point halos
 ,offy  = 1               !

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOC_CRMWORK_ARRAYS'

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------
! Orography etc for PE 0 and processing of fields read in before scatering
! fields across PEs


ALLOCATE (xcoslat_full(in_cols,in_rows) )

ALLOCATE(orog_full(in_cols,in_rows) )

IF (l_class_col) THEN
  ALLOCATE(prec_full(in_cols,in_rows) )
END IF

ALLOCATE (th_weight_full(in_cols,in_rows,mlevs) )
ALLOCATE (th_km_level_full(in_cols,in_rows,mlevs) )

! All PEs

ALLOCATE (th_km_level(local_row_len,local_rows,mlevs) )
ALLOCATE (uv_km_level(local_row_len,local_rows,mlevs) )

ALLOCATE (th_weight(local_row_len,local_rows,mlevs) )
ALLOCATE (uv_weight(local_row_len,local_rows,mlevs) )

ALLOCATE (mask(local_row_len,local_rows,mlevs) )

! Allocate index arrays
ALLOCATE(index_col(local_row_len) )
ALLOCATE(index_row(local_rows) )

! Allocate level arrays full grid

ALLOCATE(r_theta_levels(in_cols,in_rows,0:in_lev) )
ALLOCATE(r_rho_levels(in_cols,in_rows,in_lev) )

ALLOCATE(r_theta_levels_local(local_row_len,local_rows,0:in_lev) )
ALLOCATE(r_rho_levels_local(local_row_len,local_rows,in_lev) )

! initialise  index arrays
DO j = 1,local_new_x
  DO i = 1,new_res(1)
    ii = (j-1)*new_res(1)+i
    index_col(ii) = j
  END DO
END DO
! initialise  index arrays
DO j = 1,local_new_y
  DO i = 1,new_res(1)
    ii = (j-1)*new_res(1)+i
    index_row(ii) = j
  END DO
END DO

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE alloc_crmwork_arrays
END MODULE alloc_crmwork_arrays_mod
