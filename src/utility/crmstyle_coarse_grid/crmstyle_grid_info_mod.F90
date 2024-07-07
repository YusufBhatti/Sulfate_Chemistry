! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Holds information on dimensions etc for grids which are not in the namelist

MODULE crmstyle_grid_info_mod

IMPLICIT NONE
SAVE

! Description:
!   Module containing information on dimensions etc for grids
!   Plus PE resources
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: crmstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:

!------------------------------------------------------------------------------
! Information on MPP jobs resources
!------------------------------------------------------------------------------

INTEGER ::         &
  nprocs           & ! total number of processors/nodes available to job
 ,nproc_x          & ! number of processors/nodes to use in x direction
 ,nproc_y          & ! number of processors/nodes to use in y direction
 ,thread_level_set   ! indicates threading type to MPI (i.e. multiple or single)

!------------------------------------------------------------------------------
! Information controlling full and local grids at both the original fine
! resolution and the output coarse resolution.
!------------------------------------------------------------------------------

! Full input grid

!INTEGER ::         &
!  in_cols          & ! number of columns full input data
! ,in_rows            ! number of rows full input data

! Full output coarse grid

INTEGER ::         &
  full_new_x       & ! Number of columns new coarse output grid
 ,full_new_y         ! Number of rows new coarse output grid


INTEGER ::         &
  new_res_x        & ! Number of columns forming new coarse grid
 ,new_res_y          ! Number of rows  forming new coarse grid

! Local fine grid on this processor/node

INTEGER ::         &
  local_row_len    & ! Number of E-W points on this processor
 ,local_rows         ! Number of rows on this processor

! Local new coarse grid on this processor/node
INTEGER ::         &
  local_new_x      & ! Number of new x coarse points on this processor
 ,local_new_y        ! Number of new y coarse points on this processor

! Full fine grid for bcu_mask
INTEGER ::         &
  mask_full_x      & ! Number of columns new high res output grid
 ,mask_full_y        ! Number of rows new high res output grid


END MODULE crmstyle_grid_info_mod
