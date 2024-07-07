! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Description
!   This module contains the variables which contain the
!   atmosphere land sea mask - both the full field, and the
!   local subdomain on this processor.
!   This data is required for various compression/decompression
!   algorithms.
!

MODULE atm_land_sea_mask

IMPLICIT NONE

LOGICAL, ALLOCATABLE ::                                                        &
  atmos_landmask(:)                 ! Full-grid land-sea mask

LOGICAL, ALLOCATABLE ::                                                        &
  atmos_landmask_local(:)           ! Local subdomain area lsm

INTEGER ::                                                                     &
  atmos_number_of_landpts           ! total number of land points

INTEGER, ALLOCATABLE ::                                                        &
  atmos_number_of_landpts_proc(:)   ! Per PE number of land points

END MODULE atm_land_sea_mask

