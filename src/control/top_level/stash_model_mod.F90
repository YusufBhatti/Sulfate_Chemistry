! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Purpose: Module defines model-dependent quantities used by
!          data addressing and STASH
!

MODULE stash_model_mod

USE submodel_mod, ONLY: n_internal_model_max, n_submodel_partition_max

IMPLICIT NONE

SAVE
PRIVATE :: n_internal_model_max, n_submodel_partition_max

! Horizontal Atmosphere Grid Variables
REAL :: h_a_ewspace      ! E-W grid spacing (i.e. dx), (units?)
REAL :: h_a_nsspace      ! N-S spacing (i.e. dy), (units?)
REAL :: h_a_firstlat     ! Latitude of first point
REAL :: h_a_firstlong    ! Longitude of first point
REAL :: h_a_polelat      ! Latitude of pole co-ordinates
REAL :: h_a_polelong     ! Longitude of pole co-ordinates

INTEGER :: h_orog_rough  ! Comment?


! Data for each submodel data partition
!============================================
INTEGER :: len_prim (n_submodel_partition_max)  ! Tot. len for primary fields
INTEGER :: len_dump (n_submodel_partition_max)  ! Tot. len for diag fields
INTEGER :: len_secd (n_submodel_partition_max)  ! Tot. len for secondary fields
INTEGER :: nheadsub (n_submodel_partition_max)  ! Tot. number of headers
                                                ! (i.e. levels)

! Global versions of len_prim/len_dump
INTEGER :: global_len_prim (n_submodel_partition_max) 
                                                ! Tot. length of primary fields
INTEGER :: global_len_dump (n_submodel_partition_max)
                                                ! Tot. length of diag fields

INTEGER :: len_work (n_submodel_partition_max)  ! Total workspace length
INTEGER :: len_extra (n_submodel_partition_max) ! Total length of extra space


! Data for each internal model
!=============================================================
INTEGER :: len_primim (n_internal_model_max) ! Tot. len for primary fields
INTEGER :: len_dumpim (n_internal_model_max) ! Tot. len for diag fields
INTEGER :: len_secdim (n_internal_model_max) ! Tot. len for secondary fields
INTEGER :: nhead      (n_internal_model_max) ! Tot. number of headers
                                             ! (i.e. levels)

! Global versions of len_primim/len_dumpim
INTEGER :: global_len_primim(n_internal_model_max)
                                               ! Tot. length for primary fields
INTEGER :: global_len_dumpim(n_internal_model_max)
                                               ! Tot. length for diag fields

INTEGER :: item_max_req
INTEGER :: item_max_all

INTEGER :: nrecs_s
INTEGER :: ntimes_s
INTEGER :: nserblk_s
INTEGER :: nserrec_s
INTEGER :: nlevl_s
INTEGER :: nmaxlev_s
INTEGER :: npslists_s
INTEGER :: nmaxpsl_s
LOGICAL :: lstuser

CHARACTER(LEN=1) :: h_atmos
CHARACTER(LEN=1) :: h_strat
CHARACTER(LEN=1) :: h_global(n_internal_model_max)

INTEGER :: mean_number(n_internal_model_max)
INTEGER :: stlevgwdrag
INTEGER :: botvdifflev
INTEGER :: topvdifflev

END MODULE stash_model_mod
