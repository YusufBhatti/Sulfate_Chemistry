! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module containing a subroutine to initialise module variables 
!          belonging to STASH (and other things not used in the SCM) to zero.  
!          Called from scm_shell.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
! Code Description:
!   Language: Fortran 90
!   This code is written to UMDP3 v9.1 programming standards.
!
MODULE init_scm_misc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='INIT_SCM_MISC_MOD'

CONTAINS

SUBROUTINE init_scm_misc

USE stash_array_mod, ONLY:                                                     &
    n_req_items, num_pseudo_lists, totitems, num_stash_pseudo, nitems,         &
    nstash_series_records, num_level_lists, nstash_series_block,               &
    num_stash_levels, nsttims, nsttabl, nsects, n_ppxrecs

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='INIT_SCM_MISC'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise to zero
nsects      = 0
n_req_items = 0
nitems      = 0
n_ppxrecs   = 0
totitems    = 0
nsttims     = 0
nsttabl     = 0

num_stash_levels      = 0
num_level_lists       = 0
num_stash_pseudo      = 0
num_pseudo_lists      = 0
nstash_series_block   = 0
nstash_series_records = 0

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE init_scm_misc

END MODULE init_scm_misc_mod
