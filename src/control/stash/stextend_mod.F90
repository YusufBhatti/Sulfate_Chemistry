! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash

MODULE stextend_mod

USE version_mod, ONLY: nprofdp, nsectp, nitemp,  &
                       nrecdp, ntimep, nproftp, nelemp,    &
                       nlevp_s, npslistp

USE submodel_mod, ONLY: n_internal_model_max, n_submodel_partition_max

IMPLICIT NONE

! Description:
! Contains variables and arrays involved in STASH processing in the UM.

!   Output levels lists
!     List type (real/int)
CHARACTER(LEN=1) :: llistty(nprofdp)
!     Real levels
REAL ::     rlevlst_s(nlevp_s, nprofdp)
!     Integer (i.e. model) levels
INTEGER ::   levlst_s(nlevp_s, nprofdp)

!   STASH list array (extra row only for internal indexing in
!                   processing routines)
INTEGER :: list_s  (nelemp+1, nrecdp)
!   Output times tables
INTEGER :: itim_s  (ntimep, 2*nproftp+2)

!   STASH lengths and addresses
INTEGER :: in_s    (2,N_INTERNAL_Model_MAX,0:nsectp,nitemp)
!   STASH list index
INTEGER :: indx_s  (2,N_INTERNAL_Model_MAX,0:nsectp,nitemp)

!   Start addresses for pp headers
INTEGER :: ppind_s (N_INTERNAL_Model_MAX,nitemp)
!   Time series block information
!     No. of records in a block
INTEGER :: nrecs_ts(nprofdp)
!     Start position of block
INTEGER :: npos_ts (nprofdp)
!   lengths of pseudo-levels lists
INTEGER :: lenplst (npslistp)

!     Set up preliminary array for addressing D1:
!     Number of items of info needed for each object and likely maximum
!     number of objects in D1 - this can be increased if necessary

INTEGER, PARAMETER :: D1_Items_Prel = 5
INTEGER, PARAMETER :: Max_D1_Len    = 2200

! Names of items

INTEGER, PARAMETER :: d1_type       = 1 ! Prognostic, diagnostic
                                        ! or other
INTEGER, PARAMETER :: d1_im         = 2 ! Internal model id
INTEGER, PARAMETER :: d1_extra_info = 3 ! Progs and other :-
                                        ! PPXREF item no
                                        ! Diags :-
                                        ! Stash list item no
INTEGER, PARAMETER :: d1_levs       = 4 ! No of levels
INTEGER, PARAMETER :: d1_sect       = 5 ! Section No

! Types of items for d1_type

INTEGER, PARAMETER :: Prog     = 0
INTEGER, PARAMETER :: Diag     = 1
INTEGER, PARAMETER :: Seco     = 2
INTEGER, PARAMETER :: Extra_d1 = 3

! Stores number of objects in D1
INTEGER :: N_Obj_D1(N_SUBMODEL_Partition_MAX)

!     Preliminary array for addressing D1. Holds the minimum amount of
!     info required for order of objects of D1; this amount of info is
!     enough to obtain any other required info from stashlist or ppxref

INTEGER :: D1_PAddr(D1_Items_Prel,Max_D1_Len,N_SUBMODEL_Partition_MAX)

END MODULE stextend_mod
