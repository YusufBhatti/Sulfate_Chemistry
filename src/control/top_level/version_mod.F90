! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Purpose: STASH parameter definitions

MODULE version_mod

IMPLICIT NONE

! Max. no. of STASH sections per internal model
INTEGER,PARAMETER :: nsectp=99
! Max. no. of STASH items per section
INTEGER,PARAMETER :: nitemp=999
! Max. no. of STASH list records (prognostic + diagnostic)
INTEGER,PARAMETER :: nrecdp=4500
! Max. no. of output times tables in STASHC
INTEGER,PARAMETER :: ntimep=100
! Max. no. of time profiles in STASHC
INTEGER,PARAMETER :: nproftp=100
! Max. no. of domain profiles/levels lists in STASHC (used for both)
INTEGER,PARAMETER :: nprofdp=100
! Max. total no. of time series in STASHC
INTEGER,PARAMETER :: NTimSerP=1500
! Max. no. time series per domain profile
INTEGER,PARAMETER :: tsdp=250
! Max. no. of useage profiles in STASHC
INTEGER,PARAMETER :: nprofup=40
! Max. no. of packages in STASHC
INTEGER,PARAMETER :: nprofpp=30
! Max. no. of levels in a levels list
INTEGER,PARAMETER :: nlevp=50
! Max. no. of pseudo levels in a  pseudo levels list
INTEGER,PARAMETER :: npslevp=5000
! Max. no. of pseudo levels lists in STASHC
INTEGER,PARAMETER :: npslistp=150
! Max. no. non-blank records in PPXREF file
INTEGER,PARAMETER :: ndiagp=5000
INTEGER,PARAMETER :: ndiagpm=nrecdp
INTEGER,PARAMETER :: nelemp=39
INTEGER,PARAMETER :: nlevp_S=nlevp*6+1
INTEGER,PARAMETER :: nlevlstsp=nprofdp

END MODULE version_mod
