! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE rcf_ideal_pprofile_constants_mod

IMPLICIT NONE
!
! Description: Parameters for pressure balance choices
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standards in UMDP3. 

! Options for pressure balance
INTEGER, PARAMETER :: no_balance  = 0
INTEGER, PARAMETER :: dry_hydro   = 1
INTEGER, PARAMETER :: moist_hydro = 2

END MODULE rcf_ideal_pprofile_constants_mod
