! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: stochastic physics
MODULE fp_mod

 ! Global data module for storing the SPT 3-D forcing pattern psif,
 ! to be used across the different parameterization schemes to perturb
 ! their tendencies.
IMPLICIT NONE

REAL, ALLOCATABLE ::                                                    &
psif(:,:,:)                                                             &
! 1st forcing pattern for SPT slow physics
,   psif2(:,:)
! 2nd forcing pattern for SPT fast physics (180 rotation of the 1st)
END MODULE fp_mod
