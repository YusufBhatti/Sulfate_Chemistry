! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
MODULE rcf_ideal_gflat_mod

IMPLICIT NONE
! Description: module containing vertical coordinate flattening types
!  for use in idealised problems

INTEGER, PARAMETER :: gflat_old=1
INTEGER, PARAMETER :: gflat_linear=0
INTEGER, PARAMETER :: gflat_linear1=2
INTEGER, PARAMETER :: gflat_quadratic=3

END MODULE rcf_ideal_gflat_mod
