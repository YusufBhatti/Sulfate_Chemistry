! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
MODULE rcf_ideal_qprofile_mod

IMPLICIT NONE
! Description: module containing q profile types
!  for use in idealised problems

INTEGER, PARAMETER :: qp_dry=0
INTEGER, PARAMETER :: qp_qsat=1
INTEGER, PARAMETER :: qp_namelist_rh=8
INTEGER, PARAMETER :: qp_namelist=9
INTEGER, PARAMETER :: qp_dump=10
INTEGER, PARAMETER :: qp_BS1999=11

END MODULE rcf_ideal_qprofile_mod
