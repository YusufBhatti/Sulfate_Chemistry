! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
MODULE rcf_ideal_initial_profiles_mod

IMPLICIT NONE
! Description: module containing values for choosing a preset
!     baroclinic/wave test case.

INTEGER, PARAMETER :: no_preset      = 0  ! None; use tprofile_number
INTEGER, PARAMETER :: baro_inst      = 1  ! Baroclinic test case
INTEGER, PARAMETER :: rot_solid_body = 2  ! Rot. solid body solution
INTEGER, PARAMETER :: deep_baro_inst = 3  ! Deep baroclinic test case

END MODULE rcf_ideal_initial_profiles_mod
