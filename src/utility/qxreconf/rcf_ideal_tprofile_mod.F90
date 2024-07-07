! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
MODULE rcf_ideal_tprofile_mod

IMPLICIT NONE
! Description: module containing temperature profile types
!  for use in idealised problems
!
INTEGER, PARAMETER :: tp_dthetadz=1
INTEGER, PARAMETER :: tp_isothermal=2
INTEGER, PARAMETER :: tp_BruntV=3
INTEGER, PARAMETER :: tp_BV_isoth=4
INTEGER, PARAMETER :: tp_dyn_core=6
INTEGER, PARAMETER :: tp_dyn_core_lam=7
INTEGER, PARAMETER :: tp_namelist=9
INTEGER, PARAMETER :: tp_dump=10
INTEGER, PARAMETER :: tp_HD209458b_Heng_mid=11
INTEGER, PARAMETER :: tp_HD209458b_Heng_smooth_mid=12
INTEGER, PARAMETER :: tp_HD209458b_iro=13
INTEGER, PARAMETER :: tp_Y_Dwarf=14
INTEGER, PARAMETER :: tp_file=15
INTEGER, PARAMETER :: tp_HD209458b_Heng_smooth_day=16
INTEGER, PARAMETER :: tp_GJ1214b=17

END MODULE rcf_ideal_tprofile_mod
