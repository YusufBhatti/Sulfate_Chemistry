! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE tforce_mod

IMPLICIT NONE
!
! Description: Options for temperature profiles to be used
!              with idealised forcing.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: IDEALISED
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standards in UMDP3. 


INTEGER, PARAMETER :: tf_none=0
INTEGER, PARAMETER :: tf_HeldSuarez=1
INTEGER, PARAMETER :: tf_TLE=2
INTEGER, PARAMETER :: tf_EL=3
INTEGER, PARAMETER :: tf_SHJ=4
INTEGER, PARAMETER :: tf_Jupiter=5
INTEGER, PARAMETER :: tf_HD209458b_Heng=6
INTEGER, PARAMETER :: tf_HD209458b_Heng_smooth=7
INTEGER, PARAMETER :: tf_HD209458b_iro=8
INTEGER, PARAMETER :: tf_Y_Dwarf=9
INTEGER, PARAMETER :: tf_GJ1214b=10
INTEGER, PARAMETER :: tf_GJ1214b_dT800=11
INTEGER, PARAMETER :: tf_file=99
INTEGER, PARAMETER :: tf_isothermal=100

END MODULE tforce_mod
