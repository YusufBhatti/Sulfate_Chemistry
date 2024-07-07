! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE trelax_mod
IMPLICIT NONE
! Description: Module containing integer to select the
!   relaxation timescale for newtonian forcing setups
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: IDEALISED
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

      INTEGER, PARAMETER :: tr_none=0
      INTEGER, PARAMETER :: tr_HeldSuarez=1
      INTEGER, PARAMETER :: tr_EL=2
      INTEGER, PARAMETER :: tr_SHJ=3
      INTEGER, PARAMETER :: tr_Jupiter=4
      INTEGER, PARAMETER :: tr_HD209458b_Iro=5
      INTEGER, PARAMETER :: tr_Y_Dwarf=6
      INTEGER, PARAMETER :: tr_GJ1214b=7

END MODULE trelax_mod
