! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE fieldsfile_constants_mod

IMPLICIT NONE

! Description:
!    Module to contain data and routines related to the fieldsfile format.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! Fixed Header(1) - dataset format version number
INTEGER, PARAMETER :: fixed_header_1 = 20

! Fixed Header(9) - grid staggering values
INTEGER, PARAMETER :: new_dynamics = 3
INTEGER, PARAMETER :: endgame = 6


END MODULE fieldsfile_constants_mod
