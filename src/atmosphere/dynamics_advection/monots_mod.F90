! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Options for monotone vertical interpolation in the semi-Lagrangian scheme.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics Advection

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE monots_mod

IMPLICIT NONE

! Tri-Linear monotone scheme interpolation
INTEGER, PARAMETER :: triLinear             = 1

! Quasi-Cubic monotone scheme interpolation
INTEGER, PARAMETER :: mono_quasiCubic       = 2

END MODULE monots_mod
