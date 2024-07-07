! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Holds constants required for the electric scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

MODULE electric_constants_mod

USE conversions_mod, ONLY: zerodegc

IMPLICIT NONE

INTEGER :: i    ! x pointer
INTEGER :: j    ! y pointer
INTEGER :: k    ! z pointer

! Conversion terms
REAL, PARAMETER :: min2sec = 1.0 / 60.0  ! minutes in a second
REAL, PARAMETER :: gwp_thresh = 0.2      ! 200 g m-3 - threshold for a 'storm'

! Graupel Water Path method terms
REAL, PARAMETER :: minus5 = zerodegc - 5.0 ! Minus 5 level

! McCaul et al (2009) parameters
REAL, PARAMETER :: minus15 = zerodegc - 15.0 ! Minus 15 C level for McCaul
REAL, PARAMETER :: mccaul_r1 = 0.95 ! Factor for inclusion of flash1 in McCaul
                                    ! scheme
REAL, PARAMETER :: mccaul_r2 = 0.05 ! Factor for inclusion of flash2 in McCaul
                                    ! scheme
                                    ! Note: mccaul_r1 + mccaul_r2 should always
                                    !       equal 1.0 exactly.

END MODULE electric_constants_mod

