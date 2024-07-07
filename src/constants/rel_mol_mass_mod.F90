! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! A module containing constants/parameters for Relative Molecular Mass
!
MODULE rel_mol_mass_mod

IMPLICIT NONE

! Description:
!   This module contains constants for Relative Molecular Mass
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Constants
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.
!----------------------------------------------------------------------

! Relative Molecular Mass (kg/mole)
REAL, PARAMETER :: rmm_s = 3.20e-2    ! S
REAL, PARAMETER :: rmm_h2o2 = 3.40e-2 ! H2O2
REAL, PARAMETER :: rmm_o3 = 4.8e-2    ! O3
REAL, PARAMETER :: rmm_air = 2.896e-2 ! dry air
REAL, PARAMETER :: rmm_w = 1.8e-2     ! water

END MODULE rel_mol_mass_mod

