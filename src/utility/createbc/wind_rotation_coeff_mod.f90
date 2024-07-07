! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE wind_rotation_coeff_mod

IMPLICIT NONE

! Description:
!   A base class to hold wind rotation coefficients
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

TYPE, PUBLIC :: wind_rotation_coeff_type

  REAL, PUBLIC, ALLOCATABLE :: input_wind_coeff1(:)
  REAL, PUBLIC, ALLOCATABLE :: input_wind_coeff2(:)
  REAL, PUBLIC, ALLOCATABLE :: output_wind_coeff1(:)
  REAL, PUBLIC, ALLOCATABLE :: output_wind_coeff2(:)

END TYPE wind_rotation_coeff_type

END MODULE wind_rotation_coeff_mod
