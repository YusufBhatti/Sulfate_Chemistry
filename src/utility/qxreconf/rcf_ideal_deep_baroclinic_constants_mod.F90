! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Parameters for deep atmosphere baroclinic wave test

MODULE rcf_ideal_deep_baroclinic_constants_mod

USE conversions_mod,       ONLY: pi

IMPLICIT NONE

!
! Description:
!
!   Parameters for Staniforth deep atmosphere
!   baroclinic wave test
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!

REAL, PARAMETER :: u0 = 35.0        ! max wind (m/s)
REAL, PARAMETER :: eg_lapse = 0.005 ! tropo lapse rate (K/m)
REAL, PARAMETER :: xc = pi/9.0      ! lon of perturbation centre (rad)
REAL, PARAMETER :: yc = 2*pi/9.0    ! lat of perturbation centre (rad)
REAL, PARAMETER :: up = 1.0         ! Magnitude of perturbation (m/s)
REAL, PARAMETER :: taper_top = 15000.0 ! Top of u perturbation

REAL, SAVE :: t0
REAL, SAVE :: h, b, c
REAL, SAVE, ALLOCATABLE :: tau1(:), tau2(:),                          &
                           tau1_int(:), tau2_int(:)

INTEGER, SAVE :: k_const_save, b_const_save

LOGICAL, SAVE :: L_shallow_save

END MODULE rcf_ideal_deep_baroclinic_constants_mod
