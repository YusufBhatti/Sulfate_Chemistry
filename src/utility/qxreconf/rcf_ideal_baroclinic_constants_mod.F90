! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Parameters for a baroclinic wave test

MODULE rcf_ideal_baroclinic_constants_mod

USE conversions_mod,     ONLY: pi

IMPLICIT NONE

!
! Description:
!
!   Parameters for Jablonowski & Williamson (2006)
!   baroclinic wave test: QJRMS 132, 2943--2975
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UMDP3 programming standards.
!

REAL, PARAMETER :: eta_t = 0.2      ! tropopause level
REAL, PARAMETER :: eta0 = 0.252     ! zonal jet level
REAL, PARAMETER :: u0 = 35.0        ! max wind (m/s)
REAL, PARAMETER :: dt = 4.8e5       ! empirical temp diffnce in strat (K)
REAL, PARAMETER :: t0 = 288.0       ! mean surface temperature (K)
REAL, PARAMETER :: eg_lapse = 0.005 ! tropo lapse rate (K/m)
REAL, PARAMETER :: xc = pi/9.0      ! lon of perturbation centre (rad)
REAL, PARAMETER :: yc = 2*pi/9.0    ! lat of perturbation centre (rad)
REAL, PARAMETER :: up = 1.0         ! Magnitude of perturbation (m/s)

! Parameters for channel test
LOGICAL         :: L_channel
REAL, PARAMETER :: b = 2.0
REAL            :: Lx
REAL            :: Ly
REAL            :: Lp
REAL            :: f0
REAL            :: beta0
REAL            :: baro_phi0

END MODULE rcf_ideal_baroclinic_constants_mod
