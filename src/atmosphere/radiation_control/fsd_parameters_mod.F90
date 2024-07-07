! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Global module for storing parameters used in parametrization of 
! fractional standard deviation (FSD) of subgrid water content.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!----------------------------------------------------------------------

MODULE fsd_parameters_mod

IMPLICIT NONE
SAVE

! Parameters for FSD parametrization

REAL, ALLOCATABLE :: f_arr(:,:,:,:)
!   Array of regime and layer thickness dependent
!   parameters required for parametrization of FSD.

REAL, ALLOCATABLE :: f_arr_c(:,:)
!   Array of regime and layer thickness dependent
!   parameters required for parametrization of FSD.
!   Gathered version for use in microphysics.

REAL :: f_cons(3)
!   Constant parameters required for parametrization of FSD

REAL :: fsd_eff_lam ! effective delta_lambda for fsd parametrization
REAL :: fsd_eff_phi ! effective delta_phi for fsd parametrization

INTEGER, PARAMETER :: ip_fsd_constant     = 0
!   flag to use constant value of FSD/scaling factor
INTEGER, PARAMETER :: ip_fsd_param        = 1
!   flag to use FSD parametrisation described in Hill et al (2012)
!   doi: 10.1002/qj.1893
INTEGER, PARAMETER :: ip_fsd_regime        = 2
!   flag to use regime dependent FSD parametrisation
INTEGER, PARAMETER :: ip_fsd_regime_no_sh  = 3
!   flag to use regime dependent FSD parametrisation
!   doesn't include shallow convection
INTEGER, PARAMETER :: ip_fsd_boutle        = 4
!   flag to use FSD parametrisation described in Boutle et al (2013)
!   doi: 10.1002/qj.2140
INTEGER, PARAMETER :: ip_fsd_regime_smooth = 5
!   flag to use regime dependent FSD parametrisation on smoothed
!   cca field
INTEGER, PARAMETER :: ip_fsd_regime_smooth_no_sh  = 6
!   flag to use regime dependent FSD parametrisation on smoothed
!   cca field. Doesn't include shallow convection

END MODULE fsd_parameters_mod
