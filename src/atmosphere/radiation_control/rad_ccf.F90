! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Radiation Constant Configuration File

MODULE rad_ccf

! Description:
!   Module containing settings of standard constants used within the
!   radiation scheme. This replaces a number of separate constant
!   configuration files. The original file names are indicated at the
!   head of each section.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Code Description:
!   Language: Fortran 95
!   This code is written to UMDP3 v8.3 programming standards.

! The SOCRATES routines rely on each host model providing the constants
! via a module called rad_ccf. Include common UM parameters here.
USE conversions_mod,      ONLY: pi
USE planet_constants_mod, ONLY: r, repsilon, c_virtual

IMPLICIT NONE

! astron_constants_ccf
! ------------------------------------------------------------------
! Module to set values of astronomical constants.
REAL, PARAMETER :: astronomical_unit = 149597870700.0
!   Standard Astronomical Unit (mean Earth-Sun distance)

! diff_elsasser_ccf, elsass3a
! ------------------------------------------------------------------
! Module to set diffusivity for elsasser's scheme.
REAL, PARAMETER :: elsasser_factor = 1.66e+00
!   Diffusivity factor for elsasser's scheme

! ------------------------------------------------------------------
! diff_keqv_ucf, diffke3a
! ------------------------------------------------------------------
! Module to set the diffusivity factor for use with equivalent
! extinction
REAL, PARAMETER :: diffusivity_factor_minor = 1.66e+00
!   Minor diffusivity factor

! ------------------------------------------------------------------
! physical_constants_0_ccf, phycn03a
! ------------------------------------------------------------------
! Module setting physical constants for the Earth.
REAL, PARAMETER :: mol_weight_air = 28.966e-03
!   Molar weight of dry air
REAL, PARAMETER :: n2_mass_frac   = 0.781e+00
!   Mass fraction of nitrogen

! ------------------------------------------------------------------
! physical_constants_pp_ccf
! ------------------------------------------------------------------
! Module setting physical constants.
REAL, PARAMETER :: n_avogadro = 6.022045e+23
!   Avogadro's number
REAL, PARAMETER :: k_boltzmann = 1.380662e-23
!   Boltzmann's constant
REAL, PARAMETER :: r_gas = 8.3143e+00
!   Universal gas constant

! ------------------------------------------------------------------
! Exoplanet constants
! ------------------------------------------------------------------
! Physical constants for hot Jupiters.
REAL, PARAMETER :: a_h  = 0.91183e+00
!   Number fraction of H
REAL, PARAMETER :: A_He = 1.0e+00 - a_h
!   Number fraction of He

END MODULE rad_ccf
