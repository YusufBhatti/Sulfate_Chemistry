! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE surface_flux_mod

! Purpose:
!   Contains information on surface flux forcing for the idealised model.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

USE missing_data_mod, ONLY: rmdi, imdi

IMPLICIT NONE

SAVE

! Parameters
! ======================
! Maximum size of arrays that can be read through namelist
INTEGER, PARAMETER :: num_surf_max  = 100

INTEGER, PARAMETER :: zero_flux       = 1
INTEGER, PARAMETER :: diurnal_flux    = 2
INTEGER, PARAMETER :: constant_flux   = 3
INTEGER, PARAMETER :: hot_spot        = 4
INTEGER, PARAMETER :: time_varying    = 5


! Variables from idealised namelist
! ------------------------------------

INTEGER :: IdlSurfFluxSeaOption = imdi

REAL    :: IdlSurfFluxseaParams(4) = rmdi

! Time varying surface fluxes
INTEGER :: num_surface_flux_times   = imdi

REAL    :: surface_flux_time(num_surf_max) = rmdi
REAL    :: sh_flux(num_surf_max)           = rmdi
REAL    :: lh_flux(num_surf_max)           = rmdi


END MODULE surface_flux_mod
