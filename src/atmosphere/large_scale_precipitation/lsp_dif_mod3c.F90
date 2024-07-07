! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lsp_dif_mod

! Description:
! The apb and tw parameters represent diffusional growth constants
! and wet bulb temperature parameters for use in LSP.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE


! Values of reference variables
REAL,PARAMETER:: air_density0=1.0             ! kg m-3
REAL,PARAMETER:: air_viscosity0=1.717e-5      ! kg m-1 s-1
REAL,PARAMETER:: air_conductivity0=2.40e-2    ! J m-1 s-1 K-1
REAL,PARAMETER:: air_diffusivity0=2.21e-5     ! m2 s-1
REAL,PARAMETER:: air_pressure0=1.0e5          ! Pa

! Values of constants used in correction terms
REAL, PARAMETER :: tcor1 = 393.0
REAL, PARAMETER :: tcor2 = 120.0
REAL, PARAMETER :: cpwr  = 1.5

! Values of diffusional growth parameters (set in lspcon)
! Terms in deposition and sublimation
REAL :: apb1 ! =(lc+lf)**2 * repsilon /(r*air_conductivity0)
REAL :: apb2 ! =(lc+lf) / air_conductivity0
REAL :: apb3 ! =r/(repsilon*air_pressure0*air_diffusivity0)
! Terms in evap of melting snow and rain
REAL :: apb4 ! =lc**2*repsilon/(r*air_conductivity0)
REAL :: apb5 ! =lc /air_conductivity0
REAL :: apb6 ! =r/(repsilon*air_pressure0*air_diffusivity0)

! Values of numerical approximation to wet bulb temperature
! Numerical fit to wet bulb temperature
REAL,PARAMETER:: tw1=1329.31
REAL,PARAMETER:: tw2=0.0074615
REAL,PARAMETER:: tw3=0.85e5
! Numerical fit to wet bulb temperature
REAL,PARAMETER:: tw4=40.637
REAL,PARAMETER:: tw5=275.0

! Ventilation parameters
REAL,PARAMETER:: sc=0.6
! f(v)  =  vent_ice1 + vent_ice2  Sc**(1/3) * Re**(1/2)
REAL,PARAMETER:: vent_ice1=0.65
REAL,PARAMETER:: vent_ice2=0.44
! f(v)  =  vent_rain1 + vent_rain2  Sc**(1/3) * Re**(1/2)
REAL,PARAMETER:: vent_rain1=0.78
REAL,PARAMETER:: vent_rain2=0.31

END MODULE lsp_dif_mod
