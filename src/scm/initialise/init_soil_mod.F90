! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: SCM      

MODULE init_soil_mod

! Defines the soil parameters used in the SCM.
! Values in each variable are for soil types Ice, Clay, Loam, Loamy Sand

IMPLICIT NONE

INTEGER, PARAMETER :: nsoilp = 4 
! Number of possible soil parameters

REAL, PARAMETER :: b_exp_typ(nsoilp)  = (/ 0.0,    11.2,   6.6,    3.6   /)
! Single Layer :
!   (C_EAG in code) Eagleson's exponent for calc. sub surf. runoff, P253.4
! Multilayer Hydrology:
!   Exponent used in calculation of soil water suction and hydraulic 
!   conductivity (known as B_WAG in HYDROL2A)
! JULES:
!   Exponent used in calculation of soil water suction and hydraulic
!   conductivity (known as B Clapp-Hornberger exponent)

REAL, PARAMETER :: v_crit_typ(nsoilp) = (/ 0.370,  0.370,  0.332,  0.128 /)
! Volumetric soil moisture content the critical point 
! Below this value evaporation falls below its max
! (m^3 per m^3 soil)                            

REAL, PARAMETER :: v_sat_typ(nsoilp)  = (/ 0.456,  0.456,  0.458,  0.382 /)
! Volumetric soil moisture content at saturation 
! (m^3/m^3 soil) 

REAL, PARAMETER :: v_wilt_typ(nsoilp) = (/ 0.263,  0.263,  0.187,  0.045 /)
! Volumetric soil moisture content at wilting point
! (m^3/m^3)

REAL, PARAMETER :: hcap_typ(nsoilp)   = (/ 0.63E6, 1.23E6, 1.19E6, 1.23E6/)
! Soil heat capacity (J/K/m^3)

REAL, PARAMETER :: hcon_typ(nsoilp)   = (/ 0.265,  0.22,   0.23,   0.32  /)
! Soil thermal conductivity (W/m/K)

REAL, PARAMETER :: sathh_typ(nsoilp)  = (/ 0.0,    0.324,  0.397,  0.062 /)
!     Dummy
! JULES/Multilayer hydrology : 
!     Saturated soil water suction 

REAL, PARAMETER :: satcon_typ(nsoilp) = (/ 0.0,    0.0015, 0.0028, 0.0195/)
! Saturated hydrological conductivity of the soil 
! (kg/m2/s) 

END MODULE init_soil_mod
