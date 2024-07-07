! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
MODULE c_sulbdy_mod

IMPLICIT NONE
! Parameters for dry deposition of Sulphur Cycle tracers
! Rb(SO2)/Rb(H2O)
REAL,PARAMETER:: resb_so2=1.53    ! CUBRT(MolWt(SO2)/MolWt(H2O))

! Rb(NH3)/Rb(H2O)
REAL,PARAMETER:: resb_nh3=0.981

! Rb(SO4 modes)/Rb(H2O)
REAL,PARAMETER:: resb_so4_ait=94.9
REAL,PARAMETER:: resb_so4_acc=2530.0
REAL,PARAMETER:: resb_so4_dis=0.0

! Rs(SO2)/Rs(H2O)
REAL,PARAMETER:: ress_so2=1.89    ! SQRT(MolWt(SO2)/MolWt(H2O))

! Rs(NH3)/Rs(H2O)
REAL,PARAMETER:: ress_nh3=0.972

! Rb(SO4 modes)/Rb(H2O)
REAL,PARAMETER:: ress_so4_ait=0.0    ! Valid for small particles
REAL,PARAMETER:: ress_so4_acc=0.0
REAL,PARAMETER:: ress_so4_dis=0.0

! low limit for canopy conductance
REAL,PARAMETER:: cond_lim=1.0e-3

! res to dry dep SO2 over snow
REAL,PARAMETER:: r_snow=1.0e3     ! s/m, from Pasdro et al 1993

! parameter for snow fraction calcn
REAL,PARAMETER:: asnow=0.2         ! m2/kg, from radiation doc.

! C_SULBDY end

END MODULE c_sulbdy_mod
