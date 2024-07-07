! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
MODULE c_nitrate_bdy_mod

IMPLICIT NONE
!----------------C_NITRATE_BDY------------------------------------
! Parameters for dry deposition of ammonium nitrate tracers.
! Use same vales as for sulphate for the moment.
!
      ! Rb(NH4NO3 modes)/Rb(H2O)
REAL,PARAMETER:: resb_nitr_acc=2530.0
REAL,PARAMETER:: resb_nitr_diss=0.0
!
      ! Rs(NH4NO3 modes)/Rs(H2O)
REAL,PARAMETER:: ress_nitr_acc=0.0
REAL,PARAMETER:: ress_nitr_diss=0.0
!
! C_NITRATE_BDY end

END MODULE c_nitrate_bdy_mod
