! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
MODULE c_micro_mod

IMPLICIT NONE
!
! Description:
!
!  Contains various cloud droplet parameters, defined for
!  land and sea areas.
!
!  KPARAM_* is the ratio of the cubes of the volume-mean radius
!                                           and the effective radius;
!  DCONRE_* is the effective radius (m) for deep convective clouds;
!  DEEP_CONVECTION_LIMIT_* is the threshold depth (m) bewteen shallow
!                                          and deep convective cloud.
!
! Current Code Owner: Andy Jones
!
!
REAL,PARAMETER:: kparam_land = 0.67
REAL,PARAMETER:: kparam_sea = 0.80
REAL,PARAMETER:: dconre_land = 9.5e-06
REAL,PARAMETER:: dconre_sea = 16.0e-06
REAL,PARAMETER:: deep_convection_limit_land = 500.0
REAL,PARAMETER:: deep_convection_limit_sea = 1500.0
!
! Maximum Temp for homogenous nucleation (deg C)
! Thomo has moved back into mphys_ice_mod
! To call this in your code use the following syntax:
! USE mphys_ice_mod, ONLY: thomo

END MODULE c_micro_mod
