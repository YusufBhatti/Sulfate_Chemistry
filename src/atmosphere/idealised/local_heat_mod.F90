! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE local_heat_mod

USE missing_data_mod, ONLY: rmdi, imdi

IMPLICIT NONE
!
! Description: Local heating options and parameters
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: IDEALISED
!
! Code description:
!   Language: Fortran 90.
!   This code is written to programming standards in UMDP3. 

INTEGER, PARAMETER :: omit      = 0   ! no local heating
INTEGER, PARAMETER :: analytic  = 1   ! Heating centred on a point like
                                      ! the heating used by Oliver Halliday
                                      ! in an analytic study

INTEGER :: local_heat_option = imdi

REAL    :: local_heat_xoffset = rmdi
REAL    :: local_heat_yoffset = rmdi
REAL    :: local_heat_amp     = rmdi
REAL    :: local_heat_sigma   = rmdi
REAL    :: local_heat_base    = rmdi
REAL    :: local_heat_top     = rmdi
REAL    :: local_heat_period  = rmdi


END MODULE local_heat_mod
