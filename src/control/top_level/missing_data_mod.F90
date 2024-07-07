! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

! Data module containing missing data indicies

! Migrated from include file c_mdi.h

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level
!
! Code Description:
!   Language: FORTRAN 90
!

MODULE missing_data_mod

IMPLICIT NONE

     ! PP missing data indicator (-1.0E+30)
REAL, PARAMETER    :: rmdi_pp  = -1.0e+30

! Old real missing data indicator (-32768.0)
REAL, PARAMETER    :: rmdi_old = -32768.0

! New real missing data indicator (-2**30)
REAL, PARAMETER    :: rmdi     = -32768.0*32768.0

! Integer missing data indicator
INTEGER, PARAMETER :: imdi     = -32768

END MODULE missing_data_mod
