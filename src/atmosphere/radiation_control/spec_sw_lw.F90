! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module defines the LW & SW spectral file data for each call
!   to radiation.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Code Description:
!   Language: FORTRAN 90

MODULE spec_sw_lw

! Define the type of the spectrum

USE max_calls
USE def_spectrum

IMPLICIT NONE

TYPE (StrSpecData), SAVE :: sw_spectrum(npd_swcall)
TYPE (StrSpecData), SAVE :: lw_spectrum(npd_lwcall)

END MODULE spec_sw_lw
