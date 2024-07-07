! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Source for Gravity Wave (Ultra-Simple Spectral Parametrization) Scheme.

MODULE cgw_input_grab_mod

! Description:
!       Make convection data available for input to gw_ussp scheme.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Deprecated
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.
!   Documentation: Deprecated
!
! No longer required as the data is now passed round timestep via D1.
! Once testing is complete, will want to remove module entirely.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CGW_INPUT_GRAB_MOD'

END MODULE cgw_input_grab_mod
