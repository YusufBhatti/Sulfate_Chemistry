! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Module containing error codes

MODULE Err_Mod

! Description:
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle coarse grid
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: IOSTAT_END

IMPLICIT NONE

! Error codes
INTEGER, PARAMETER :: StatusOK      =  0
INTEGER, PARAMETER :: StatusWarning = -9
INTEGER, PARAMETER :: StatusFatal   =  9
INTEGER, PARAMETER :: EndofFile     = IOSTAT_END

END MODULE Err_Mod
