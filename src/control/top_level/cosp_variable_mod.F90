! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  COSP

!  Description:

! Fortran module to control cosp variables at a top level for memory use
! efficiency.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

! Code Description:
!   Language: FORTRAN 90

MODULE cosp_variable_mod

USE cosp_types_mod, ONLY: cosp_gridbox

IMPLICIT NONE

! COSP gridbox information. Input for COSP
TYPE(cosp_gridbox) :: cosp_gbx

END MODULE cosp_variable_mod
