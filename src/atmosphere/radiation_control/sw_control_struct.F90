! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Description:
!   This module declares the controlling structure for SW radiation.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

MODULE sw_control_struct

USE max_calls,   ONLY: npd_swcall
USE def_control, ONLY: StrCtrl

IMPLICIT NONE

TYPE (StrCtrl), SAVE :: sw_control(npd_swcall)

END MODULE sw_control_struct
