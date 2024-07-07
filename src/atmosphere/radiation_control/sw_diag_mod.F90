! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE sw_diag_mod

USE max_calls, ONLY: npd_swcall
USE def_diag, ONLY: StrDiag

IMPLICIT NONE

! Declaration of SW diagnostics.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

! Those calculated within the radiation code and not used to
! advance the integration in any way are contained within
! a structure defined in def_diag.

TYPE (StrDiag), SAVE :: sw_diag(npd_swcall)

END MODULE sw_diag_mod
