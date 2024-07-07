! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lw_diag_mod

USE max_calls, ONLY: npd_lwcall
USE def_diag, ONLY: StrDiag

IMPLICIT NONE

! Declaration of LW diagnostics.
! Those quantities which are purely diagnostic are bundled into
! a structure.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

TYPE (StrDiag), SAVE :: lw_diag(npd_lwcall)

END MODULE lw_diag_mod
