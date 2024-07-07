! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE field_types

IMPLICIT NONE

! FLDTYPE definitions for the different field types recognised on the
! decomposition
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

INTEGER, PARAMETER :: Nfld_max           = 8

! maximum number of field types

INTEGER, PARAMETER :: fld_type_unknown   = -1
! non-standard grid

INTEGER, PARAMETER :: fld_type_p         = 1
! grid on P points

INTEGER, PARAMETER :: fld_type_u         = 2
! grid on U points

INTEGER, PARAMETER :: fld_type_v         = 3
! grid on V points

INTEGER, PARAMETER :: fld_type_r         = 7
! grid on river points

INTEGER, PARAMETER :: fld_type_w         = 8
! Only for Endgame

END MODULE field_types
