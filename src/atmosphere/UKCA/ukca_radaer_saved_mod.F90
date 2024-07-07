! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Module file ukca_radaer_saved_mod
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA

! UKCA_RADAER:
!
! ukca_radaer is saved here

MODULE ukca_radaer_saved_mod

USE ukca_radaer_struct_mod, ONLY: &
    ukca_radaer_struct

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UKCA_RADAER_STRUCT_MOD'

! Structure for UKCA/radiation interaction
TYPE (ukca_radaer_struct), SAVE :: ukca_radaer

END MODULE ukca_radaer_saved_mod
