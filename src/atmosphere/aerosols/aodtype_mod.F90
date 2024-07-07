! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
MODULE aodtype_mod

IMPLICIT NONE
! Start AODTYPE

! Description:
!   Sets aerosol type numbers for aerosol optical depth diags.
!
! Current Code Owner: Nicolas Bellouin
!
! A given aerosol type may gather several aerosol components
! (see AERCMP3A)
!
INTEGER, PARAMETER :: ip_type_allaod   = 0
INTEGER, PARAMETER :: ip_type_sulphate = 1
INTEGER, PARAMETER :: ip_type_dust     = 2
INTEGER, PARAMETER :: ip_type_seasalt  = 3
INTEGER, PARAMETER :: ip_type_soot     = 4
INTEGER, PARAMETER :: ip_type_biomass  = 5
INTEGER, PARAMETER :: ip_type_biogenic = 6
INTEGER, PARAMETER :: ip_type_ocff     = 7
INTEGER, PARAMETER :: ip_type_delta    = 8
INTEGER, PARAMETER :: ip_type_nitrate  = 9
INTEGER, PARAMETER :: ip_type_twobdust = 10
! End AODTYPE

END MODULE aodtype_mod
