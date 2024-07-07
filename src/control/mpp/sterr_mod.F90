! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: MPP
MODULE sterr_mod

IMPLICIT NONE
!    Purpose: PARAMETER names for STASH processing error codes;
!             fatal errors have positive codes, warnings negative.
!
!    Programming standard: UM Doc Paper 3, version 2 (7/9/90)
!
!    Logical components covered: D70
!
!    Project task: D7
!
!    External documentation:
!      Unified Model Doc Paper C4 - Storage handling and diagnostic
!                                   system (STASH)
!
! Warning codes
!
INTEGER :: st_upper_less_lower ! warning code for bad domain
PARAMETER(st_upper_less_lower=-1)

INTEGER :: st_not_supported ! warning code for unsupported routine
PARAMETER(st_not_supported=-2)
INTEGER :: st_no_data,st_nd ! indicates no data on a processor
PARAMETER(st_no_data=-3,st_nd=-3)
!
! Error codes
!
INTEGER :: st_bad_array_param ! error code for dodgy array params
PARAMETER(st_bad_array_param=1)

INTEGER :: st_bad_address     ! error code for address violation
PARAMETER(st_bad_address=2)

INTEGER :: st_unknown ! error code for unknown option
PARAMETER(st_unknown=3)

INTEGER :: st_bad_wraparound ! error code for illegal wraparound
PARAMETER(st_bad_wraparound=4)

INTEGER :: st_illegal_weight ! error code for illegal weighting
PARAMETER(st_illegal_weight=9)

INTEGER :: unknown_weight ! error code for an unknown weight
PARAMETER(unknown_weight=10)

INTEGER :: unknown_mask ! error code for an unknown mask
PARAMETER(unknown_mask=11)

INTEGER :: unknown_processing ! error code for unknown processing
PARAMETER(unknown_processing=12)

INTEGER :: nonsense ! error code for general nonsense request
PARAMETER(nonsense=13)

END MODULE sterr_mod
