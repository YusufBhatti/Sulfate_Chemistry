! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
MODULE o3intp_mod

IMPLICIT NONE
!  ---------------------------------------------------------------------
!  Module to specify allowed methods of interpolating from the
!  ancillary file.
!- ---------------------------------------------------------------------
!
INTEGER, PARAMETER :: io3_3dspec = 1
!       Ozone is provided as a full 3D field.
INTEGER, PARAMETER :: io3_2dspec = 2
!       Ozone is expanded from a 2D field by direct copying.
INTEGER, PARAMETER :: io3_2dmasscon = 3
!       Ozone is expanded from a 2D field with conservation of mass.
INTEGER, PARAMETER :: io3_trop_map = 4
!       Ozone mixing ratios are set by mapping each height at each
!       grid-point in the real profile to a height in the ancillary
!       profile and using the mixing ratio there. The mapping is
!       set using the height of the tropopause
INTEGER, PARAMETER :: io3_trop_map_masscon = 5
!       Ozone mixing ratios are set as above, but scaled so as to
!       preserve the vertically integrated column ozone
!
! ----------------------------------------------------------------------

END MODULE o3intp_mod
