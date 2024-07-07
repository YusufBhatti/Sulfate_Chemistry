! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
!  Holds information on pp header i.e. dates and grid

MODULE crmstyle_pp_data_mod

IMPLICIT NONE
SAVE

! Description:
!   Module containing information on pp headers for output
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Crmstyle_coarse_grid
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3  programming standards.
!
! Declarations:

INTEGER ::      &
  iyear         & ! year   validity time from input fields
 ,imon          & ! month
 ,iday          & ! day
 ,ihour         & ! hour
 ,imin          & ! min
 ,isec            ! second

INTEGER ::       &
  isyear         & ! year    data time
 ,ismon          & ! month
 ,isday          & ! day
 ,ishour         & ! hour
 ,ismin          & ! min
 ,issec            ! second

INTEGER ::       &
  itim             ! time indicator includes calendar info

REAL                :: &
  bzy                  &  ! Full input grid info
 ,bzx                  &
 ,bdy                  &
 ,bdx                  &
 ,pseudo_lat           &
 ,pseudo_lon

! new output grid info
REAL                :: &
  new_bzy              &
 ,new_bzx              &
 ,new_bdy              &
 ,new_bdx              &
 ,origin_y1            &
 ,origin_x1

! high resolution output grid info - only setup if l_bcu_mask is .TRUE.
REAL                :: &
  mask_bzy             &
 ,mask_bzx             &
 ,mask_bdy             &
 ,mask_bdx

END MODULE crmstyle_pp_data_mod
