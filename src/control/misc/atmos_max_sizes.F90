! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc.
!
! Description
!   This module provides parameters giving the maximum likely sizes
!   of key UM resolution variables, useful for sizing static arrays.

MODULE atmos_max_sizes

IMPLICIT NONE

! Maximum sector size for I/O

#if defined(UTILIO)
! Small execs needs a larger value to cope with some of the wet
! model types
INTEGER, PARAMETER :: row_length_max   = 8000 ! maximum row length
INTEGER, PARAMETER :: rows_max         = 6000 ! max no of rows
#else
INTEGER, PARAMETER :: row_length_max   = 4096 ! maximum row length
INTEGER, PARAMETER :: rows_max         = 3200 ! max no of rows
#endif

! maximum permitted size of a halo
INTEGER, PARAMETER :: Max_Halo_Size    = 20
INTEGER, PARAMETER :: model_levels_max = 250 ! max no of total levels


INTEGER, PARAMETER :: horiz_dim_max=MAX(row_length_max,rows_max)
INTEGER, PARAMETER :: Max2DFieldSize   = row_length_max*rows_max
INTEGER, PARAMETER :: MaxHaloArea      = horiz_dim_max*Max_Halo_Size

END MODULE atmos_max_sizes
