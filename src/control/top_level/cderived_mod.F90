! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top_level

MODULE cderived_mod

IMPLICIT NONE

! Description:
!   This file contains declarations for derived constants and
!   grid definition coords in radians within the atmospheric model.

!   The derived constants are calculated in the routine SETCONA.
!
      ! No of cloud types ie low/med/high
INTEGER, PARAMETER :: NUM_CLOUD_Types = 3

! derived constants:
INTEGER :: LOW_BOT_Level      ! Bottom level of lowest cloud type
INTEGER :: LOW_TOP_Level      ! Top      "    "   "       "    "
INTEGER :: MED_BOT_Level      ! Bottom   "    "  med      "    "
INTEGER :: MED_TOP_Level      ! Top      "    "   "       "    "
INTEGER :: HIGH_BOT_Level     ! Bottom   "    "  top      "    "
INTEGER :: HIGH_TOP_Level     ! Top      "    "   "       "    "

! height values to split model levels into l/m/h cloud
REAL ::    h_split(NUM_CLOUD_Types+1)

LOGICAL :: elf                ! T if atmosphere model on LAM grid

! Grid definition co-ordinates in radians
REAL:: Delta_lambda       ! EW (x) grid spacing in radians
REAL:: Delta_phi          ! NS (y) grid spacing in radians
REAL:: Base_phi           ! Lat of first theta point in radians
REAL:: Base_lambda        ! Long of first theta point in radians
REAL:: lat_rot_NP         ! Real lat of 'pseudo' N pole in radians
REAL:: long_rot_NP        ! Real long of 'pseudo' N pole in radians

CONTAINS

SUBROUTINE h_split_defaults()

! height values to split model levels into l/m/h cloud

IMPLICIT NONE

! Default settings for h_split
!  low:middle:high cloud model levels =(1)->(2):(2)->(3):(3)->(4)
h_split(1) =   111.0        ! ICAO 1000mb height (m)
h_split(2) =  1949.0        ! ICAO  800mb height (m)
h_split(3) =  5574.0        ! ICAO  500mb height (m)
h_split(4) = 13608.0        ! ICAO  150mb height (m)

END SUBROUTINE h_split_defaults

END MODULE cderived_mod
