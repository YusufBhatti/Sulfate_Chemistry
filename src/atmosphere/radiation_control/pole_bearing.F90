! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE pole_bearing_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'POLE_BEARING_MOD'
CONTAINS

SUBROUTINE pole_bearing(row_length, rows, lat_rot_NP, &
    long_rot_NP, bear_rot_NP)

USE trignometric_mod, ONLY: true_latitude, true_longitude
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Description:
!   Calculates bearing of Grid North (from True North) for gridpoints
!   in a Local Area Model rotated grid.
!
! Method:
!   Spherical trig.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.

! Subroutine arguments

INTEGER, INTENT(IN) :: &
     row_length, rows         ! grid size
REAL, INTENT(IN) ::  &
     lat_rot_NP,     &        ! Real latitude and longitude
     long_rot_NP              !  of 'pseudo' N pole in radians.

REAL, INTENT(OUT) :: &
     bear_rot_NP(row_length,rows)  ! Bearing of 'pseudo' N pole (rads)

! Local variables
REAL :: long_diff(row_length,rows)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='POLE_BEARING'

! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate bearing of Grid North
long_diff = true_longitude - long_rot_NP

bear_rot_NP = ATAN2( -COS(lat_rot_NP)*SIN(long_diff),  &
      COS(true_latitude)*SIN(lat_rot_NP) -                  &
      SIN(true_latitude)*COS(lat_rot_NP)*COS(long_diff) )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE pole_bearing
END MODULE pole_bearing_mod
