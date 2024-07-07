! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

PROGRAM main

!    Purpose:  Calculates latitude and longitude on standard grid
!              from input arrays of latitude and longitude on
!              equatorial latitude-longitude (eq) grid used
!              in regional models. Both input and output latitudes
!              and longitudes are in degrees.


IMPLICIT NONE

REAL :: plon,plat,rlon,rlat

READ (5,*)plon,plat,rlon,rlat

CALL EQTOLL(RLAT,RLON,RLAT,RLON,PLAT,PLON,1)

WRITE(6,'(2(F8.3))')RLON,RLAT
STOP
END
