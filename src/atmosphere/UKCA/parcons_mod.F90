! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
MODULE parcons_mod

IMPLICIT NONE
!     *PARAMETER* OF GLOBAL CONSTANTS.
!
REAL, PARAMETER :: g = 9.806
REAL, PARAMETER :: pi = 3.14159265358978
REAL, PARAMETER :: circ = 40000000.0
REAL, PARAMETER :: zpi  = 2.0*pi
REAL, PARAMETER :: rad  = pi/180.0
REAL, PARAMETER :: deg  = 180.0/pi
REAL, PARAMETER :: r    = circ/zpi
!
!      VARIABLE.   TYPE.     PURPOSE.
!     ---------   -------   --------
!     *G*         REAL      ACCELERATION OF GRAVITY.
!     *PI*        REAL      PI.
!     *CIRC*      REAL      EARTH CIRCUMFERENCE (METRES).
!     *RAD*       REAL      PI / 180.
!     *DEG*       REAL      180. / PI.
!     *ZPI*       REAL      2. * PI.
!     *R*         REAL      EARTH RADIUS        (METRES).
!

END MODULE parcons_mod
