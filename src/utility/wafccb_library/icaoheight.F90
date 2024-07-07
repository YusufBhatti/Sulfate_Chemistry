! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine to convert pressure fields to ICAO height
MODULE icaoheight_mod

IMPLICIT NONE 

CONTAINS

SUBROUTINE icaoheight(PField, rmdi)

! Description:
!   Convert pressure (Pa) to height (kft) using ICAO standard atmosphere
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: WAFC CB Library
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

IMPLICIT NONE

! P field for conversion
REAL, INTENT(INOUT)  :: PField(:, :)  
REAL, INTENT(IN)     :: rmdi

! Local Constants:
REAL, PARAMETER :: ft2m = 0.3048
REAL, PARAMETER :: mtokft = 1/(ft2m * 1000.0)
REAL, PARAMETER :: Lapse_RateL = 6.5e-03   ! For levels below 11,000 gpm
REAL, PARAMETER :: Lapse_RateU = -1.0e-03  ! For levels above 11,000 gpm
REAL, PARAMETER :: Press_Bot = 101325.0    ! ICAO std: surface pressure
REAL, PARAMETER :: Press_Mid = 22632.0     !      pressure @ 11,000 gpm
REAL, PARAMETER :: Press_Top = 5474.87     !      pressure @ 20,000 gpm
REAL, PARAMETER :: Temp_Bot = 288.15       ! Surface temperature
REAL, PARAMETER :: Temp_Top = 216.65       ! Temperature of isotherm
REAL, PARAMETER :: Gpm1 = 11000.0          ! Ht limit (gpm) for std lower
                                           ! lapse rate
REAL, PARAMETER :: Gpm2 = 20000.0          ! Ht (gpm) of top of isothermal layer
REAL, PARAMETER :: g_over_r = 9.80665 / 287.05 ! g over earth radius

! Local Variables:
INTEGER :: i,j      ! loop counters
REAL    :: Pressure    ! Local pressure
REAL    :: zp1, zp2    ! Exponents used for calculation

! End of header --------------------------------------------------------

zp1 = Lapse_RateL/G_over_R ! Exponents used for
zp2 = Lapse_RateU/G_over_R ! calculation

DO j = 1, SIZE(PField, 2)
  DO i = 1, SIZE(PField, 1)
    pressure = PField(i,j)
    IF (pressure == rmdi) THEN
      CYCLE
    END IF
    IF ( (pressure <= 1000) .AND. (pressure >= 0.0) ) THEN
      pressure = 1000.0
    END IF
    IF ( pressure > Press_Bot) THEN
      pressure = Press_Bot
    END IF

    IF (pressure > Press_Mid) THEN ! Hts up to 11,000 GPM
      pressure = pressure/Press_Bot
      pressure = 1.0 - pressure**zp1
      PField(i,j) = pressure*Temp_Bot/Lapse_RateL

    ELSE IF (pressure > Press_Top) THEN ! Hts between 11,000
                                        !     and     20,000 GPM
      pressure = pressure/Press_Mid
      pressure = -LOG(pressure)
      PField(i,j) = Gpm1 + pressure*Temp_Top/G_over_R

    ELSE                                ! Hts above 20,000 GPM
      pressure = pressure/Press_Top
      pressure = 1.0 - pressure**zp2
      PField(i,j) = Gpm2 + pressure*Temp_Top/Lapse_RateU

    END IF

  END DO
END DO

WHERE ( PField /= rmdi )
  PField = PField * MtoKft
END WHERE

END SUBROUTINE icaoheight

END MODULE icaoheight_mod
