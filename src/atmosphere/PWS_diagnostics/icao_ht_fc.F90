! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routine to convert pressure fields to ICAO height
MODULE icao_ht_fc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = "ICAO_HT_FC_MOD"

PRIVATE

PUBLIC :: ICAOHeight_FC

CONTAINS

SUBROUTINE ICAOHeight_FC(PField,       &  ! in 
                       IHField,         &  ! out
                       numcols,numrows)    ! inout

! Description:
!   Convert pressure (Pa) to height (kft) using ICAO standard atmosphere
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE planet_constants_mod, ONLY: g_over_r
USE missing_data_mod, ONLY:  imdi, rmdi
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: numcols, numrows
REAL, INTENT(IN)    :: PField(numcols,numrows)       !P field for conversion
REAL, INTENT(INOUT)   :: IHField(numcols,numrows)     !ICAO height in kft

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "ICAOHEIGHT_FC"
REAL, PARAMETER :: ft2m = 0.3048
REAL, PARAMETER :: mtokft = 1/(ft2m * 1000.0)
REAL, PARAMETER :: Lapse_RateL = 6.5e-03       ! For levels below 11,000 gpm
REAL, PARAMETER :: Lapse_RateU = -1.0e-03      ! For levels above 11,000 gpm
REAL, PARAMETER :: Press_Bot = 101325.0        ! ICAO std: surface pressure
REAL, PARAMETER :: Press_Mid = 22632.0         !      pressure @ 11,000 gpm
REAL, PARAMETER :: Press_Top = 5474.87         !      pressure @ 20,000 gpm
REAL, PARAMETER :: Temp_Bot = 288.15           ! Surface temperature
REAL, PARAMETER :: Temp_Top = 216.65           ! Temperature of isothermal layer
REAL, PARAMETER :: Gpm1 = 11000.0              ! Ht limit (gpm) for std lower
                                               ! lapse rate
REAL, PARAMETER :: Gpm2 = 20000.0       ! Ht (gpm) of top of isothermal layer

! Local Variables:
INTEGER :: i,j      ! loop counters
REAL :: Pressure    ! Local pressure
REAL :: zp1, zp2    ! Exponents used for calculation

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! End of header --------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

zp1 = Lapse_RateL/G_over_R ! Exponents used for
zp2 = Lapse_RateU/G_over_R ! calculation


DO j = 1,NumRows
  DO i = 1,NumCols
    pressure = PField(i,j)
    IF ( (pressure <= 1000.0) .AND. (pressure >= 0.0) ) THEN
      pressure = 1000.0
    END IF
    IF ( pressure > Press_Bot) THEN
      pressure = Press_Bot
    END IF

    IF (pressure == rmdi) THEN
      IHField(i,j) = rmdi
    ELSE IF (pressure > Press_Mid) THEN ! Hts up to 11,000 GPM
      pressure = pressure/Press_Bot
      pressure = 1.0 - pressure**zp1
      IHField(i,j) = pressure*Temp_Bot/Lapse_RateL

    ELSE IF (pressure > Press_Top) THEN ! Hts between 11,000
                                        !     and     20,000 GPM
      pressure = pressure/Press_Mid
      pressure = -LOG(pressure)
      IHField(i,j) = Gpm1 + pressure*Temp_Top/G_over_R

    ELSE                                ! Hts above 20,000 GPM
      pressure = pressure/Press_Top
      pressure = 1.0 - pressure**zp2
      IHField(i,j) = Gpm2 + pressure*Temp_Top/Lapse_RateU

    END IF

  END DO
END DO

WHERE ( IHField /= rmdi )
  IHField = IHField * MtoKft
END WHERE

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE ICAOHeight_FC

END MODULE icao_ht_fc_mod
