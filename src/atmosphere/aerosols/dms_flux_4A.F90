! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE dms_flux_mod_4A

IMPLICIT NONE

! Indicators for sea-air flux calculation methods
INTEGER, PARAMETER :: i_liss_merlivat = 1
INTEGER, PARAMETER :: i_wanninkhof = 2
INTEGER, PARAMETER :: i_nightingale = 3
INTEGER, PARAMETER :: i_blomquist = 4
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DMS_FLUX_MOD_4A'

CONTAINS

SUBROUTINE dms_flux_4A(                                            &
  ! Arguments IN
  row_length, rows,                                                &
  wind_10m, tstar, land_fract, dms_conc,                           &
  i_dms_flux,  &
  ! Arguments OUT
  f_dms )
!---------------------------------------------------------------------
! Purpose: To calculate the flux of DMS (as kg m-2 s-1 of sulphur)
!          from the ocean surface as a function of its concentration
!          in seawater and of windspeed. The sea-air exchange can
!          be determined according to one of three commonly-used
!          parametrization schemes, those of Liss & Merlivat (1986),
!          Wanninkhof (1992) or Nightingale et al. (2000). The routine
!          is called by Aero_Ctl.

! Method:  The Schmidt number
!          for DMS is calculated as in Saltzman et al. (1993), and
!          used with the windspeed to determine the mass transfer (or
!          "piston") velocity according to the desired parametrization.
!          This is then used to determine the sea-air mass flux of DMS
!          as a function of sea-water DMS concentration. High surface
!          temperatures (caused by the land portion of a gridbox when
!          coastal tiling is not active) cause negative Sc values which
!          would give a floating-point error in the k_DMS calculation,
!          so the Tstar values are capped. This shouldn't be a problem
!          when coastal tiling is on as then the Tstar values passed in
!          are those for sea only.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols

! Code Description:
!  Language: Fortran 90
!  This code is written to UMDP3 v8 programming standards

!---------------------------------------------------------------------

USE conversions_mod, ONLY: zerodegc
USE water_constants_mod, ONLY: tfs

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE umPrintMgr, ONLY: umPrint, umMessage
USE ereport_mod, ONLY: EReport

IMPLICIT NONE

! Arguments with intent IN:
INTEGER :: row_length
INTEGER :: rows

REAL    :: wind_10m(row_length, rows)   ! 10m windspeed (ms-1)
REAL    :: tstar(row_length, rows)      ! Surface temperature (K)
REAL    :: land_fract(row_length, rows) ! Fraction of land in gridbox
REAL    :: dms_conc(row_length, rows)   ! Concentration of DMS in seawater (nmol l-1)

! Switch to determine which scheme to use to calculate mass transfer velocity
INTEGER :: i_dms_flux

! Arguments with intent OUT:

REAL    :: f_dms(row_length, rows)       ! Sea-air flux of DMS (kg[S] m-2 s-1)

! Local variables:
INTEGER  :: i,j                          ! Loop counters
REAL     :: sc(row_length, rows)         ! Schmidt number
REAL     :: k_dms(row_length, rows)      ! Piston velocity of DMS (cm h-1)
REAL     :: t_c                          ! Surface temperature in degrees Celsius
! Piston velocities for gases with Schmidt numbers of 600 & 660 resp. (cm h-1)
REAL     :: k_600
REAL     :: k_660
REAL     :: n                            ! Schmidt number exponent
REAL, PARAMETER :: t_max=47.0            ! Max T to avoid breaking the Sc fit (C)

INTEGER  :: errcode

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DMS_FLUX_4A'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Calculate the Schmidt number (Sc):

DO j = 1, rows
  DO i = 1, row_length
    t_c = MIN((MAX(tstar(i, j), tfs) - zerodegc), t_max)
!!  Saltzmann et al 1993   YAB 07/11/22 (918 t_c=20)
!    sc(i, j) = 2674.0 - (147.12*t_c) + (3.726*t_c**2)                          &
!      - (0.038*t_c**3)

!!  Updated to Wannkinkhof (2014)   YAB UCNZ 07/11/22 (940.6088 t_c=20)
    sc(i, j) = 2855.7+  (-177.63 + (6.0438 + (-0.11645 + 0.00094743*t_c)*t_c)*t_c)*t_c
  END DO
END DO

! Determine the mass transfer (or "piston") velocity (k_DMS) over sea
! according to the specified parametrization scheme:

SELECT CASE (i_dms_flux)

CASE (i_liss_merlivat)
  DO j = 1, rows
    DO i = 1, row_length
      IF (wind_10m(i, j)  <=  3.6) THEN
        k_600 = 0.17 * wind_10m(i, j)
        n = -2.0/3.0
      END IF
      IF (wind_10m(i, j)  >   3.6 .AND.                                        &
        wind_10m(i, j)  <=  13.0) THEN
        k_600 = (2.85*wind_10m(i, j)) - 9.65
        n = -0.5
      END IF
      IF (wind_10m(i, j)  >   13.0) THEN
        k_600 = (5.9*wind_10m(i, j)) - 49.3
        n = -0.5
      END IF
      IF (land_fract(i, j)  <   1.0) THEN
        k_dms(i, j) = k_600 * (sc(i, j)/600.0)**n
      ELSE
        k_dms(i, j) = 0.0
      END IF
    END DO
  END DO

CASE (i_wanninkhof)
  DO j = 1, rows
    DO i = 1, row_length
      k_660 = 0.251 * wind_10m(i, j)**2
      n = -0.5
      IF (land_fract(i, j)  <   1.0) THEN
        k_dms(i, j) = k_660 * (sc(i, j)/660.0)**n
      ELSE
        k_dms(i, j) = 0.0
      END IF
    END DO
  END DO

CASE (i_nightingale)
  DO j = 1, rows
    DO i = 1, row_length
      k_600 = (0.222*wind_10m(i, j)**2) + (0.333*wind_10m(i, j))
      n = -0.5
      IF (land_fract(i, j)  <   1.0) THEN
        k_dms(i, j) = k_600 * (sc(i, j)/600.0)**n
      ELSE
        k_dms(i, j) = 0.0
      END IF
    END DO
  END DO

CASE (i_blomquist)
  DO j = 1, rows
    DO i = 1, row_length
      k_600 = 0.7432 * (wind_10m(i, j)**1.33)
      n = -0.5
      IF (land_fract(i, j)  <   1.0) THEN
        k_dms(i, j) = k_600 * (sc(i, j)/660.0)**n
      ELSE
        k_dms(i, j) = 0.0
      END IF
    END DO
  END DO
 

CASE DEFAULT
  WRITE(umMessage,'(A,I4,A)') 'DMS_FLUX_4A, i_dms_flux=',i_dms_flux,  &
                      ' is not a valid choice.'
  CALL umPrint(umMessage,src='dms_flux_4a')
  errcode = 1
  CALL EReport('dms_flux_4a', errcode, 'i_dms_flux not valid')

END SELECT

! Finally, calculate the sea-air flux of DMS as a function of k_DMS
! and dissolved DMS concentration. The former requires a conversion
! from cm hour-1 to ms-1, and the latter from nanomoles per litre to
! kg[S] m-3, to return the flux in kg[S] m-2 sec-1.

DO j = 1, rows
  DO i = 1, row_length
    f_dms(i, j) = (k_dms(i, j) / 3.6e5)                                        &
      * (dms_conc(i, j) * 32.0e-9)
  END DO
END DO

!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE dms_flux_4A

END MODULE
