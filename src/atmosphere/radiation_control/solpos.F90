! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE solpos_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLPOS_MOD'
CONTAINS

SUBROUTINE solpos (day, year, seconds, timestep, l_planet_obs,    &
                   eqt, sindec, scs, sindec_obs, eqt_obs)

!   Subroutine SOLPOS   ----------------------------------------------
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Radiation Control

!   Purpose :
!    Calculations of the earth's orbit described in the first page of
!    the "Calculation of incoming insolation" section of UMDP 23, i.e.
!    from the day of the year (and, in forecast mode, whether it is a
!    leap year) and the orbital "constants" (which vary over
!    "Milankovitch" timescales) it calculates the sin of the solar
!    declination and the inverse-square scaling factor for the solar
!    "constant".  It is thus intrinsically scalar.  The FORTRAN code
!    present depends on whether *DEF CAL360 is set during UPDATE: this
!    replaces the Julian calendar with the climate-mode 360-day calendar

!   ------------------------------------------------------------------


USE conversions_mod, ONLY: pi, rsec_per_day
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE nlstcall_mod, ONLY: lcal360
USE rad_input_mod, ONLY: l_sec_var
USE planet_constants_mod, ONLY: l_planet_orbit,                   &
  planet_e, planet_lph, planet_oblq, planet_a, planet_M,          &
  planet_de, planet_dlph, planet_doblq, planet_da, planet_dM,     &
  planet_ha, planet_dha, planet_epoch, i_eqt,                     &
  ip_smart, ip_mueller, planet_obs_lon, planet_obs_lat
USE umPrintMgr, ONLY: umPrint

USE orbprm_mod, ONLY: orbprm
IMPLICIT NONE

INTEGER, INTENT(IN) :: day        ! Day-number in the year
INTEGER, INTENT(IN) :: year       ! Calendar year
REAL, INTENT(IN)    :: seconds    ! Seconds since midnight GMT
REAL, INTENT(IN)    :: timestep   ! Length of timestep in seconds

LOGICAL, INTENT(IN) :: l_planet_obs
!       Calculate angles towards distant observer

REAL, INTENT(OUT)   :: sindec
!       Sin(solar declination)
REAL, INTENT(OUT)   :: scs
!       Solar constant scaling factor
REAL, INTENT(OUT)   :: eqt
!       The equation of time, specified as an hour angle in radians.
REAL, INTENT(OUT)   :: eqt_obs
!       The equation of time, specified as an hour angle in radians for the
!       position of the observer rather than the sun.
REAL, INTENT(OUT)   :: sindec_obs
!       Sin(observer declination)

! Mathematical constants:
REAL, PARAMETER :: twopi = 2.0 * pi

! Parameters of the Earth's orbit:
REAL :: e       ! Eccentricity of the orbit
REAL :: gamph   ! Supplement of the longitude of the perihelion
REAL :: oblq    ! Obliquity of the orbit

! Derived orbital constants:
REAL :: e1, e2, e3, e4, y, y2, p

REAL :: m       ! Mean anomaly: positional angle of a "mean" Earth
                !  rotating around the sun with a constant angular
                !  speed equal to 2pi/T and counted counterclock-
                !  wise from the perihelion
REAL :: v       ! True anomaly: positional angle of Earth in its
                !  orbit, counted counterclockwise from perihelion

REAL :: ha
!       Hour angle of the planet (in radians) at Earth midnight

REAL :: day_number_at_midnight
!       Number of days from epoch to Earth midnight at the beginning of
!       the current Earth day.
REAL :: day_number
!       Number of days from epoch to the middle of the current timestep

REAL, PARAMETER :: jd1 = 1721425.5
!       Julian day number for the beginning of the first day (year 1) of
!       the Gregorian calendar.

! Cartesian coordinates of distant observer in orbital and equatorial
! reference frames:
REAL :: obs_orb_x, obs_orb_y, obs_orb_z
REAL :: obs_eqt_x, obs_eqt_y, obs_eqt_z

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=50000)          :: lineBuffer
CHARACTER(LEN=*), PARAMETER   :: RoutineName='SOLPOS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_planet_orbit) THEN

  ! Calculate number of days from epoch to the middle of the 
  ! current timestep. The conversion to Julian day number assumes
  ! the Gregorian calendar for the calculation of leap years.
  day_number_at_midnight = jd1 + REAL( 365*(year-1)    &
    + (year-1)/4 - (year-1)/100 + (year-1)/400 )        &
    + REAL(day) - 1.0 - planet_epoch
  day_number = day_number_at_midnight                   &
    + (seconds+timestep/2.0)/rsec_per_day

  m    = planet_M    + planet_dM    * day_number
  e    = planet_e    + planet_de    * day_number
  oblq = planet_oblq + planet_doblq * day_number
  gamph = pi - (planet_lph + planet_dlph * day_number)

  ha = planet_ha + planet_dha * day_number_at_midnight

ELSE

  ! Determine the orbital parameters for the Earth
  CALL orbprm(day, year, l_sec_var, lcal360, e, gamph, oblq, m)

END IF

! Calculate the coefficients in the equation of the centre and
! thence derive the true anomaly.
e1 = e * ( 2.0 - .25 * e*e )
e2 = 1.25 * e*e
e3 = e*e*e * 13.0 / 12.0

! True anomaly, equation 87 in Smart on p. 120 (UMDP23 Eq 3.1.2)
v  = m + e1*SIN(m) + e2*SIN(2.0*m) + e3*SIN(3.0*m)

! Solar constant scaling factor (UMDP23 Eq 3.1.4)
e4  = ( (1.0 + e*e*.5) / (1.0 - e*e) )**2
scs = e4 * ( 1.0 + e * COS(v) ) **2
IF (l_planet_orbit) THEN
  scs = scs / (planet_a + planet_da*day_number)**2
END IF

! sin(solar declination) (UMDP23 Eq 3.1.5)
! The solar declination is related to
!  the true longitude of the earth (lambda) by:
!  sindec = sin(obliquity) * sin(lambda)
! Lambda is counted counterclockwise from the vernal equinox
!  and is related to v (the true anomaly) through
!  lambda = v + (longitude of perihelion)
sindec = SIN(oblq) * SIN (v - gamph)


! Calculate the equation of time.
SELECT CASE (i_eqt)
CASE (ip_smart)
  ! Use equation (29)on page 149 of Smart (1944).
  ! (Recall the factor of 5/4 in the definition of e2).
  y   = ( TAN ( 0.5*oblq ) )**2
  eqt = y * SIN(2.0 * ( m - gamph ))                &
        - 2.0*e * SIN(m)                            &
        + 4.0*e*y * SIN(m) * COS(2.0*( m - gamph )) &
        - 0.5*y*y * SIN(4.0*( m - gamph ))          &
        - e2 * SIN(2.0*m)
CASE (ip_mueller)
  ! M. Mueller, Acta Physica Polonica A 88 Supplement, S-49 (1995)
  y   = ( TAN ( 0.5*oblq ) )**2
  y2  = y*y
  e2  = e*e 
  p   = 0.5*pi - gamph
  eqt = - y * (1.0-4.0*e2) * SIN(2.0*(m+p)) &
        - 2.0  * e * SIN(m)                 &
        + 2.0  * e * y  * SIN(m+2.0*p)      &
        - 2.0  * e * y  * SIN(3.0*m+2.0*p)  &
        - 0.5  * y2 * SIN(4.0*(m+p))        &
        - 1.25 * e2 * SIN(2.0*m)            &
        + 2.0  * e * y2 * SIN(3.0*m+4.0*p)  &
        - 2.0  * e * y2 * SIN(5.0*m+4.0*p)  &
        - 3.25 * e2* y  * SIN(4.0*m+2.0*p)  &
        - y2   * y * SIN(6.0*(m+p))/3.0
CASE DEFAULT
  eqt=0.0e+00
END SELECT

IF (l_planet_orbit) THEN
  ! Add on a correction for the planets hour angle at Earth midnight
  eqt = eqt + MODULO(ha + pi, twopi)
END IF


IF (l_planet_obs) THEN
  ! Calculations needed to output emission spectra towards a distant observer
  WRITE(lineBuffer,'(A,F18.14)') 'Phase angle of orbit =', &
    MODULO(v - gamph + pi - planet_obs_lon, twopi)
  CALL umPrint(lineBuffer,src='solpos')

  ! The observer declination is calculated using coordinate transformation
  ! from the planet's orbital to the planet's equatorial reference frame:
  obs_orb_x = COS(planet_obs_lat)*COS(planet_obs_lon)
  obs_orb_y = COS(planet_obs_lat)*SIN(planet_obs_lon)
  obs_orb_z = SIN(planet_obs_lat)
  obs_eqt_x = obs_orb_x
  obs_eqt_y = COS(oblq)*obs_orb_y - SIN(oblq)*obs_orb_z
  obs_eqt_z = SIN(oblq)*obs_orb_y + COS(oblq)*obs_orb_z
  sindec_obs = obs_eqt_z
  
  ! We add the position of the observer on to the equation of time
  ! noting that the equation of time is in a clockwise sense.
  ! First add the direction of the vernal equinox (the zero point for
  ! the longitude of the observer in the planet's equatorial frame):
  eqt_obs = eqt + ATAN2( COS(oblq) * SIN(v - gamph), COS(v - gamph) )
  ! Then subtract the longitude of the observer in the planet's equatorial
  ! reference frame (reduces to planet_obs_lon if oblq=0):
  eqt_obs = eqt_obs - ATAN2(obs_eqt_y, obs_eqt_x)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE solpos
END MODULE solpos_mod
