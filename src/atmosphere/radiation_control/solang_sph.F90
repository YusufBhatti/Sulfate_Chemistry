! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE solang_sph_mod
IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLANG_SPH_MOD'
CONTAINS

! Calculation of path of solar beam using spherical geometry
!
! Purpose :
!  Calculations of the earth's orbit described in the second page of
!  the "Calculation of incoming insolation" section of UMDP 23, i.e.
!  from the sin of the solar  declination, the position of each point
!  and the time limits it calculates how much sunlight, if any, it
!  receives.
!    Here the calculation is done for each layer of the atmosphere,
!  determining whether each layer is lit assuming spherical geometry
!  for the path of the incoming solar beam.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
SUBROUTINE solang_sph(sindec, t, dt, eqt, lat, longit, n_profile, n_layer,     &
     r_layer_in, r_level, slope_aspect, slope_angle,                           &
     horiz_ang, horiz_aspect, n_horiz_layer, n_horiz_ang, l_orog, l_skyview,   &
     lit, cosz, sol_azimuth)

USE conversions_mod, ONLY: pi
USE planet_constants_mod, ONLY: s2r, planet_dha,                               &
  l_fix_solang, solar_zenith_angle, solar_azimuth_angle, planet_radius
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE


INTEGER, INTENT(IN) :: n_profile, n_layer, n_horiz_layer, n_horiz_ang

REAL, INTENT(IN) ::                                                            &
  sindec,                                                                      &
!   Sin(solar declination)
  t, dt,                                                                       &
!   Start time (GMT) & timestep
  lat(n_profile), longit(n_profile),                                           &
!   latitude & longitude of each point
  eqt,                                                                         &
!   The equation of time (the difference between true and mean solar time)
  r_layer_in(n_profile, n_layer),                                              &
!   Height of layers from centre of planet
  r_level(n_profile, 0:n_layer),                                               &
!   Height of levels (boundaries between layers) from centre of planet
  slope_aspect(n_profile),                                                     &
!   The direction faced by the mean surface slope - i.e. the
!   bearing of the slope normal projected on the surface
!   measured in radians clockwise from true north.
  slope_angle(n_profile),                                                      &
!   Angle of the mean surface slope measured in radians from the horizontal.
  horiz_ang(n_profile, n_horiz_layer, n_horiz_ang),                            &
!   Angle in radians measured from the local zenith to the obscuring terrain.
  horiz_aspect(n_profile, n_horiz_ang)
!   The local bearing for each horizon angle measured in radians clockwise
!   from true north.

LOGICAL, INTENT(IN) :: l_orog, l_skyview

REAL, INTENT(OUT) ::                                                           &
  lit(n_profile, 0:n_layer+1),                                                 &
!   Sunlit fraction of the timestep
  cosz(n_profile, 0:n_layer+1),                                                &
!   Mean cos(solar zenith angle) during the sunlit fraction
  sol_azimuth(n_profile)
!   Mean solar azimuth angle (radians clockwise from grid north) during the
!   sunlit fraction

INTEGER :: i, j, k
!   Loop integers

REAL ::                                                                        &
  mean_omega(n_profile),                                                       &
!   Mean hour angle over the timestep measured in radians west of local noon
  sinsin(n_profile), coscos(n_profile),                                        &
!   Products of the sines and cosines of solar declination and latitude
  omega_sunrise, omega_sunset,                                                 &
!   Hour-angle of sunrise and sunset
  omega_90_rising, omega_90_setting,                                           &
!   Hour-angle of rising and setting sun at zenith angle of 90 degrees
  hat(n_profile),                                                              &
!   Local hour angle at the start time.
  omegab, omegae, omegas,                                                      &
!   Local hour angle at the beginning and end of the timestep and sunset.
!   All measured in radians after local sunrise, not from local noon as the
!   true hour angle is.
  omega1, omega2, omega3, omega4,                                              &
!   Local hour angle at the boundaries over which cosz is integrated.
  difsin, diftim,                                                              &
!   A difference-of-sines intermediate value and the corresponding time period
  trad, dtrad,                                                                 &
!   These are the start-time and length of the timestep (T & DT)
!   converted to radians after midday GMT, or equivalently, hour
!   angle of the mean sun on the Greenwich meridian.
  sinlat(n_profile), dec,                                                      &
!   Working variables for azimuth calculation.
  cos_horizon_ang,                                                             &
!   Cosine of the horizon for a given level, measured from the zenith
  cos_omega_horizon,                                                           &
!   Cosine of the hour angle when the sun is on the horizon for a given level
  cos_omega_90,                                                                &
!   Cosine of the hour angle when the sun is perpendicular to the zenith
  impact(n_profile, 0:n_layer+1),                                              &
!   Impact parameter: height at which beam is tangential to the surface
  r_layer(n_profile, 0:n_layer+1),                                             &
!   Height of layers (plus surface and TOA) from centre of planet
  rel_horiz_aspect(n_profile, n_horiz_ang)
!   Bearing of horizon angles relative to the solar azimuth

INTEGER :: min_sub(n_profile), max_sub(n_profile)
!   Subscripts for horizon angles bracketing the solar azimuth

REAL, PARAMETER :: twopi = 2.0*pi
REAL, PARAMETER :: eps = EPSILON(1.0)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLANG_SPH'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Incorporate height of surface and TOA into r_layer so things can be dealt
! with in a single loop
r_layer(:,0) = r_level(:,0)
r_layer(:,1:n_layer) = r_layer_in
r_layer(:,n_layer+1) = r_level(:,n_layer)

IF (ABS(planet_dha) < eps) THEN

  ! Planet is tidally locked.
  ! The longitude of the overhead sun may vary with the equation of time
  ! to allow for an eccentric orbit and a constant rotation rate.
  mean_omega = MODULO(longit + eqt - pi, twopi)

  ! The solar declination is allowed to vary with the orbital parameters.
  dec = ASIN(sindec)
  sol_azimuth = MODULO(ATAN2( -COS(dec)*SIN(mean_omega),                       &
    COS(lat)*sindec - SIN(lat)*COS(dec)*COS(mean_omega) ), twopi)
  DO k=0, n_layer+1
    DO j = 1, n_profile
      cosz(j,k) = sindec*SIN(lat(j)) + COS(dec)*COS(lat(j))*COS(mean_omega(j))
      ! Impact parameter
      impact(j,k)=r_layer(j,k)*SQRT(1.0-cosz(j,k)**2)
      IF (cosz(j,k) >= eps .OR. impact(j,k) > planet_radius) THEN
        ! Sun is above the horizon (planet_radius for now).
        lit(j,k) = 1.0
      ELSE
        lit(j,k) = 0.0
      END IF
    END DO
  END DO

ELSE IF (l_fix_solang) THEN

  ! Fix the sun at a particular zenith and azimuth angle for all points
  ! for idealised tests.
  sol_azimuth = solar_azimuth_angle
  DO k=0, n_layer+1
    DO j = 1, n_profile
      cosz(j,k) = COS(solar_zenith_angle)
      ! Impact parameter
      impact(j,k)=r_layer(j,k)*SQRT(1.0-cosz(j,k)**2)
      IF (cosz(j,k) >= eps .OR. impact(j,k) > planet_radius) THEN
        ! Sun is above the horizon (planet_radius for now).
        lit(j,k) = 1.0
      ELSE
        lit(j,k) = 0.0
      END IF
    END DO
  END DO

ELSE

  trad = t * s2r - pi
  dtrad = dt * s2r
  dec = ASIN(sindec)
  DO j = 1, n_profile
    sinlat(j) = SIN(lat(j))
    sinsin(j) = sindec * sinlat(j)
    coscos(j) = SQRT( (1.0-sindec**2) * (1.0-sinlat(j)**2) )

    ! Since T and DT represent mean time and all solar calculations
    ! are done using the variable HAT, it suffices to add the
    ! equation of time on to HAT.
    hat(j) = longit(j) + trad + eqt

    ! Calculate the solar azimuth angle here to allow horizon angles to
    ! be calculated for each layer incorporating the orography (to be done).
    mean_omega(j) = hat(j) + dtrad/2.0
    sol_azimuth(j) = MODULO(ATAN2( -COS(dec)*SIN(mean_omega(j)),               &
      COS(lat(j))*sindec - sinlat(j)*COS(dec)*COS(mean_omega(j)) ), twopi)
  END DO
  IF (l_skyview) THEN
    ! Calculate cosine of terrain horizon angle
    rel_horiz_aspect = MODULO( horiz_aspect - &
      SPREAD(sol_azimuth,2,n_horiz_ang), pi*2.0)
    min_sub=MINLOC(rel_horiz_aspect, dim=2)
    max_sub=MAXLOC(rel_horiz_aspect, dim=2)
  END IF
  DO k=0, n_layer
    DO j = 1, n_profile
      IF (l_skyview .AND. n_horiz_layer > k) THEN
        cos_horizon_ang = COS( (horiz_ang(j, k+1, min_sub(j))                  &
        *(pi*2.0 - rel_horiz_aspect(j, max_sub(j)))                            &
        +horiz_ang(j, k+1, max_sub(j))*rel_horiz_aspect(j, min_sub(j)))        &
        /(pi*2.0 - rel_horiz_aspect(j, max_sub(j))                             &
        +rel_horiz_aspect(j, min_sub(j))) )
      ELSE
        cos_horizon_ang = -SQRT(MAX(1.0-planet_radius**2/r_layer(j,k)**2,0.0))
      END IF
      cos_omega_horizon = (cos_horizon_ang - sinsin(j)) / coscos(j)
      IF ( cos_omega_horizon > 1.0 ) THEN
        ! Perpetual night
        diftim = 0.0
      ELSE
        IF ( cos_omega_horizon < -1.0 ) THEN
          ! Perpetual day: hour angles are start and end of timestep
          omega1 = hat(j)
          omega2 = hat(j) + dtrad
          difsin = SIN(omega2) - SIN(omega1)
          diftim = dtrad
        ELSE
          ! At this latitude some points are sunlit, some not. Different ones
          ! need different treatment.
          ! The logic seems simplest if one takes all "times" - actually hour
          ! angles - relative to sunrise, but they must be kept in the
          ! range 0 to 2pi for the tests on their orders to work.
          omega_sunset  = ACOS(cos_omega_horizon)
          omega_sunrise = -omega_sunset
          omegab = MODULO(hat(j) - omega_sunrise, twopi)
          omegae = MODULO(omegab + dtrad, twopi)
          omegas = omega_sunset - omega_sunrise
          IF (omegab <= omegas .OR. omegab < omegae) THEN
            omega1 = omegab + omega_sunrise
          ELSE
            ! Sun rises during timestep
            omega1 = omega_sunrise
          END IF
          IF (omegae <= omegas) THEN
            omega2 = omegae + omega_sunrise
          ELSE
            ! Sun sets during timestep
            omega2 = omega_sunset
          END IF
          IF (omegae > omegab .AND. omegab > omegas) THEN
            ! Sun does not rise during timestep
            diftim = 0.0
          ELSE IF (omega1 > omega2) THEN
            ! Sun sets and then rises again within the timestep.
            difsin = SIN(omega_sunset) - SIN(omega1) + &
                     SIN(omega2) - SIN(omega_sunrise) 
            diftim = omega_sunset - omega1 + omega2 - omega_sunrise
          ELSE
            difsin = SIN(omega2) - SIN(omega1)
            diftim = omega2 - omega1
          END IF
        END IF
      END IF
      IF (diftim < eps) THEN
        cosz(j,k) = 0.0
        lit(j,k) = 0.0
      ELSE
        cosz(j,k) = difsin*coscos(j)/diftim + sinsin(j)
        lit(j,k) = diftim / dtrad
      END IF
    END DO
  END DO

  ! At TOA only consider radiation from below the horizontal
  k=n_layer+1
  DO j = 1, n_profile
    cos_omega_horizon = &
      (-SQRT(1.0-planet_radius**2/r_layer(j,k)**2)-sinsin(j)) / coscos(j)
    cos_omega_90 = -sinsin(j) / coscos(j)
    
    IF ( cos_omega_horizon > 1.0 ) THEN
      ! Perpetual night
      diftim = 0.0
    ELSE IF ( cos_omega_90 < -1.0 ) THEN
      ! Sun will always be above the horizontal here (perpetual day)
      diftim = 0.0
    ELSE IF ( cos_omega_horizon < -1.0 .AND. cos_omega_90 > 1.0 ) THEN
      ! Perpetual twilight
      omega1 = hat(j)
      omega2 = hat(j) + dtrad
      difsin = SIN(omega2) - SIN(omega1)
      diftim = dtrad
    ELSE
      ! At this latitude some points are lit from below, some not.
      ! There are two periods when the sun is below the horizontal:
      ! just after sunrise and just before sunset. We first calculate the
      ! hour angle boundaries of these two periods along with the hour angles
      ! for the beginning and end of the timestep. The omega1,2,3,4 hour
      ! angles can then be determined to give the boundaries of when these
      ! periods intersect with timestep.
      ! The angles are taken relative to sunrise in the range 0 to 2pi to
      ! allow tests on their orders to work.

      ! First calculate the hour angles where the sun is on the horizon
      IF ( cos_omega_horizon < -1.0 ) THEN
        ! Sun is never below the horizon (but sometimes above the horizontal)
        ! so take sunrise and sunset to both be at midnight.
        omega_sunset  = pi
      ELSE
        omega_sunset  = ACOS(cos_omega_horizon)
      END IF
      omega_sunrise = -omega_sunset
      omegas = omega_sunset - omega_sunrise

      ! Now calculate the hour angles where the sun is horizontal
      IF ( cos_omega_90 > 1.0 ) THEN
        ! Sun is always below the horizontal but sometimes in twilight
        ! so take these angles to both be at midday.
        omega_90_setting = 0.0
        omega_90_rising = 0.0
      ELSE
        omega_90_setting = MODULO(ACOS(cos_omega_90) - omega_sunrise, twopi)
        omega_90_rising = MODULO(-ACOS(cos_omega_90) - omega_sunrise, twopi)
      END IF

      ! Beginning and end of timestep relative to sunrise
      omegab = MODULO(hat(j) - omega_sunrise, twopi)
      omegae = MODULO(omegab + dtrad, twopi)

      IF (omegab < omega_90_rising) THEN
        omega1 = omegab + omega_sunrise
      ELSE IF (omegae < omegab) THEN
        omega1 = omega_sunrise
      ELSE
        omega1 = omega_90_rising + omega_sunrise
      END IF
      IF (omegae < omega_90_rising) THEN
        omega2 = omegae + omega_sunrise
      ELSE
        omega2 = omega_90_rising + omega_sunrise
      END IF
      IF (omegab >= omega_90_setting .AND. omegab < omegas) THEN
        omega3 = omegab + omega_sunrise
      ELSE IF ( (omegab < omega_90_setting .AND. omegae > omega_90_setting) &
           .OR. (omegab < omega_90_setting .AND. omegae < omegab) &
           .OR. (omegab > omegas .AND. omegae > omega_90_setting .AND. &
                 omegae < omegab) ) THEN
        omega3 = omega_90_setting + omega_sunrise
      ELSE
        omega3 = omegas + omega_sunrise
      END IF
      IF (omegae > omega_90_setting .AND. omegae < omegas) THEN
        omega4 = omegae + omega_sunrise
      ELSE
        omega4 = omegas + omega_sunrise
      END IF
      difsin = SIN(omega2) - SIN(omega1) + SIN(omega4) - SIN(omega3)
      diftim = omega2 - omega1 + omega4 - omega3
    END IF
    IF (diftim < eps) THEN
      cosz(j,k) = 0.0
      lit(j,k) = 0.0
    ELSE
      cosz(j,k) = difsin*coscos(j)/diftim + sinsin(j)
      lit(j,k) = diftim / dtrad
    END IF
  END DO

END IF

DO j = 1, n_profile
  IF (cosz(j,0) <= 0.0) THEN
    ! The surface can only receive radiation from above
    lit(j,0) = 0.0
  ELSE IF (l_orog) THEN
    IF (slope_angle(j) /= 0.0) THEN
      ! Adjust the lit fraction for the area of the sloping surface
      ! aligned with the solar beam
      lit(j,0) = lit(j,0) * MAX( 1.0 + TAN( ACOS(cosz(j,0)) ) * &
        TAN(slope_angle(j)) * COS( sol_azimuth(j) - slope_aspect(j) ), 0.0 )
    END IF
  END IF
  ! At TOA only consider radiation from below
  IF (cosz(j,n_layer+1) >= 0.0) lit(j,n_layer+1) = 0.0
END DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE solang_sph
END MODULE solang_sph_mod
