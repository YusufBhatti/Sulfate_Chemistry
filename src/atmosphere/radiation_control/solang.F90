! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculation of solar zenith/azimuth angles and sunlit fraction.
!
! Purpose :
!  Calculations of the earth's orbit described in the second page of
!  the "Calculation of incoming insolation" section of UMDP 23, i.e.
!  from the sin of the solar  declination, the position of each point
!  and the time limits it calculates how much sunlight, if any, it
!  receives.
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE solang_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOLANG_MOD'
CONTAINS

SUBROUTINE solang (sindec, t, dt, eqt, lat, longit, k,                         &
     lit, cosz, sol_azimuth, cosz_beg, cosz_end)

USE conversions_mod, ONLY: pi
USE planet_constants_mod, ONLY: s2r, planet_dha,                               &
  l_fix_solang, solar_zenith_angle, solar_azimuth_angle
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE


INTEGER, INTENT(IN) ::                                                         &
  k
!   Number of points

REAL, INTENT(IN) ::                                                            &
  sindec,                                                                      &
!   Sin(solar declination)
  t, dt,                                                                       &
!   Start time (GMT) & timestep
  lat(k), longit(k),                                                           &
!   latitude & longitude of each point
  eqt
!   The equation of time (the difference between true and mean solar time)

REAL, INTENT(OUT) ::                                                           &
  lit(k),                                                                      &
!   Sunlit fraction of the timestep
  cosz(k),                                                                     &
!   Mean cos(solar zenith angle) during the sunlit fraction
  sol_azimuth(k),                                                              &
!   Mean solar azimuth angle (radians clockwise from grid north) during the
!   sunlit fraction
  cosz_beg(k), cosz_end(k)
!   Cosine of the solar zenith angle at the beginning and end of the
!   period over which cosz is integrated


INTEGER ::                                                                     &
  j
!   Loop counter over points

REAL ::                                                                        &
  mean_omega(k),                                                               &
!   Mean hour angle over the timestep measured in radians west of local noon
  sinsin, coscos,                                                              &
!   Products of the sines and cosines of solar declination and latitude
  hld, coshld,                                                                 &
!   Half-length of the day in radians (equal to the hour-angle of sunset,
!   and minus the hour-angle of sunrise) & its cosine.
  hat,                                                                         &
!   Local hour angle at the start time.
  omegab, omegae, omega1, omega2, omegas,                                      &
!   Local hour angle at the beginning and end of the timestep and
!   of the period over which cosz is integrated, and sunset - all measured in
!   radians after local sunrise, not from local noon as the true hour angle is.
  difsin, diftim,                                                              &
!   A difference-of-sines intermediate value and the corresponding time period
  trad, dtrad,                                                                 &
!   These are the start-time and length of the timestep (T & DT)
!   converted to radians after midday GMT, or equivalently, hour
!   angle of the mean sun on the Greenwich meridian.
  sinlat(k), x(k), dec
!   Working variables for bearing calculation.

REAL, PARAMETER :: twopi = 2.0*pi
REAL, PARAMETER :: eps = EPSILON(1.0)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SOLANG'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ABS(planet_dha) < eps) THEN

  ! Planet is tidally locked.
  ! The longitude of the overhead sun may vary with the equation of time
  ! to allow for an eccentric orbit and a constant rotation rate.
  mean_omega = MODULO(longit + eqt - pi, twopi)

  ! The solar declination is allowed to vary with the orbital parameters.
  dec = ASIN(sindec)
  cosz = sindec * SIN(lat) + COS(dec)*COS(lat)*COS(mean_omega)
  WHERE (cosz < eps)
    cosz = 0.0
    lit = 0.0
  ELSEWHERE
    lit = 1.0
  END WHERE
  cosz_beg = cosz
  cosz_end = cosz
  sol_azimuth = MODULO(ATAN2( -COS(dec)*SIN(mean_omega),                       &
    COS(lat)*sindec - SIN(lat)*COS(dec)*COS(mean_omega) ), twopi)

ELSE IF (l_fix_solang) THEN

  ! Fix the sun at a particular zenith and azimuth angle for all points
  ! for idealised tests.
  cosz = COS(solar_zenith_angle)
  WHERE (cosz < eps)
    cosz = 0.0
    lit = 0.0
  ELSEWHERE
    lit = 1.0
  END WHERE
  cosz_beg = cosz
  cosz_end = cosz
  sol_azimuth = solar_azimuth_angle

ELSE

  trad = t * s2r - pi
  dtrad = dt * s2r
  !DIR$ IVDEP
  DO j = 1, k                          ! Loop over points
    coshld = 0.0
    hld = 0.0                                ! Logically unnecessary
    ! statement without which the CRAY compiler will not vectorize this code
    sinlat(j) = SIN(lat(j))
    sinsin = sindec * sinlat(j)
    coscos = SQRT( (1.0-sindec**2) * (1.0-sinlat(j)**2) )
    IF ( sinsin  <   -coscos ) THEN                   ! Perpetual nig
      lit(j) = 0.0
      cosz(j) = 0.0
      cosz_beg(j) = 0.0
      cosz_end(j) = 0.0
      mean_omega(j) = 0.0
    ELSE
      !         Since T and DT represent mean time and all solar calculations
      !         are done using the variable HAT, it suffices to add the
      !         equation of time on to HAT.
      hat = longit(j) + trad + eqt         ! (3.2.2)
      IF ( sinsin  >   coscos ) THEN
        omega1 = hat                      ! angles for (3.2.3) are
        omega2 = hat + dtrad              ! start & end of timestep
      ELSE                                !   At this latitude some
        coshld = sinsin / coscos
        ! points are sunlit, some not.  Different ones need different treatment.
        hld = ACOS(-coshld)               ! (3.2.4)
        ! The logic seems simplest if one takes all "times" - actually hour
        ! angles - relative to sunrise (or sunset), but they must be kept in the
        ! range 0 to 2pi for the tests on their orders to work.
        omegab = MODULO(hat + hld, twopi)
        omegae = MODULO(omegab + dtrad, twopi)
        omegas = 2.0 * hld
        ! Now that the start-time, end-time and sunset are set in terms of hour
        ! angle, can set the two hour-angles for (3.2.3).  The simple cases are
        ! start-to-end-of-timestep, start-to-sunset, sunrise-to-end and sunrise-
        ! -to-sunset, but two other cases exist and need special treatment.
        IF (omegab <= omegas .OR. omegab <  omegae) THEN
          omega1 = omegab - hld
        ELSE
          omega1 = - hld
        END IF
        IF (omegae <= omegas) THEN
          omega2 = omegae - hld
        ELSE
          omega2 = omegas - hld
        END IF
        IF (omegae >  omegab .AND. omegab >  omegas) omega2=omega1
        !  Put in an arbitrary marker for the case when the sun does not rise
        !  during the timestep (though it is up elsewhere at this latitude).
        !  (Cannot set COSZ & LIT within the ELSE ( COSHLD < 1 ) block
        !  because 3.2.3 is done outside this block.)
      END IF           ! This finishes the ELSE (perpetual day) block
      difsin = SIN(omega2) - SIN(omega1)             ! Begin (3.2.3)
      diftim = omega2 - omega1
      mean_omega(j) = (omega1 + omega2)/2.0
      ! Next, deal with the case where the sun sets and then rises again
      ! within the timestep.  There the integration has actually been done
      ! backwards over the night, and the resulting negative DIFSIN and DIFTIM
      ! must be combined with positive values representing the whole of the
      ! timestep to get the right answer, which combines contributions from
      ! the two separate daylit periods.  A simple analytic expression for the
      ! total sun throughout the day is used.  (This could of course be used
      ! alone at points where the sun rises and then sets within the timestep)
      IF (diftim <  0.0) THEN
        difsin = difsin + 2.0 * SQRT(1.0-coshld**2)
        diftim = diftim + 2.0 * hld
        mean_omega(j) = mean_omega(j) + pi
      END IF
      IF (mean_omega(j) >  pi) mean_omega(j)=mean_omega(j)-twopi
      IF (diftim == 0.0) THEN
        ! Pick up the arbitrary marker for night points at a partly-lit latitude
        cosz(j) = 0.0
        cosz_beg(j) = 0.0
        cosz_end(j) = 0.0
        lit(j) = 0.0
      ELSE
        cosz(j) = difsin*coscos/diftim + sinsin     ! (3.2.3)
        cosz_beg(j) = coscos*COS(omega1) + sinsin
        cosz_end(j) = coscos*COS(omega2) + sinsin
        lit(j) = diftim / dtrad
      END IF
    END IF            ! This finishes the ELSE (perpetual night) block
  END DO

  dec = ASIN(sindec)
  sol_azimuth = MODULO(ATAN2( -COS(dec)*SIN(mean_omega),                       &
      COS(lat)*sindec - sinlat*COS(dec)*COS(mean_omega) ), twopi)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE solang
END MODULE solang_mod
