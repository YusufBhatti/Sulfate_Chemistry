! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Calculate geopotential height for baroclinic wave test.

MODULE rcf_ideal_baro_geo_p_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_IDEAL_BARO_GEO_P_MOD'

CONTAINS

REAL FUNCTION baro_geo_p(xi2, eta)

USE parkind1,            ONLY: jpim, jprb       !DrHook
USE yomhook,             ONLY: lhook, dr_hook   !DrHook
USE planet_constants_mod, ONLY: r, planet_radius, g, omega
USE conversions_mod,     ONLY: pi
USE rcf_ideal_baroclinic_constants_mod, ONLY: &
    eta_t, eta0, u0, dt, t0, eg_lapse, l_channel, b, Ly, f0, beta0

IMPLICIT NONE
!
! Description:
!
! Calculates geopotential height at location given by
! xi2=latitude in radians and eta=normalised
! pressure coordinate.
!
!
! Method:
!
! Baroclinic instability test case (QJRMS 132, 2943--2975).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Idealised
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

! Arguments

REAL, INTENT(IN) :: xi2 ! Latitude (radians)
REAL, INTENT(IN) :: eta ! p/p_surface

! Local variables
REAL :: eta_v
REAL :: av_geo
REAL :: const
REAL :: wind_term
REAL :: rotation_term
REAL :: phi_prime

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BARO_GEO_P'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

eta_v = (eta - eta0) * pi / 2.0

! Zonally averaged geopoential height in tropopause
av_geo = (t0*g/eg_lapse)*(1.0-eta**(r*eg_lapse/g))

IF ( L_channel ) THEN
  const = 2.0*pi*xi2/Ly

  phi_prime = 0.5*u0*((f0-beta0*Ly/2.0)*(xi2-Ly/2.0-Ly/(2.0*pi)*SIN(const)) &
            + 0.5*beta0*(xi2**2-Ly*xi2/pi*SIN(const)                        &
            - Ly**2/(2.0*pi**2)*COS(const) - Ly**2/3.0 - Ly**2/(2.0*pi**2)))

  baro_geo_p = av_geo + phi_prime*LOG(eta)*EXP(-(LOG(eta)/b)**2)
ELSE
  ! Zonally averaged geopotential height in stratosphere
  IF (eta < eta_t) THEN
    av_geo = av_geo - r*dt*( (LOG(eta/eta_t)+137.0/60.0)*eta_t**5       &
                               - 5.0*eta_t**4 * eta                     &
                               + 5.0*eta_t**3 * eta**2                  &
                               - (10.0/3.0)*eta_t**2 * eta**3           &
                               + (5.0/4.0)*eta_t * eta**4               &
                               - (1.0/5.0) * eta**5 )
  END IF

  const =  u0 * COS(eta_v)**(3.0/2.0)

  wind_term = (-2.0*SIN(xi2)**6 * (COS(xi2)**2 + 1.0/3.0) + 10.0/63.0) *  &
              u0*COS(eta_v)**(3.0/2.0)

  rotation_term = ( (8.0/5.0)*COS(xi2)**3 * (SIN(xi2)**2 + 2.0/3.0)       &
                    - pi/4.0 ) * planet_radius * omega

  baro_geo_p = av_geo + const*(wind_term + rotation_term)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END FUNCTION baro_geo_p
END MODULE rcf_ideal_baro_geo_p_mod
