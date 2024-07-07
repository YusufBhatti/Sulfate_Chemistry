! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Invert relationship between aerosol, visibility and water content

! Description:
!   Invert relationship between aerosol, visibility and water
!   content. This is needed for fog probability calculation.

! Documentation:
!    Wright, B. J., 1997: Improvements to the Nimrod Visibility
!       Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!    Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!       for Nimrod. Met. Office FR Tech Rep., No. 222.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Atmosphere Service

! Code description:
!   This code is written to UMDP3 standards.

! Subroutine Interface:
SUBROUTINE vistoqt                                                      &
           (visibility,                                                 &
                                                      !INPUT
            qs,                                                         &
                                                      !INPUT
            aerosol,                                                    &
                                                      !INPUT
            l_murk,                                                     &
                                                      !INPUT
            npoints,                                                    &
                                                      !INPUT
            qt )

! Modules
USE visbty_constants_mod, ONLY: m0, aero0, n0, power, radius0, a0, b0,  &
                                visfactor, recipvisair, fourthirds
USE conversions_mod, ONLY: pi
USE water_constants_mod, ONLY: rho_water
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

!  -------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                  &
        npoints             ! IN NO. points in field.
REAL, INTENT(IN) ::                                                     &
        visibility(npoints),                                            &
                               ! IN visibility
        qs(npoints),                                                    &
                                !  Saturated humidity mixing ratio
        aerosol(npoints)    ! IN Aerosol mixing ratio(ug/kg)
LOGICAL, INTENT(IN) ::                                                  &
        l_murk                    ! IN : Aerosol present
!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
REAL, INTENT(OUT) ::                                                    &
        qt(npoints)         ! OUT Total water mixing ratio (kg/kg)
!---------------------------------------------------------------------
! Local parameter variables
REAL, PARAMETER :: qt_limit = 0.0001

! Local workspace variables
INTEGER ::                                                              &
 point            !  Loop variable for points

REAL ::                                                                 &
  ql,                                                                   &
                       !  Liquid water mixing ratio (Kg/Kg).
  radius_dry,                                                           &
                       !  Dry particle radius for aerosol (m)
  radius,                                                               &
                       !  Radius of fog droplets (m)
  radius_star,                                                          &
                       !  Water droplet radius divided by dry radius
  radius_act,                                                           &
                       !  Activation droplet radius
  radius_star_act,                                                      &
                       !  Activation droplet radius divided by the dry
  radius_star_used,                                                     &
                       !  Water droplet radius divided by the dry
                       !   radius actually used for the relative
                       !   humidity calculation
  rh,                                                                   &
                       !  Relative humidity derived from visibility
  a,                                                                    &
                       !  A0 divided by the dry radius
  m_over_m0,                                                            &
                       !  Ratio of the aerosol mass mixing ratio and
                       !   the standard aerosol mass mixing ratio
  n                !  Number density of aerosol particles (/m3)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='VISTOQT'

! *

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)
  !---------------------------------------------------------------------
  !  1.  Calculate the ratio of the aerosol mass mixing ratio to the
  !      standard mass mixing ratio, m_over_m0, and the aerosol number
  !      density, N:
  !
  !                       (1-3p)
  !                   (m )
  !            N = N0 (--)
  !                   (m0)
  !
  !      And the dry radius, radius_dry:
  !                       p
  !                   (m )
  !            r = r0 (--)
  !             d     (m0)
  !
  !      And A (A0 divided by the dry radius).
  !
  !      And the activation radius:
  !
  !                              1/2
  !                   (       3 )
  !                   ( 3 B0 r  )
  !            r    = ( ------d-)
  !             act   (   A0    )
  !
  ! N.B. AEROSOL is in ug/kg, m in kg/kg
  ! If not available, use 10 ug/kg
  !-------------------------------------------------------------------

! If l_murk=false, this calculation doesn't depend on points,
! so can be done outside of loop over points.
! Must be kept identical with code inside loop.
IF (.NOT. l_murk) THEN
  m_over_m0 = MAX(10.0/m0*1.0e-9, 0.0001)

  n = n0 * m_over_m0**(1.0-3*power)

  radius_dry = radius0 * (m_over_m0)**power
  a = a0 / radius_dry

  radius_act = SQRT( (3 * b0 * radius_dry**3) / a0 )
  radius_star_act =  radius_act/radius_dry
END IF

DO point = 1, npoints

  IF (l_murk) THEN
    ! Ensure that assumed aerosol conc. is at least Aero0:
    m_over_m0 = MAX(aerosol(point)/m0*1.0e-9,                           &
                             aero0/m0*1.0e-9, 0.0001)

    n = n0 * m_over_m0**(1.0-3*power)

    radius_dry = radius0 * (m_over_m0)**power
    a = a0 / radius_dry

    radius_act = SQRT( (3 * b0 * radius_dry**3) / a0 )
    radius_star_act =  radius_act/radius_dry
  END IF

  !----------------------------------------------------------------------
  !  2.  Calculate a water droplet radius, from the visibility:
  !
  !                                 1/2
  !                ( ln( epsilon ) )
  !            r = (---------------)
  !                (  Vis N Beta0  )
  !
  !     (An extra term RecipVisAir is included in the recipical of
  !      visibility to limit visibilities to 100km in clean air).
  !----------------------------------------------------------------------

  radius =(visfactor/n) *                                               &
     (1.0/visibility(point) - recipvisair)
  IF (radius  >   0.0) THEN
    radius = SQRT( radius )
  ELSE
    radius = radius_dry
  END IF

  !----------------------------------------------------------------------
  !  3.  Provided the diagnosed radius is greater than the dry radius,
  !      calculate the normalised droplet radius, and the saturated
  !      humidity mixing ratio.
  !----------------------------------------------------------------------

  IF ( radius  >   radius_dry ) THEN

    radius_star = radius / radius_dry

    !----------------------------------------------------------------------
    !  5.  Calculate the corresponding liquid water mixing ratio:
    !
    !
    !                4            (  3     3 )
    !           qL = - Pi rho_w N ( r  - r   )
    !                3            (       d  )
    !
    !----------------------------------------------------------------------

    ql = fourthirds * pi * rho_water * n  *                             &
         ( radius**3 - radius_dry**3 )

    !----------------------------------------------------------------------
    !  6.  Calculate the relative humidity:
    !
    !                   ( A        B0   )
    !           RH = exp( --  -  ------ )
    !                   ( r       3     )
    !                   (  *     r  - 1 )
    !                   (         *     )
    !
    !----------------------------------------------------------------------

    IF ( radius_star  <   radius_star_act ) THEN
      rh = EXP( a/radius_star                                           &
                - b0 /( radius_star **3 - 1.0 ) )
    ELSE
      rh = EXP( a/radius_star_act                                       &
                  - b0 /( radius_star_act **3 - 1.0 ) )
    END IF

    !----------------------------------------------------------------------
    !  7.  Calculate the total water mixing ratio:  qt = RH * qs(T) + qL
    !----------------------------------------------------------------------

    qt(point) = MAX( rh * qs(point) + ql, qt_limit )

    !----------------------------------------------------------------------
    !  8. If the droplet radius is less than the dry radius, then set the
    !     total water mixing ratio to the minimum value.
    !----------------------------------------------------------------------

  ELSE

    qt(point) = qt_limit

  END IF

END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE vistoqt
