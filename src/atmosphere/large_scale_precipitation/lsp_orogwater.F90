! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Large-scale precipitation scheme. Subgrid orographic water calculation

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Subroutine Interface:
MODULE lsp_orogwater_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_OROGWATER_MOD'

CONTAINS

SUBROUTINE lsp_orogwater(                                                     &
  points,                                                                     &
                     ! Number of points
  hmteff,                                                                     &
                     ! Effective peak-to-trough mountain height (m)
  zb,                                                                         &
                     ! Blocked layer depth (m)
  r_theta_levels_c,r_theta_surf_c, fv_cos_theta_latitude_c,                   &
                     ! For obtaining altitude
  p, t, q, qsl, esw, qcl, cfliq,                                              &
                     ! Grid-box thermodynamics
  ql_orog, cf_orog                                                            &
                     ! Orographic cloud water MR (kg/kg) and fraction
  )

USE um_types,             ONLY: real_lsprec
USE lsprec_mod,           ONLY: zerodegc, pi, delta_lambda, delta_phi, lc,    &
                                qcfmin, repsilon, r, cp, g, rv,               &
                                zero, one, two, half,                         &
! Horizontal scale of sub-grid sinusoidal ridge is (nscalesf * grid-length)
! in the vertical decay equation - longer scales decay more slowly
                                 nscalesf

! Dr Hook Modules
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------------
! Purpose:
!   Calculates sub-grid orographic water MR and cloud fraction
!   No standard model variables are modified here.
!   Assumes moist neutral flow over a 2D ridge of height determined by sd_orog
!   Its horizontal wavelength is assumed equal to the grid-spacing, which
!   determines the rate at which orographic displacements decay with altitude.
!   If air is already saturated and some cloud water is present then the LCL
!   is assumed to lie below the model level at the distance required to
!   evaporate qcl,  Otherwise the LCL is above the model level by an amount
!   determined by the difference in the T and Td lapse rates.

!   Documentation: UMDP 026

!------------------------------------------------------------------------------

! Subroutine Arguments

INTEGER, INTENT(IN) :: points
!   Number of points calculation is done over

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  hmteff(points),                                                             &
!         Effective mountain peak-to-trough height (metres)
  zb(points),                                                                 &
!         Sub-grid orographic blocked layer depth (metres)
  p(points),                                                                  &
!         Air pressure at this level (Pa).
  qcl(points),                                                                &
!         Cloud liquid water (kg water per kg air).
  cfliq(points),                                                              &
!         Liquid cloud fraction (0 to 1).
  t(points),                                                                  &
!         Temperature at this level (K).
  q(points),                                                                  &
!         Specific humidity at this level (kg water per kg air).
  qsl(points),                                                                &
!         Saturation specific humidity (kg water per kg air).
  esw(points),                                                                &
!         Saturation vapour pressure wrt water at all T (Pa).
  r_theta_levels_c(points),                                                   &
                     ! Distance from centre of the Earth
  r_theta_surf_c(points),                                                     &
                     ! ...and near surface (k=0)
  fv_cos_theta_latitude_c(points)
                     ! Finite volume cosine of latitude.

REAL (KIND=real_lsprec), INTENT(OUT) ::                                       &
  ql_orog(points),                                                            &
!         Grid-box mean orographic cloud water MR (kg/kg)
  cf_orog(points)
!        Orographic cloud fraction

!------------------------------------------------------------------------------
! Local Variables

INTEGER :: i

REAL (KIND=real_lsprec) :: altitude(points)

REAL (KIND=real_lsprec) :: gridspacing,waveno,tdew,dzcond,salr,               &
        xlim,ydis,dqldz,maxdisp,gridspacing_surf,                             &
        tmidcld,pmidcld,esatmidcld,qsatmidcld,                                &
        dewlr,satdz,evapdist,satdzmin,hmin,tc,tdewc,tmidcldc

REAL (KIND=real_lsprec),PARAMETER ::                                          &
! Constant in Tv equation used within tmidcld equation
    d1 = 0.61,                                                                &
! Dry Adiabatic Lapse Rate (K/km)
    dalr = 9.8e-3,                                                            &
! Constants in dewpoint and saturation vapour pressure equations as
! suggested by Alduchov and Eskridge (1996)
    a1 = 17.625,                                                              &
    b1 = 243.04,                                                              &
! Multiple in saturation vapour pressure estimate (Pa)
    c1 = 610.94

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_OROGWATER'

!==============================================================================
! Start of calculations

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise all temporary variables used
  gridspacing = zero
  gridspacing_surf  = zero
  waveno = zero
  maxdisp = zero
  dqldz = zero
  evapdist = zero
  xlim = zero
  ydis = zero
  tdew = zero
  tdewc = zero
  tc = zero
  dzcond = zero
  salr = zero
  tmidcld = zero
  tmidcldc = zero
  pmidcld = zero
  esatmidcld = zero
  qsatmidcld = zero
  dewlr = zero
  satdz = zero
  hmin = zero
  satdzmin = zero

! Initialise arrays passed back to lsp_ice
  DO i = 1, points
    ql_orog(i) = zero
    cf_orog(i) = zero
  END DO  ! For each cloudy datapoint

! Minimum subgrid orography able to enhance rain
  hmin = one
  satdzmin = 0.1_real_lsprec

! Calculate sub-grid orographic water MR when needed
! (Need rain or snow to be able to enhance accretion)
  DO i = 1, points

!   Only create extra water if above blocked layer depth zb
    altitude(i) = MAX( zero, r_theta_levels_c(i)-r_theta_surf_c(i) )

!   Only give non-zero ql_orog if above the blocking level
!   Do NOT proceed if negative Q has been passed in as this will
!    produce NaNs and cause the model to crash
    IF ( altitude(i) > zb(i) .AND. (half*hmteff(i)) > hmin  .AND.             &
                     q(i) > zero ) THEN

!         Horizontal gridspacing at this altitude for averaging
          gridspacing = SQRT ( r_theta_levels_c(i) * delta_lambda             &
                             * r_theta_levels_c(i) * delta_phi                &
                             * fv_cos_theta_latitude_c(i)  )

!         Horizontal surface gridspacing to characterise orography
          gridspacing_surf = SQRT ( r_theta_surf_c(i) * delta_lambda          &
                                  * r_theta_surf_c(i) * delta_phi             &
                                  * fv_cos_theta_latitude_c(i)  )
          waveno = two * pi / gridspacing_surf

!         Max vertical orographic displacement accounting for exp decay
!         Effective hill height is peak-to-trough value so halve it
!         Horizontal wavelength used by decay function is scaled by
!         nscalesf (longer scales decay more slowly aloft)

          maxdisp = half*hmteff(i) * EXP(-altitude(i)*waveno/nscalesf)


!         Cloudy grid-box:  saturated with water present
          IF ( qcl(i) > qcfmin .AND. q(i)/qsl(i) > one )  THEN

!           Rate of change of Qsat per metre of saturated descent
            dqldz = MAX( zero,  g*(one + lc*q(i)/r/t(i))/                     &
               (cp + repsilon*lc*lc*q(i)/r/t(i)/t(i))                         &
               *(repsilon + qsl(i))*qsl(i)*lc/r/t(i)/t(i)                     &
               -qsl(i)*p(i)*g/(p(i) - esw(i))/r/t(i)  )

            evapdist = qcl(i) / dqldz

            IF (maxdisp > evapdist) THEN

              xlim = (one/waveno)*ACOS(-evapdist/maxdisp)
              ydis = maxdisp*two*SIN(waveno*xlim)/waveno                      &
                       - two*evapdist*((half*gridspacing)-xlim)
              ql_orog(i) = MAX(zero, dqldz*ydis/gridspacing)
!             Make sure orog water is no larger than available vapour
              ql_orog(i) = MIN( ql_orog(i), q(i) )

!             Cloud fraction reduced by descent, assuming cfliq=1
              cf_orog(i) = MAX(zero, two*xlim/gridspacing)

            ELSE
              ql_orog(i) = zero
              cf_orog(i) = zero
            END IF

          ELSE     ! Clear gridbox - no resolved cloud

!           Convert temperature to C for use in tdew
            tc = t(i)-zerodegc

!           Dewpoint temperature (C)
            tdewc = b1*( LOG(q(i)/qsl(i)) + (a1*tc/(b1 + tc)) ) /             &
                ( a1 - LOG(q(i)/qsl(i)) - (a1*tc/(b1 + tc)) )

!           Convert dewpoint temperature to K
            tdew = tdewc + zerodegc

!           Dewpoint temperature adiabatic lapse rate
            dewlr = (tdew * tdew * g * rv)/(lc * r * t(i))

!           Condensation level (ascent required to reach saturation)
            dzcond = MAX( (t(i) - tdew)/(dalr - dewlr), zero)

!           Only carry on if there is saturated ascent
!           Will need distance to mid-cloud level satdz
            satdz = half * (maxdisp - dzcond)

            IF (satdz > satdzmin) THEN

!             Estimate T and SALR at cloud base - use q which remains
!             constant during dry ascent and = (qsat at LCL)
              tmidcld = t(i) - (dzcond * dalr)

              salr = g * (one + (lc*q(i))/(r*tmidcld) ) /                     &
                     (cp + (lc*lc*q(i)*repsilon)/(r*tmidcld*tmidcld))

!             Find variables at mid-cloud level - further decrease T
              tmidcld = tmidcld - (satdz*salr)

!             New pressure from total dry+moist ascent (use mean Tv)
              pmidcld = p(i) * EXP( -( (dzcond + satdz) * g)/                 &
                          (r*half*(t(i) + tmidcld)*(one + d1*q(i)) ) )

!             Convert mid-cloud temperature to C
              tmidcldc = tmidcld - zerodegc

!             New esat and qat to get dqldz
              esatmidcld = c1 * EXP( a1*tmidcldc/(tmidcldc + b1) )

              qsatmidcld = repsilon * esatmidcld / pmidcld

              salr = g * (one + (lc*qsatmidcld)/(r*tmidcld) ) /               &
                     (cp + (lc*lc*qsatmidcld*repsilon)/                       &
                     (r*tmidcld*tmidcld))

              dqldz = MAX( zero, salr*(repsilon + qsatmidcld)                 &
                        *qsatmidcld*lc/r/tmidcld/tmidcld                      &
                        -qsatmidcld*pmidcld*g/                                &
                        (pmidcld - esatmidcld)/r/tmidcld )

              xlim = (one/waveno)*ACOS(dzcond/maxdisp)
              ydis = maxdisp*two*SIN(waveno*xlim)/waveno                      &
                        - (two*xlim*dzcond)

              ql_orog(i) = MAX(zero, dqldz*ydis/gridspacing)
              ql_orog(i) = MIN( ql_orog(i), q(i) )

              cf_orog(i) = MAX(zero, two*xlim/gridspacing)

            ELSE
              ql_orog(i) = zero
              cf_orog(i) = zero
            END IF ! Water streamline above corrected dzcond

          END IF  ! Clear or cloudy gridbox

    ELSE   !if sub-grid orog NOT present

        ql_orog(i) = zero
        cf_orog(i) = zero

    END IF

  END DO  ! For each cloudy datapoint

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_orogwater
END MODULE lsp_orogwater_mod
