! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Reflectivity calculation

MODULE mphys_reflec_mod
! Description:
! Calculates radar reflectivity in dBZ for all available
! hydrometer species.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MPHYS_REFLEC_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE mphys_reflec( points, rho, t, qgraup, qcf, qcf2, qrain, qcl, ndrop,&
                         cfice, cfliq, cfrain, tcg, tcgc,                     &
                         dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_l, dbz_r )

!Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod,           ONLY: cx, constp, kliq, kice, mm6m3, ref_lim,       &
                                ref_lim_lin, mr_lim, cf_lim, nd_lim, ref_mom, &
                                zerodegc,                                     &
                                zero, one, ten

!- logicals and integers
USE mphys_inputs_mod,     ONLY: l_psd, l_mcr_qgraup, l_no_cf
USE mphys_bypass_mod,     ONLY: l_crystals

USE lsp_moments_mod,      ONLY: lsp_moments

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,             ONLY: real_lsprec

! Dr Hook Modules
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

IMPLICIT NONE

!------------------------------------------------------------------------------
! Purpose:
!   Calculates any radar reflectivity required for diagnostics
!   Documentation: UMDP 26.

!------------------------------------------------------------------------------
! Subroutine Arguments

INTEGER, INTENT(IN) :: points      ! Number of points calculation is done over

REAL (KIND=real_lsprec), INTENT(IN) :: rho(points)
  ! Air density [kg m-3]
REAL (KIND=real_lsprec), INTENT(IN) :: t(points)
  ! Temperature [K]
REAL (KIND=real_lsprec), INTENT(IN) :: qgraup(points)
  ! Graupel mixing ratio [kg kg-1]
REAL (KIND=real_lsprec), INTENT(IN) :: qcf(points)
  ! Ice Agg mixing ratio [kg kg-1]
REAL (KIND=real_lsprec), INTENT(IN) :: qcf2(points)
  ! Ice Cry mixing ratio [kg kg-1]
REAL (KIND=real_lsprec), INTENT(IN) :: qrain(points)
  ! Rain mixing ratio [kg kg-1]
REAL (KIND=real_lsprec), INTENT(IN) :: qcl(points)
  ! Cloud liquid mixing ratio [kg kg-1]
REAL (KIND=real_lsprec), INTENT(IN) :: ndrop(points)
  ! Cloud droplet number [m-3]
REAL (KIND=real_lsprec), INTENT(IN) :: cfice(points)
  ! Ice cloud fraction
REAL (KIND=real_lsprec), INTENT(IN) :: cfrain(points)
  ! Rain fraction
REAL (KIND=real_lsprec), INTENT(IN) :: cfliq(points)
! Liquid cloud fraction

REAL (KIND=real_lsprec), INTENT(IN) :: tcg(points)
! Agg temperature intercept function
REAL (KIND=real_lsprec), INTENT(IN) :: tcgc(points)
! Cry temperature intercept function

REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_tot(points)
  ! Total reflectivity [dBZ]
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_g(points)
  ! Reflectivity due to graupel [dBZ]
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_i(points)
  ! Reflectivity due to ice aggregates [dBZ]
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_r(points)
  ! Reflectivity due to rain [dBZ]
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_i2(points)
  ! Reflectivity due to ice crystals [dBZ]
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_l(points)
  ! Reflectivity due to liquid cloud [dBZ]

!------------------------------------------------------------------------------
! Local Variables

INTEGER :: i

REAL (KIND=real_lsprec) :: one_over_cfice(points)
  ! 1.0 / ice cloud fraction

REAL (KIND=real_lsprec) :: kgraup
  ! Reflectivity prefactor due to graupel
REAL (KIND=real_lsprec) :: krain
  ! Reflectivity prefactor due to rain
REAL (KIND=real_lsprec) :: kclw
  ! Reflectivity prefactor due to cloud liquid water
REAL (KIND=real_lsprec) :: kice_a
  ! Reflectivity prefactor due to ice aggregates
REAL (KIND=real_lsprec) :: kice_c
  ! Reflectivity prefactor due to ice crystals
REAL (KIND=real_lsprec) :: gwc(points)
  ! Graupel water content [kg m-3]
REAL (KIND=real_lsprec) :: iwc(points)
  ! Ice water content: Aggregates [kg m-3]
REAL (KIND=real_lsprec) :: iwc2(points)
  ! Ice water content: Aggregates [kg m-3]
REAL (KIND=real_lsprec) :: rwc(points)
  ! Rain water content [kg m-3]
REAL (KIND=real_lsprec) :: lwc(points)
  ! Liquid water content [kg m-3]

REAL (KIND=real_lsprec) :: lamr(points)
  ! Lambda (slope parameter) for rain [m-1]
REAL (KIND=real_lsprec) :: lamg(points)
  ! Lambda (slope parameter) for graupel [m-1]
REAL (KIND=real_lsprec) :: lamic(points)
  ! Lambda (slope parameter) for ice cry [m-1]
REAL (KIND=real_lsprec) :: lamia(points)
  ! Lambda (slope parameter) for ice agg [m-1]
REAL (KIND=real_lsprec) :: ze_g(points)
  ! Linear reflectivity due to graupel [mm6 m-3]
REAL (KIND=real_lsprec) :: ze_i(points)
  ! Linear reflectivity due to ice agg [mm6 m-3]
REAL (KIND=real_lsprec) :: ze_i2(points)
  ! Linear reflectivity due to ice cry [mm6 m-3]
REAL (KIND=real_lsprec) :: ze_r(points)
  ! Linear reflectivity due to rain    [mm6 m-3]
REAL (KIND=real_lsprec) :: ze_l(points)
  ! Linear reflectivity due to liq cld [mm6 m-3]
REAL (KIND=real_lsprec) :: ze_tot(points)
  ! Total linear reflectivity [mm6 m-3]

REAL (KIND=real_lsprec) :: mom4(points)
  ! 4th moment of the ice particle size distribution

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='MPHYS_REFLEC'

!==============================================================================
! Start of calculations
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise all output dBZ variables. These are initialised to the
! minimum reflectivity as this is usually less than zero.
dbz_tot(:) = ref_lim
dbz_g(:)   = ref_lim
dbz_i(:)   = ref_lim
dbz_r(:)   = ref_lim
dbz_i2(:)  = ref_lim
dbz_l(:)   = ref_lim

! Initialise local variables:
ze_g(:)   = zero
ze_i(:)   = zero
ze_i2(:)  = zero
ze_r(:)   = zero
ze_l(:)   = zero
ze_tot(:) = zero

one_over_cfice(:) = zero

kgraup = kice / 0.93_real_lsprec
krain  = kliq / 0.93_real_lsprec
kclw   = kliq / 0.93_real_lsprec
kice_a = kice / 0.93_real_lsprec
kice_c = kice / 0.93_real_lsprec

! Generate one_over_cfice - used in a few places.

IF (l_no_cf) THEN

  one_over_cfice(:) = one

ELSE ! l_no_cf

  DO i = 1, points
    IF (cfice(i) > cf_lim) one_over_cfice(i) = one / cfice(i)
  END DO

END IF ! l_no_cf

IF (l_psd) THEN

  !Convert to correct precision

  CALL lsp_moments(points, rho, t, qcf, one_over_cfice,                       &
                          ref_mom, mom4)

  DO i = 1, points

    IF ( qcf(i) > mr_lim ) THEN

      ! This will already have an in-cloud value due to one_over_cfice
      ! being passed into lsp_moments above. However, need to ensure
      ! cloud fraction is dealt with appropriately below.

      IF (l_no_cf) THEN
        ze_i(i) = mm6m3 * cx(90) * mom4(i)
      ELSE
        ze_i(i) = cfice(i) * mm6m3 * cx(90) * mom4(i)
      END IF
      ! Here cx(90) = 0.224 * (6.0 * ai/pi/900)**2

    END IF ! qcf(i) > mr_lim

  END DO ! points

ELSE ! not l_psd

  ! Two ice categories (qcf and qcf2) need to be treated separately

  DO i = 1, points

    ! Non-psd Crystals.
    IF (l_crystals .AND. qcf2(i) > mr_lim .AND.                               &
        one_over_cfice(i) /= zero ) THEN

      ! Work out ice crystal water content
      IF (l_no_cf) THEN
        iwc2(i)  = qcf2(i) * rho(i)
      ELSE
        iwc2(i)  = qcf2(i) * rho(i) * one_over_cfice(i)
      END IF ! l_no_cf

      ! Determine the slope parameter, Lambda
      ! Note that constp(101) does not include the temperature
      ! intercept from lsp_init. Therefore this must be added
      ! at this stage - see UMDP26 section on inferring PSD relations
      ! from mixing ratios and fluxes
      ! N.B. cx(94) =  1.0 / (bic + 1.0 + x4ic - x2ic)

      lamic(i) = ( ( tcgc(i) * constp(101) ) / iwc2(i) ) ** cx(94)

      ! Calculate linear radar reflectivity
      ! N.B. cx(97)      = -1.0 * (1.0 + x4ic + (2.0*bic) - x2ic)
      !      constp(104) = x1ic * (aic**2) * GAMMA(1.0 + x4ic + ( 2.0 * bic))
      !                         * ( 6.0 /(pi * rho_i2))**2
      ! rho_i2 = 900.0 kg m-3

      IF (l_no_cf) THEN
        ze_i(i) = mm6m3 * kice_c * constp(104) * lamic(i)**cx(97)
      ELSE
        ze_i(i) = mm6m3 * kice_c * cfice(i) * constp(104) * lamic(i)**cx(97)
      END IF

    END IF ! l_crystals

    ! Non-psd Aggregates.
    IF (qcf(i) > mr_lim .AND. one_over_cfice(i) /= zero ) THEN

      ! Work out ice crystal water content
      IF (l_no_cf) THEN
        iwc(i) = qcf(i) * rho(i)
      ELSE
        iwc(i) = qcf(i) * rho(i) * one_over_cfice(i)
      END IF ! l_no_cf

      ! Determine the slope parameter, Lambda
      ! Note that constp(102) does not include the temperature
      ! intercept from lsp_init. Therefore this must be added
      ! at this stage - see UMDP26 section on inferring PSD relations
      ! from mixing ratios and fluxes
      ! N.B. cx(95) =  1.0 / (bi + 1.0 + x4i - x2i)

      lamia(i) = ( (tcg(i) * constp(102)) / iwc(i) ) ** cx(95)

      ! Calculate linear radar reflectivity
      ! N.B. cx(98)      = -1.0 * (1.0 + x4i + (2.0*bi) - x2i)
      !      constp(105) = x1i  * (ai**2) * GAMMA(1.0 + x4i + ( 2.0 * bi))
      !                         * ( 6.0 /(pi * rho_i))**2
      ! rho_i = 900.0 kg m-3

      IF (l_no_cf) THEN
        ze_i(i)  = mm6m3 * kice_a * constp(105) * lamia(i)**cx(98)
      ELSE
        ze_i(i)  = mm6m3 * kice_a * cfice(i) * constp(105) * lamia(i)**cx(98)
      END IF ! l_no_cf

    END IF ! qcf(i) > mr_lim

  END DO ! points

END IF ! l_psd

DO i = 1, points

  ! Graupel

  ! In the UM graupel does not have a cloud fraction (rightly or wrongly)
  ! i.e. the assumption is that cloud fraction is always 1.0

  ! Thus there is no need to calculate an in-cloud amount. Where the ice
  ! cloud fraction is zero (e.g. below the freezing level, the graupel
  ! reflectivity will still occur).

  IF (l_mcr_qgraup .AND. qgraup(i) > mr_lim ) THEN

    ! Work out graupel water content (no cloud fraction assumed)
    gwc(i)  = qgraup(i) * rho(i)

    ! Determine the slope parameter, Lambda
    ! (  pi/6 rhog * x1r * GAMMA(4 + x4g) ) ** (1.0 / (bg  + 1.0 + x4g  - x2g))
    ! (---------------------------------- )
    ! (               GWC                 )

    lamg(i) = (constp(100) / gwc(i) ) ** cx(93)
    ! cx(93) = (1.0 / (bg  + 1.0 + x4g  - x2g))

    ! Calculate linear radar reflectivity (no cloud fraction assumed).
    ze_g(i) = mm6m3 * kgraup * constp(103) * lamg(i)**cx(96)
    ! constp(103) = x1g  * (ag**2)  * gref1x4g  * ( 6.0 /(pi * rho_g))**2
    ! cx(96)      = -1.0 * (1.0 + x4g  + (2.0*bg)  - x2g)

  END IF ! l_mcr_qgraup

  ! Rain
  IF (qrain(i) > mr_lim .AND. cfrain(i) > cf_lim ) THEN

    ! Work out rain water content
    IF (l_no_cf) THEN
      rwc(i) = qrain(i) * rho(i)
    ELSE
      rwc(i) = qrain(i) * rho(i)/cfrain(i)
    END IF

    ! define a lamr, equivalent to the following equation

    ! (  pi/6 rhow * x1r * GAMMA(4 + x4r) ) ** (1.0 / (4 + x4r -x2r))
    ! ( ----------------------------------)
    ! (              RWC                  )

    lamr(i) = (constp(106) / rwc(i) ) ** cx(52)
    ! cx(52) = (1.0 / (4 + x4r -x2r))

    ! Then define the reflectivity in mm6 m-3.
    ! Reminder:
    ! constp(107) = gamma(1.0 + (2.0 * 3.0) +x4r) * x1r

    IF (l_no_cf) THEN
      ze_r(i) = mm6m3 * krain * constp(107) * lamr(i)**cx(92)
    ELSE
      ze_r(i) = cfrain(i) * mm6m3 * krain * constp(107) * lamr(i)**cx(92)
    END IF

  END IF !qrain > mr_lim

  ! Liquid Cloud
  IF (qcl(i) > mr_lim .AND. ndrop(i) > nd_lim .AND. cfliq(i) > cf_lim) THEN

    ! Work out liquid water content
    IF (l_no_cf) THEN
      lwc(i) = qcl(i) * rho(i)
    ELSE
      lwc(i) = qcl(i) * rho(i)/cfliq(i)
    END IF

    ! Equation A12 of Stein et al (2014), Monthly Weather Review.
    IF (l_no_cf) THEN
      ze_l(i) = mm6m3 * kclw * (201.6_real_lsprec/ndrop(i))                   &
                * (lwc(i)*cx(91))**2
    ELSE
      ze_l(i) = cfliq(i) * mm6m3 * kclw * (201.6_real_lsprec/ndrop(i))        &
                * (lwc(i)*cx(91))**2
    END IF

  END IF ! qcl > mr_lim

  ! Compute total reflectivity
  ze_tot(i) = ze_g(i) + ze_i(i) + ze_i2(i) + ze_r(i) + ze_l(i)

  ! Convert from linear (mm^6 m^-3) to dBZ
  IF (ze_tot(i) > ref_lim_lin) THEN
    dbz_tot(i) = MAX(ref_lim, ten * LOG10(ze_tot(i)))
  END IF

  IF (ze_g(i) > ref_lim_lin) THEN
    dbz_g(i) = MAX(ref_lim, ten * LOG10(ze_g(i)))
  END IF

  IF (ze_i(i) > ref_lim_lin) THEN
    dbz_i(i) = MAX(ref_lim, ten * LOG10(ze_i(i)))
  END IF

  IF (ze_i2(i) > ref_lim_lin) THEN
    dbz_i2(i) = MAX(ref_lim, ten * LOG10(ze_i2(i)))
  END IF

  IF (ze_r(i) > ref_lim_lin) THEN
    dbz_r(i) = MAX(ref_lim, ten * LOG10(ze_r(i)))
  END IF

  IF (ze_l(i) > ref_lim_lin) THEN
    dbz_l(i)   = MAX(ref_lim, ten * LOG10(ze_l(i)))
  END IF

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mphys_reflec
END MODULE mphys_reflec_mod
