! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Initialisation of variables
! Subroutine Interface:
MODULE lsp_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_INIT_MOD'

CONTAINS

SUBROUTINE lsp_init(                                                          &
  points,                                                                     &
                                     ! Number of points
  timestep, timestep_mp,                                                      &
                                     ! Timesteps
  t, p, cttemp,                                                               &
                                     ! Temperature and pressure
  deltaz,rhodz,rhodz_dry,rhodz_moist,                                         &
                                           ! Air density information
  rho, rhor,                                                                  &
  dhi, dhir,                                                                  &
                                     ! Tstep and layer thick. ratios
  q, qcl, qcf, qcf2, qrain, qgraup, rainrate,                                 &
                                                   ! Water contents
  qcf_agg, qcf_cry, qcf_tot, frac_agg, frac_cry_dep, frac_agg_dep,            &
  qs, qsl, esi, esw,                                                          &
                                     ! Saturation water quantities
  psdep, psaut, psacw, psacr,                                                 &
                                     ! Aggregate transfer rates
  psaci, psmlt, psmltevp,psfall,                                              &
  pifrw, pifrr, piprm, pidep, piacw,                                          &
                                     ! Ice transfer rates
  piacr, pimlt, pimltevp,pifall,                                              &
  praut, pracw, prevp, prfall,                                                &
                                     ! Rain transfer rates
  plset, plevpset,                                                            &
                                     ! Droplet settling transfers
  pgaut, pgacw, pgacs, pgmlt,                                                 &
                                     ! Graupel transfer rates
  pgfall,                                                                     &
  cf_transfer_diag, cfl_transfer_diag,                                        &
  cff_transfer_diag, rf_transfer_diag,                                        &
                                     ! Cloud and rain fraction
                                     ! transfer diagnostics
  snow_cry, snow_agg, snowt_cry,                                              &
                                     ! Precipitation rates
  snowt_agg, rainratet, graupratet,                                           &
  lheat_correc_liq, lheat_correc_ice,                                         &
                                           ! Latent heat corrections
  corr, corr2, rocor,                                                         &
                                     ! Fall speed and diffusivity
                                     ! corrections
  tcg,tcgi,tcgc,tcgci,tcgg,tcggi,                                             &
                                           ! Ice PSD intercepts
  cficekeep, vm_cry, vm_agg, l_use_agg_vt, vm_used,                           &
       ! Variables for different fallspeeds with generic psd
  niters_mp )

!Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod,         ONLY: t_scaling, m0, qcf0, t_agg_min, cx, constp,     &
                              tcor1, tcor2, zerodegc, lc, lf, cp, r,          &
                              repsilon, recip_epsilon, cpwr
USE lsprec_mod,         ONLY: zero, half, one

  ! Microphysics modules- logicals and integers
USE mphys_constants_mod, ONLY: l_cry_agg_dep
USE mphys_inputs_mod,    ONLY: l_psd, l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup,  &
                               l_diff_icevt
USE mphys_bypass_mod,    ONLY: l_crystals

! General and Atmospheric Modules- logicals and integers
USE gen_phys_inputs_mod,  ONLY: l_mr_physics

! Use in KIND for large scale precip, used for compressed variables passed down
! from here
USE um_types,             ONLY: real_lsprec

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    qsat_wat_new     => qsat_wat,                       &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_lsprec !Currently defaults to FALSE

! Dr Hook Modules
USE yomhook,             ONLY: lhook, dr_hook
USE parkind1,            ONLY: jprb, jpim

! Mathematical modules
USE vectlib_mod,   ONLY: powr_v => powr_v_interface,                          &
                         exp_v  => exp_v_interface

! Large scale precip modules
USE lsp_moments_mod, ONLY: lsp_moments

! Module for bug fixes
USE science_fixes_mod,    ONLY: l_fix_mphys_diags_iter

IMPLICIT NONE

! Purpose:
!   Perform initialization of variables required for the
!   microphysics scheme.

! Method:
!   Set variables to their initial values.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation


! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.



! Subroutine Arguments


INTEGER, INTENT(IN) ::                                                        &
  points,                                                                     &
                     ! Number of points
  niters_mp
                     ! Iterations of microphysics

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  timestep,                                                                   &
                     ! Model physics timestep / s
  timestep_mp,                                                                &
                     ! Timestep of each microphysics iteration / s
  t(points),                                                                  &
                     ! Temperature / K
  p(points),                                                                  &
                     ! Pressure / N m-2
  q(points),                                                                  &
                     ! Vapour content / kg kg-1
  qcl(points),                                                                &
                     ! Liquid water content / kg kg-1
  qcf(points),                                                                &
                     ! Ice aggregate content / kg kg-1
  qcf2(points),                                                               &
                     ! Ice crystal content / kg kg-1
  cttemp(points),                                                             &
                     ! Cloud top temperature / K
  rainrate(points),                                                           &
                       ! Rainrate into layer / kg m-2 s-1
  deltaz(points),                                                             &
                     ! Layer thickness / m
  rhodz(points),                                                              &
                     ! Air density * deltaz based on hydrostatic
                     ! assumption / kg m-2
  rhodz_dry(points),                                                          &
                         ! Air density*deltaz for dry air based on
                         ! non-hydrostatic assumption / kg m-2
  rhodz_moist(points),                                                        &
                         ! Air density*deltaz for moist air based on
  cficekeep(points)
                     ! input ice cloud fraction

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  qgraup(points),                                                             &
                     ! Graupel content / kg kg-1
  snow_cry(points),                                                           &
                        ! Ice crystal precip into layer / kg m-2 s-1
  snow_agg(points)  ! Ice agg. precip into layer / kg m-2 s-1

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
      ! Microphysical process rate diagnostics / kg kg-1 s-1
  psdep(points),                                                              &
                     ! Deposition of vapour to snow aggregates
  psaut(points),                                                              &
                     ! Autoconversion of aggregates from crystals
  psacw(points),                                                              &
                     ! Accretion of liq. water by snow aggregates
  psacr(points),                                                              &
                     ! Collection of rain by snow aggregates
  psaci(points),                                                              &
                     ! Collection of ice crystals by aggregates
  psmlt(points),                                                              &
                     ! Melting of snow aggregates
  psmltevp(points),                                                           &
                        ! Evaporation of melting aggregates
  psfall(points),                                                             &
                     ! Fall-out of snow aggregates

  praut(points),                                                              &
                     ! Autoconversion of cloud drops to rain
  pracw(points),                                                              &
                     ! Accretion of liq. water by rain
  prevp(points),                                                              &
                     ! Evaporation of rain
  prfall(points),                                                             &
                     ! Fall-out of rain

  plset(points),                                                              &
                     ! Droplet settling of liquid water
  plevpset(points),                                                           &
                        ! Evaporated settled droplets

  pgaut(points),                                                              &
                     ! Autoconversion of graupel from aggregates
  pgacw(points),                                                              &
                     ! Accretion of liq. water by graupel
  pgacs(points),                                                              &
                     ! Collection of snow aggregates by graupel
  pgmlt(points),                                                              &
                     ! Melting of graupel
  pgfall(points),                                                             &
                     ! Fall-out of graupel

  pifrw(points),                                                              &
                     ! Homogeneous freezing nucleation
  pifrr(points),                                                              &
                     ! Homogeneous freezing nucleation of rain
  piprm(points),                                                              &
                     ! Heterogeneous (primary) nucleation
  pidep(points),                                                              &
                     ! Deposition of vapour to ice crystals
  piacw(points),                                                              &
                     ! Accretion of liq. water by ice crystals
  piacr(points),                                                              &
                     ! Collection of rain by ice crystals
  pimlt(points),                                                              &
                     ! Melting of ice crystals
  pimltevp(points),                                                           &
                        ! Evaporation of melting ice crystals
  pifall(points),                                                             &
                     ! Fall-out of ice crystals
  frac_agg(points)
                     ! Fraction of ice that is aggregates

REAL (KIND=real_lsprec), INTENT(OUT) ::                                       &
  snowt_cry(points),                                                          &
                        ! Ice crystal precip out of layer / kg m-2 s-1
  snowt_agg(points),                                                          &
                        ! Ice agg. precip out of layer / kg m-2 s-1
  rainratet(points),                                                          &
                        ! Rain precip rate out of layer / kg m-2 s-1
  graupratet(points),                                                         &
                        ! Graupel precip rate out of layer/ kg m-2 s-1
  qcf_agg(points),                                                            &
                        ! Ice aggregate mixing ratio / kg kg-1
  qcf_cry(points),                                                            &
                        ! Ice crystal mixing ratio / kg kg-1
  qcf_tot(points),                                                            &
                        ! Total ice content / kg kg-1
  frac_cry_dep(points),                                                       &
                        ! Fraction of supersaturation that can be
                        ! removed by crystals
  frac_agg_dep(points),                                                       &
                        ! Fraction of supersaturation that can be
                        !removed by aggregates
  qrain(points),                                                              &
                        ! Rain water content / kg kg-1
  qs(points),                                                                 &
                        ! Saturated humidity wrt ice / kg kg-1
  qsl(points),                                                                &
                        ! Saturated humidity wrt liquid / kg kg-1
  rho(points),                                                                &
                     ! Air density / kg m-3
  rhor(points),                                                               &
                     ! 1 / air density / m kg-1
  esi(points),                                                                &
                     ! Vapour pressure wrt ice / N m-2
  esw(points),                                                                &
                     ! Vapour pressure wrt liquid / N m-2
  lheat_correc_liq(points),                                                   &
                               ! Latent heat correction for
                               ! liquid (no units)
  lheat_correc_ice(points),                                                   &
                               ! Latent heat correction for
                               ! ice (no units)
  dhi(points),                                                                &
                        ! microphysics timestep / deltaz / s m-1
  dhir(points),                                                               &
                        ! deltaz / microphysic timestep / m s-1
  corr(points),                                                               &
                        ! Fall speed correction factor (no units)
  corr2(points),                                                              &
                        ! Diffusivity correction factor (no units)
  rocor(points),                                                              &
                        ! sqrt(rho*corr*corr2) (no units)
  tcg(points),                                                                &
                        ! Temperature dependent aggregate PSD
                        ! intercept factor (no units)
  tcgi(points),                                                               &
                        ! 1 / tcg (no units)
  tcgc(points),                                                               &
                        ! Temperature dependent crystal PSD
                        ! intercept factor (no units)
  tcgci(points),                                                              &
                        ! 1 / tcgc (no units)
  tcgg(points),                                                               &
                        ! Temperature dependent graupel PSD
                        ! intercept factor (no units)
  tcggi(points),                                                              &
                        ! 1 / tcgg (no units)

  cf_transfer_diag(points),                                                   &
                        ! Diagnostic for change in cloud frac
  cfl_transfer_diag(points),                                                  &
                        ! Diagnostic for change in liq cloud frac
  cff_transfer_diag(points),                                                  &
                        ! Diagnostic for change in ice cloud frac
  rf_transfer_diag(points),                                                   &
                          ! Diagnostic for change in rain fraction
  vm_cry(points),                                                             &
                   ! Mass-weighted mean fallspeed crystal parameters
  vm_agg(points),                                                             &
                   ! Mass-weighted mean fallspeed aggregate parameters
  vm_used(points)
                   ! Selected mass-weighted mean fallspeed


LOGICAL, INTENT(OUT) ::                                                       &
  l_use_agg_vt(points)
                 ! At each point defines which branch of the
! fallspeed relation is used. If .true. then use aggregate fallspeed
! parameters; else use crystal parameters.

! Local Variables

REAL (KIND=real_lsprec) ::                                                    &
  corr_in(points),                                                            &
  corr2_in(points),                                                           &
  rocor_in(points),                                                           &
  frac_agg_out(points),                                                       &
  tcg_in(points),                                                             &
  tcgc_in(points),                                                            &
  tcgi_in(points),                                                            &
  tcgci_in(points)

INTEGER ::                                                                    &
  i              ! Loop counter

REAL (KIND=real_lsprec), PARAMETER ::                                         &
  rho1 = one
                                   ! Ref. air density / kg m-3

! Variables for use when different fallspeed relations are used
! for crystals and aggregates with the generic psd
REAL (KIND=real_lsprec) ::                                                    &
  m_bic_dic(points),                                                          &
         ! Psd moment that determines vertical ice mass flux when
! crystal fallspeed parameters are used
    m_bi_di(points)
           ! Psd moment that determines vertical ice mass flux when
! aggregate fallspeed parameters are used

REAL (KIND=real_lsprec) ::                                                    &
  cfice(points),                                                              &
         ! ice cloud fraction
  cficei(points)
         ! reciprocal ice cloud fraction

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_INIT'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

DO i = 1, points

    !-----------------------------------------------
    ! Set fluxes to zero
    !-----------------------------------------------
  snowt_cry(i)  = zero
  snowt_agg(i)  = zero
  rainratet(i)  = zero
  graupratet(i) = zero

   !--------------------------------------------------------------
   ! Set transfer variables (only used in lsp_ice and not actually
   ! coupled to STASH diagnostics) to zero
   !--------------------------------------------------------------
  cf_transfer_diag(i)  = zero
  cfl_transfer_diag(i) = zero
  cff_transfer_diag(i) = zero
  rf_transfer_diag(i)  = zero

   !----------------------------------------------
   ! Initialise fallspeed variables
   !----------------------------------------------
  vm_cry(i) = zero
  vm_agg(i) = zero
  vm_used(i) = zero
  l_use_agg_vt(i) = .FALSE.
  cfice(i) = zero
  cficei(i) = zero

END DO

IF ( .NOT. l_fix_mphys_diags_iter ) THEN

  ! N.B. IF REMOVING THIS FIX, PLEASE ALSO REMOVE THE
  ! DIAGNOSTICS BELOW FROM THE ARGUMENT LIST OF THIS
  ! SUBROUTINE

 DO i = 1, points

    !-----------------------------------------------
    ! Initialize transfer diagnostics to zero
    !-----------------------------------------------
    pifrw(i) = zero
    pifrr(i) = zero
    piprm(i) = zero
    pidep(i) = zero
    piacw(i) = zero
    piacr(i) = zero
    pimlt(i) = zero
    pimltevp(i) = zero
    pifall(i)= zero

    psdep(i) = zero
    psaut(i) = zero
    psacw(i) = zero
    psacr(i) = zero
    psaci(i) = zero
    psmlt(i) = zero
    psmltevp(i) = zero
    psfall(i)= zero

    praut(i) = zero
    pracw(i) = zero
    prevp(i) = zero
    prfall(i)= zero

    plset(i) = zero
    plevpset(i) = zero

    pgaut(i) = zero
    pgacw(i) = zero
    pgacs(i) = zero
    pgmlt(i) = zero
    pgfall(i)= zero

  END DO

END IF ! l_fix_mphys_diags_iter

    !-----------------------------------------------
    ! Set mixing ratios of graupel to zero if not used.
    ! If not done then (with high optimisation) this can be any
    ! old rubbish from memory and you get d_qgraup_dt being
    ! calculated as nan-nan in later routines, which causes a crash
    !-----------------------------------------------
IF (.NOT. l_mcr_qgraup) THEN

  DO i = 1, points

    qgraup(i)   = zero

  END DO

END IF  ! .not. l_mcr_qgraup

IF (l_mcr_qcf2) THEN
      ! If l_mcr_qcf2 is true then there are two ice/snow
      ! prognostics (qcf and qcf2), so copy them to
      ! qcf_cry and qcf_agg

  DO i = 1, points

    qcf_cry(i) = qcf2(i)
    qcf_agg(i) = qcf(i)

  END DO

ELSE IF (.NOT. l_crystals) THEN
      ! All ice is placed in a single ice category (aggregates).

  DO i = 1, points

    frac_agg(i) = one
    qcf_cry(i)  = zero
    qcf_agg(i)  = qcf(i)
    snow_cry(i) = zero

  END DO

ELSE
      ! Split the one ice/snow
      ! prognostic (qcf) diagnostically into ice crystals
      ! (qcf_cry) and snow aggregates (qcf_agg)

  DO i = 1, points

    frac_agg(i) = -t_scaling*MAX((t(i)-cttemp(i)),zero)                       &
                  *MAX(qcf(i)*qcf0,zero)
  END DO
  CALL exp_v(points,frac_agg,frac_agg_out)
  DO i = 1, points
    frac_agg(i) = MAX(one-frac_agg_out(i) , zero)
        ! Allocate ice content to crystals and aggregates
    qcf_cry(i) = qcf(i) * (one-frac_agg(i))
    qcf_agg(i) = qcf(i) * frac_agg(i)

        ! Assume falling snow is partitioned into crystals and
        ! aggregates. Snow_agg contains total snow on input
    snow_cry(i) = snow_agg(i) * (one-frac_agg(i))
    snow_agg(i) = snow_agg(i) * frac_agg(i)

  END DO

END IF ! l_mcr_qcf2

    !-----------------------------------------------
    ! Calculate total ice content
    !-----------------------------------------------
DO i = 1, points

  qcf_tot(i) = qcf_cry(i) + qcf_agg(i)

END DO

    !-----------------------------------------------
    ! Calculate the fraction of ice in the crystals and
    ! aggregates category that is allowed to remove
    ! supersaturation
    !-----------------------------------------------
IF (l_cry_agg_dep) THEN

  IF (l_mcr_qcf2) THEN

        ! Use the partition given by the prognostic ice categories.
        ! Here we add a small amount of ice to the crystal category
        ! to ensure no problems when the ice contents are zero.

    DO i = 1, points

      frac_cry_dep(i) = (qcf_cry(i)+m0)/(MAX(qcf_tot(i),zero)+m0)
      frac_agg_dep(i) = qcf_agg(i)/(MAX(qcf_tot(i),zero)+m0)

    END DO

  ELSE  ! l_mcr_qcf2
        ! Use the diagnostic partition function

    DO i = 1, points

      frac_cry_dep(i) = one - frac_agg(i)
      frac_agg_dep(i) = frac_agg(i)

    END DO

  END IF

ELSE
      ! Set the fractions to 1 to maintain bit reproducibility

  DO i = 1, points

    frac_cry_dep(i)=one
    frac_agg_dep(i)=one

  END DO

END IF

    !-----------------------------------------------
    ! If rain is a diagnostic, convert flux (kg m-2 s-1)
    ! to mass (kg kg-1)
    !-----------------------------------------------
IF (.NOT. l_mcr_qrain) THEN

  ! Rain is a diagnostic quantity

  DO i = 1, points

    IF (rainrate(i)  >   zero) THEN

      IF (l_mr_physics) THEN

        ! Mixing ratio formulation
        qrain(i) = rainrate(i) * timestep / rhodz_dry(i)

      ELSE ! l_mr_physics

        ! Specific humidity formulation
        qrain(i) = rainrate(i) * timestep / rhodz_moist(i)

      END IF  ! l_mr_physics

    ELSE    ! rainrate > 0

      qrain(i) = zero

    END IF  ! rainrate > 0

  END DO ! points

END IF  ! .not. l_mcr_qrain

    !-----------------------------------------------
    ! Calculate saturation specific humidities
    !-----------------------------------------------

IF ( l_new_qsat_lsprec ) THEN
  IF (l_mr_physics) THEN
      ! Qsat with respect to ice
      CALL qsat_mix_new(qs,t,p,points)
      ! Qsat with respect to liquid water
      CALL qsat_wat_mix_new(qsl,t,p,points)
  ELSE
      ! Qsat with respect to ice
      CALL qsat_new(qs,t,p,points)
      ! Qsat with respect to liquid water
      CALL qsat_wat_new(qsl,t,p,points)
  END IF
ELSE
  ! Qsat with respect to ice
! DEPENDS ON: qsat_mix
  CALL qsat_mix(qs,t,p,points,l_mr_physics)
  ! Qsat with respect to liquid water
! DEPENDS ON: qsat_wat_mix
  CALL qsat_wat_mix(qsl,t,p,points,l_mr_physics)
END IF

!When working in 32-bit, a Cray compiler bug breaks PROC comparability. 
!Conditionally using NOVECTOR makes this go away. Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAY_FORTRAN) && (CRAY_FORTRAN <8004000)
!DIR$ NOVECTOR
#endif
DO i = 1, points

    !-----------------------------------------------
    ! Calculate saturation vapour pressures
    !-----------------------------------------------
  esi(i) = qs (i) * p(i) * recip_epsilon
  esw(i) = qsl(i) * p(i) * recip_epsilon

    !-----------------------------------------------
    ! Calculate density of air
    !-----------------------------------------------
  IF (l_mr_physics) THEN

    ! rho is the dry density
    rho(i) = rhodz_dry(i) / deltaz(i)

  ELSE

    ! rho is the moist density
    rho(i) = rhodz_moist(i) / deltaz(i)

  END IF  ! l_mr_physics

      ! Calculate the inverse of the air density
  rhor(i) = one / rho(i)

      !-----------------------------------------------
      ! Estimate latent heat correction to rate of evaporation
      !-----------------------------------------------
  lheat_correc_liq(i) = one/(one+repsilon*lc**2*qsl(i)                        &
                           /(cp*r*t(i)**2))
  lheat_correc_ice(i) = one/(one+repsilon*(lc+lf)**2*qs(i)                    &
                           /(cp*r*t(i)**2))

      !-----------------------------------------------
      ! Calculate CFL timestep divided by level separation
      !-----------------------------------------------

      ! Use the formulation based on the heights of the levels
  dhi(i) = timestep_mp/deltaz(i)

      ! Calculate the inverse
  dhir(i)       = one/dhi(i)

END DO ! points

    !-----------------------------------------------
    ! Correction factors due to air density and temperature
    !-----------------------------------------------

DO i = 1, points

      ! Correction of fall speeds
  IF (l_mr_physics) THEN

    corr_in(i) = rho1*deltaz(i) / rhodz_moist(i)

  ELSE

    corr_in(i) = rho1*rhor(i)

  END IF

      ! Correction factor in viscosity etc. due to temperature

  corr2_in(i) = (t(i) * (one / zerodegc))

END DO

CALL powr_v( points, corr_in,  0.4_real_lsprec,   corr  )
CALL powr_v( points, corr2_in, cpwr , corr2 )

! Determine which set of fallspeed parameters
! to use in process rate calculations, based on
! mass-weighted mean fallspeed

IF ( l_diff_icevt .AND. l_psd ) THEN
  DO i = 1, points

    ! -------------------------------------------------------------
    ! Check that ice cloud fraction is sensible. This should mimic
    ! exactly what is done in lsp_ice prior to the process
    ! rate calculations.
    ! -------------------------------------------------------------

    cfice(i)  = MAX( cficekeep(i), 0.001_real_lsprec )
    cficei(i) = one / cfice(i)

  END DO

  !-------------------------------------------------------------
  ! Calculate the psd moments required for crystal and
  ! aggregate fallspeeds
  !-------------------------------------------------------------


  CALL lsp_moments(points,rho,t,qcf,cficei,cx(182),m_bic_dic)
    ! ice mass flux moment with crystal parameters

  CALL lsp_moments(points,rho,t,qcf,cficei,cx(82),m_bi_di)
    ! ice mass flux moment with aggregate parameters

  DO i = 1, points

     !-------------------------------------------------------------
     ! Calculate mass-weighted mean ice fallspeeds for the
     ! two vt-D relations
     !-------------------------------------------------------------
    IF (qcf(i) >  m0) THEN
        ! Only calculate a fallspeed if qcf exceeds minimum value
        ! set by the nucleation mass. This is how the
        ! fallspeeds are calculated in lsp_fall
      vm_cry(i) = constp(182)* corr(i) * m_bic_dic(i)                         &
                / (rho(i) * qcf(i) * cficei(i))
      vm_agg(i) = constp(82) * corr(i) * m_bi_di(i)                           &
                / (rho(i) * qcf(i) * cficei(i))
      !-------------------------------------------------------------
      ! Determine which set of fallspeed parameters to use
      ! in process rate calculations
      !-------------------------------------------------------------
      IF ( vm_cry(i)  <=  vm_agg(i) ) THEN
        l_use_agg_vt(i) = .FALSE.
                   ! Use crystal vt-D parameters
      ELSE
        l_use_agg_vt(i) = .TRUE.
                   ! Use aggregate vt-D parameters
      END IF
    END IF
    IF (qcf(i)  <=   m0) THEN
        ! Use crystal fallspeeds for very small qcf values
      l_use_agg_vt(i) = .FALSE.
    END IF

  END DO

END IF !  l_diff_icevt .AND. l_psd

DO i = 1, points

  corr2(i) = corr2(i)*( tcor1 /(t(i)+ tcor2))

  ! Combined correction factor

  IF (l_mr_physics) THEN

    rocor_in(i) = rhodz_moist(i) / deltaz(i) * corr(i) * corr2(i)

  ELSE

    rocor_in(i) = rho(i)*corr(i)*corr2(i)

  END IF

END DO

CALL powr_v( points, rocor_in, half, rocor )

    !-----------------------------------------------
    ! Calculate ice particle size distributions
    !-----------------------------------------------

DO i = 1, points

  ! Calculate a temperature factor for N0aggregates
  tcg_in(i)  = -cx(32)*MAX(t(i)-zerodegc,t_agg_min)

  ! Define inverse of TCG values
  tcgi_in(i)  = -tcg_in(i)

END DO

CALL exp_v( points, tcg_in,   tcg   )
CALL exp_v( points, tcgi_in,  tcgi  )

IF ( l_crystals ) THEN

  DO i = 1, points
    ! Calculate a temperature factor for N0crystals
    tcgc_in(i) = -cx(12)*MAX(t(i)-zerodegc,t_agg_min)

    ! Define inverse of TCGC values
    tcgci_in(i) = -tcgc_in(i)

  END DO

  CALL exp_v( points, tcgc_in,  tcgc  )
  CALL exp_v( points, tcgci_in, tcgci )

END IF ! l_crystals

    !-----------------------------------------------
    ! Calculate graupel size distributions
    !-----------------------------------------------
DO i = 1, points

  tcgg(i)  = one
  tcggi(i) = one ! (Equivalent to 1.0/tcgg, so will always be 1.0)

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_init
END MODULE lsp_init_mod
