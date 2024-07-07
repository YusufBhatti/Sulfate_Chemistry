! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Precipitation microphysics calculations.
! Subroutine Interface:
MODULE lsp_ice_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='LSP_ICE_MOD'

CONTAINS

SUBROUTINE lsp_ice(                                                           &
  p,rhodz,deltaz,rhodz_dry,rhodz_moist,                                       &
  points, rhcpt, land_fract, hmteff,zb,                                       &
  qcf,qcl,q,qcf2,qrain,qgraup, n_drop_tpr, n_drop_out,                        &
  rainrate, vf_rain, snow_agg, vf_agg,                                        &
  snow_cry, vf_cry, grauprate, vf_graup, droplet_flux,                        &
  frac_ice_above, frac_agg, cttemp, rainfrac, rainfrac_impr,                  &
  t,cfkeep,cfliqkeep,cficekeep,bland,                                         &
  psdep,psaut,psacw,psacr,psaci,psmlt,psmltevp,                               &
  praut,pracw,prevp,                                                          &
  pgaut,pgacw,pgacs,pgmlt,                                                    &
  pifrw,pifrr,piprm,pidep,piacw,piacr,pimlt,pimltevp,                         &
  pifall,psfall,prfall,pgfall,plset,plevpset,                                 &
  dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_l, dbz_r,                                &
  sfwater, sfrain, sfsnow,                                                    &
  niters_mp, uk, vk, ukp1, vkp1,                                              &
  r_theta_levels_c, fv_cos_theta_latitude_c,                                  &
  r_theta_surf_c,                                                             &
  f_arr1, f_arr2, f_arr3,                                                     &
  vm_cry, vm_agg, vtbranch_flag, vm_used                                      &
      )

!Use in reals in lsprec precision, both microphysics related and general atmos
USE lsprec_mod,           ONLY: lcrcp, lfrcp, timestep, m0, t_scaling,        &
                                qcf0, timestep_mp, zero, one

USE um_types,             ONLY: real_lsprec

  ! General atmosphere modules- logicals and integers
USE gen_phys_inputs_mod,  ONLY: l_mr_physics

! Microphysics modules
USE mphys_inputs_mod,     ONLY: l_psd, l_sr2graup,              &
                                l_mcr_qcf2, l_mcr_qrain,                      &
                                l_mcr_qgraup, l_shape_rime,                   &
                                l_orograin,l_orogrime,                        &
                                not_generic_size_dist
USE mphys_bypass_mod,     ONLY: l_crystals, l_ref_diag

! Extra USE statements for Seeder Feeder scheme
USE water_constants_mod,   ONLY: lc, lf
USE mphys_ice_mod,         ONLY: qcfmin
USE planet_constants_mod,  ONLY: rv

! Dr Hook Modules
USE yomhook,              ONLY: lhook, dr_hook
USE parkind1,             ONLY: jprb, jpim

USE lsp_accretion_mod,    ONLY: lsp_accretion
USE lsp_autoc_mod,        ONLY: lsp_autoc
USE lsp_capture_mod,      ONLY: lsp_capture
USE lsp_collection_mod,   ONLY: lsp_collection
USE lsp_deposition_mod,   ONLY: lsp_deposition
USE lsp_evap_mod,         ONLY: lsp_evap
USE lsp_evap_snow_mod,    ONLY: lsp_evap_snow
USE lsp_fall_mod,         ONLY: lsp_fall
USE lsp_graup_autoc_mod,  ONLY: lsp_graup_autoc
USE lsp_init_mod,         ONLY: lsp_init
USE lsp_melting_mod,      ONLY: lsp_melting
USE lsp_nucleation_mod,   ONLY: lsp_nucleation
USE lsp_riming_mod,       ONLY: lsp_riming, lsp_riming_sphere
USE lsp_settle_mod,       ONLY: lsp_settle
USE lsp_snow_autoc_mod,   ONLY: lsp_snow_autoc
USE lsp_subgrid_mod,      ONLY: lsp_subgrid
USE lsp_tidy_mod,         ONLY: lsp_tidy
USE mphys_reflec_mod,     ONLY: mphys_reflec
USE lsp_orogwater_mod,    ONLY: lsp_orogwater

IMPLICIT NONE

! Description:
!   Updates ice, liquid and vapour contents, temperature and
!   cloud fractions as a result of microphysical processes.

! Method:
!   Calculates transfers of water between vapour, ice/snow,
!   cloud liquid water, rain and graupel.

!   Processes included are:
!   - Fall of hydrometeor into and out of the layer (sedimentation)
!   - Homogenous and heterogenous nucleation of ice;
!   - Deposition and sublimation of ice/snow;
!   - Autoconversion of ice->snow, snow->graupel, liquid->rain
!   - Collection processes
!   - Melting and evaporation

!   This is described in Unified Model Documentation Paper 26.

!   There are a number of different options for prognostic
!   hydrometeor variables. Each is independent of the others.
!   - Second prognostic cloud ice variables
!      Active if l_mcr_qcf2=.True. (in CNTLATM namelist)
!      The code supports the use of a second cloud ice prognostic
!      variable so that both cloud ice aggregates (QCF/QCF_AGG)
!      and cloud ice pristine crystals (QCF2/QCF_CRY) can be
!      represented and advected separately.
!      If False, then there is a diagnostic split within each level
!      at each timestep into ice crystals and snow aggregates.
!   - Prognostic rain
!      Active if l_mcr_qrain=.True. (in CNTLATM namelist)
!      If False, then rain is treated diagnostically.
!   - Prognostic graupel
!      Active if l_mcr_qgraup=.True. (in CNTLATM namelist)
!      If False, then graupel is not represented.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.

! Declarations:

! Subroutine arguments

INTEGER, INTENT(IN) ::                                                        &
  points,                                                                     &
!         Number of points to be processed.
  niters_mp
!         Number of iterations of microphysics


LOGICAL, INTENT(IN) ::                                                        &
  bland(points)
!         Land/sea mask

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  p(points),                                                                  &
!         Air pressure at this level (Pa).
    rhodz(points),                                                            &
!         Air mass p.u.a. in this layer (kg per sq m).
    deltaz(points),                                                           &
                           ! Thickness of layer (m)
    rhodz_dry(points),                                                        &
                           ! Dry air density * layer thickness (kg m-2)
    rhodz_moist(points),                                                      &
                           ! Moist air density * layer thick.  (kg m-2)
    land_fract(points),                                                       &
!         Land fraction
    hmteff(points),                                                           &
!         Effective mountain height forcing ascent for seeder feeder
    zb(points),                                                               &
!         Blocked layer depth for seeder feeder
    rhcpt(points),                                                            &
!         Critical relative humidity of all points for cloud formation.
    n_drop_tpr(points),                                                       &
!         droplet number (/ m 3) for autoconversion and settling
    f_arr1(points),                                                           &
    f_arr2(points),                                                           &
    f_arr3(points)

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  cfliqkeep(points),                                                          &
!         Liquid cloud fraction in this layer (no units).
    cficekeep(points),                                                        &
!         Frozen cloud fraction in this layer (no units).
    cfkeep(points),                                                           &
!         Total cloud fraction in this layer (no units).
    q(points),                                                                &
!         Specific humidity at this level (kg water per kg air).
    qcf(points),                                                              &
!         Cloud ice (kg water per kg air).
    qcl(points),                                                              &
!         Cloud liquid water (kg water per kg air).
    qcf2(points),                                                             &
!         Second cloud ice (kg water per kg air).
    qrain(points),                                                            &
!         Rain (kg water per kg air).
    qgraup(points),                                                           &
!         Graupel (kg water per kg air).
    t(points)
!         Temperature at this level (K).

!     Hydrometeor flux between layers (kg m-2 s-1)
!     On input: Hydrometeor flux entering this layer from above.
!     On output: Hydrometeor flux leaving this layer.
!     Note: If only one ice prognostic is active (l_mcr_qcf2=.F.)
!     then SNOW_AGG contains all the ice/snow and SNOW_CRY is zero.

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  snow_agg(points),                                                           &
                            ! snow aggregates
  snow_cry(points),                                                           &
                            ! ice crystals
  rainrate(points),                                                           &
                            ! rain
  grauprate(points),                                                          &
                            ! graupel
  droplet_flux(points)
                            ! droplets

!     Hydrometeor mass-weighted fall velocities (m s-1)
!     On input: Fall velocity of hydrometeor entering layer.
!     On Output: Fall velocity of hydrometeor leaving layer.
!     Note: If only one ice prognostic is active, then only
!     VF_AGG is used.

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  vf_agg(points),                                                             &
                         ! snow aggregates
  vf_cry(points),                                                             &
                         ! ice crystals
  vf_rain(points),                                                            &
                         ! rain
  vf_graup(points)   ! graupel

!     Cloud/precipitation sub-grid fraction variables

REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  cttemp(points),                                                             &
!         Ice cloud top temperature (K)
    rainfrac(points),                                                         &
!         Rain fraction (no units)
    frac_ice_above(points),                                                   &
!         Fraction of ice in layer above (no units)
    frac_agg(points)
!         Aggregate fraction

REAL (KIND=real_lsprec), INTENT(IN) :: rainfrac_impr(points)

! Microphysical process rate diagnostics
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
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
  psmltevp(points)  ! Evaporation of melting aggregates
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  praut(points),                                                              &
                     ! Autoconversion of cloud drops to rain
  pracw(points),                                                              &
                     ! Accretion of liq. water by rain
  prevp(points)  ! Evaporation of rain
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  pgaut(points),                                                              &
                     ! Autoconversion of graupel from aggregates
  pgacw(points),                                                              &
                     ! Accretion of liq. water by graupel
  pgacs(points),                                                              &
                     ! Collection of snow aggregates by graupel
  pgmlt(points)  ! Melting of graupel
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
  pifrw(points),                                                              &
                     ! Homogeneous freezing nucleation
  pifrr(points),                                                              &
                     ! Homogeneous freezing of rain
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
  pimltevp(points)  ! Evaporation of melting ice crystals
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
 pifall(points),                                                              &
                     ! Sedimentation of ice crystals
 psfall(points),                                                              &
                     ! Sedimentation of aggregates
 prfall(points),                                                              &
                     ! Sedimentation of rain
 pgfall(points) ! Sedimentation of graupel
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
 plset(points),                                                               &
                     ! Droplet settling of liquid water
 plevpset(points) ! Evaporated settled droplets

! Seeder feeder water and precip rate diagnostics
REAL(KIND=real_lsprec), INTENT(INOUT) ::                                      &
 sfwater(points),                                                             &
    ! Subgrid orographic water mixing ratio
 sfrain(points),                                                              &
    ! Subgrid orographic rain creation rate
 sfsnow(points)
    ! Subgrid orographic snow creation rate

REAL (KIND=real_lsprec), INTENT(IN) ::                                        &
  uk  (points),                                                               &
                     ! U wind at level k
  vk  (points),                                                               &
                     ! V wind at level k
  ukp1(points),                                                               &
                     ! U wind at level k+1
  vkp1(points),                                                               &
                     ! V wind at level k+1
  r_theta_levels_c(points),                                                   &
                     ! Distance from centre of the Earth
  r_theta_surf_c(points),                                                     &
                     ! ...and near surface (k=0)
  fv_cos_theta_latitude_c(points)
                     ! Finite volume cosine of latitude.

REAL (KIND=real_lsprec), INTENT(INOUT) :: n_drop_out(points)
                    ! output droplet number from autoconversion

! Variables for defining different fallspeed relations for crystals
! and aggregates when the generic psd is used and l_diff_icevt = .true.
REAL (KIND=real_lsprec), INTENT(INOUT) ::                                     &
 vm_agg(points),                                                              &
    ! mass-weighted mean fallspeed using aggregate fallspeed relation
 vm_cry(points),                                                              &
    ! mass-weighted mean fallspeed using crystal fallspeed relation
 vtbranch_flag(points),                                                       &
    ! 0=crystal vt-D relation; 1=aggregate vt-D relation
 vm_used(points)
    ! mass-weigthed fallspeed used

REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_tot(points)
    ! Total reflectivity (dBZ)
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_g(points)
    ! Graupel reflectivity (dBZ)
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_i(points)
    ! Ice Agg. reflectivity (dBZ)
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_i2(points)
    ! Ice Cry. reflectivity (dBZ)
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_l(points)
    ! Cloud liquid reflectivity (dBZ)
REAL (KIND=real_lsprec), INTENT(OUT) :: dbz_r(points)
    ! Rain reflectivity (dBZ)

!  Local scalars and dynamic arrays
INTEGER ::  i  !  Loop counter (horizontal field index).

REAL (KIND=real_lsprec) :: lsrcp

REAL (KIND=real_lsprec) ::                                                    &
  qs(points),                                                                 &
!         Saturated sp humidity for (T,p) in layer (kg kg-1)
    qsl(points)
!         Saturated sp humidity for (T,p) in layer
!         wrt water at all temps (kg kg-1)

!     Cumulative fall out of hydrometeor within iterations (kg m-2 s-1)
REAL (KIND=real_lsprec) ::                                                    &
  snowt_agg(points),                                                          &
  snowt_cry(points),                                                          &
  rainratet(points),                                                          &
  graupratet(points)

! Subgrid orographic water MR
REAL (KIND=real_lsprec) ::                                                    &
        ql_orog(points),                                                      &
        cf_orog(points),                                                      &
        ql_orog1(points),                                                     &
        t1a(points),                                                          &
        t1b(points),                                                          &
        qrain1a(points),                                                      &
        qrain1b(points),                                                      &
        qsnow1a(points),                                                      &
        qsnow1b(points),                                                      &
        pracw1a(points),                                                      &
        pracw1b(points),                                                      &
        psacw1a(points),                                                      &
        psacw1b(points),                                                      &
        dqrain(points),                                                       &
        dqsnow(points),                                                       &
        rain_liq_orog(points),                                                &
        rain_mix_orog(points),                                                &
        area_mix_orog(points),                                                &
        area_liq_orog(points),                                                &
        altitude(points)

LOGICAL ::  l_enh_rime(points) !whether riming is enhanced


REAL (KIND=real_lsprec) ::                                                    &
  rho(points),                                                                &
!         Density of air in the layer (kg m-3).
    rhor(points),                                                             &
!         1.0/RHO to speed up calculations (kg-1 m3).
    esi(points),                                                              &
!         saturation vapour pressure (wrt ice below zero Celsius)(Pa)
    esw(points)
!         saturation vapour pressure (wrt water at all temperatures)(Pa)

REAL (KIND=real_lsprec) ::                                                    &
  cfice(points),                                                              &
!         fraction of ice inferred for the microphysics (no units).
    cficei(points),                                                           &
!         inverse of CFICE (no units)
    cfliq(points),                                                            &
!         liquid cloud fraction for the microphysics
    cf(points),                                                               &
!         total cloud fraction for the microphysics (no units)
    dhi(points),                                                              &
!         CFL limit (s m-1)
    dhir(points)
!         1.0/DHI (m s-1)

REAL (KIND=real_lsprec) ::                                                    &
    corr(points),                                                             &
!         density correction for fall speed (no units)
    rocor(points)
!         density correction for fall speed (no units)

REAL (KIND=real_lsprec) ::                                                    &
    corr2(points),                                                            &
!         Temperature correction of viscosity etc. (no units)
    tcg(points),                                                              &
!         Temperature Factor for aggregate size distribution (no units)
    tcgi(points),                                                             &
!         Inverse of TCG (no units)
    tcgc(points),                                                             &
!         Temperature Factor for crystal size distribution (no units)
    tcgci(points),                                                            &
!         Inverse of TCGC (no units)
    tcgg(points),                                                             &
!         Temperature Factor for graupel size distribution (no units)
    tcggi(points)
!         Inverse of TCGC (no units)

REAL (KIND=real_lsprec) ::                                                    &
  area_liq(points),                                                           &
!         Liquid only area of gridbox (no units)
    area_mix(points),                                                         &
!         Mixed phase area of gridbox (no units)
    area_ice(points),                                                         &
!         Ice only area of gridbox (no units)
    area_clear(points),                                                       &
!         Cloud free area of gridbox (no units)
    areamix_over_cfliq(points),                                               &
!         area_mix / cfliq (no units) - for perturbation sensitivity.
    rain_liq(points),                                                         &
!         Overlap fraction of gridbox between rain and liquid cloud
    rain_mix(points),                                                         &
!         Overlap fraction of gridbox between rain and mixed phase cloud
    rain_ice(points),                                                         &
!         Overlap fraction of gridbox between rain and ice cloud
    rain_clear(points),                                                       &
!         Overlap fraction of gridbox between rain and no cloud
    q_ice(points),                                                            &
!         Vapour content in the ice only part of the grid box (kg kg-1)
    q_clear(points),                                                          &
!         Vapour content in the cloud free part of the grid box(kg kg-1)
    qcf_agg(points),                                                          &
!         QCF in the form of aggregates (kg kg-1)
    qcf_cry(points)
!         QCF in the form of crystals (kg kg-1)

REAL (KIND=real_lsprec) ::                                                    &
    qcf_tot(points),                                                          &
!         Total amount of ice (crystals+aggregates)
    frac_cry_dep(points),                                                     &
!         Fraction of supersaturation that can be removed by crystals
    frac_agg_dep(points)
!         Fraction of supersaturation that can be removed by aggregates

REAL (KIND=real_lsprec) ::                                                    &
  lheat_correc_liq(points),                                                   &
!         Reduction factor in evaporation limits because of latent heat
    lheat_correc_ice(points),                                                 &
!         Reduction factor in evaporation limits because of latent heat
    q_ice_1(points), q_ice_2(points),                                         &
    area_ice_1(points), area_ice_2(points),                                   &
!         Subgrid splitting for deposition term
    qcft(points)   ! Holds total ice content qcf_cry + qcf_agg


! Cloud and rain fraction transfer rate diagnostics
REAL (KIND=real_lsprec) ::                                                    &
 cf_transfer_diag(points),                                                    &
                                ! Dummy to receive cf increment
 cfl_transfer_diag(points),                                                   &
                                ! Dummy to receive cfl increments
 cff_transfer_diag(points),                                                   &
                                ! Dummy to receive cff increments
 rf_transfer_diag(points)   ! Dummy to receive rainfrac increment

REAL (KIND=real_lsprec) :: one_over_tsi        ! 1.0/(timestep_mp*iterations)

! Variables for calls to lsp_collection
LOGICAL ::                                                                    &
 l_use_area,                                                                  &
!       Use the ice partition to calculate transfer rates rather than
!       assuming a uniform distribution through the gridbox
   l_no_t_check
!       Do not check that the temperature is below zero degrees Celsius

LOGICAL ::                                                                    &
 l_use_agg_vt(points)              ! For each point defines which
! branch of the fallspeed relation is used. .true.=aggregate params.
! Set in lsp_init

INTEGER ::                                                                    &
 ice_type1, ice_type2
!       Category of ice involved in collision process:
!       0 - crystals; 1 - aggregates; 3 - graupel.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='LSP_ICE'

!- End of header

! ======================================================================
!       Initialize variables and calculate microphysical quantities
! ======================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL lsp_init(points, timestep,                                               &
              timestep_mp,                                                    &
              t, p, cttemp,                                                   &
              deltaz, rhodz,rhodz_dry,rhodz_moist, rho,rhor,                  &
              dhi, dhir,                                                      &
              q, qcl, qcf, qcf2, qrain, qgraup, rainrate,                     &
              qcf_agg, qcf_cry,                                               &
              qcf_tot, frac_agg, frac_cry_dep, frac_agg_dep,                  &
              qs, qsl, esi, esw,                                              &
              psdep, psaut, psacw, psacr, psaci, psmlt,                       &
              psmltevp, psfall, pifrw, pifrr, piprm, pidep,                   &
              piacw, piacr, pimlt, pimltevp, pifall,                          &
              praut, pracw, prevp, prfall, plset, plevpset,                   &
              pgaut, pgacw, pgacs, pgmlt, pgfall,                             &
              cf_transfer_diag, cfl_transfer_diag,                            &
              cff_transfer_diag, rf_transfer_diag,                            &
              snow_cry, snow_agg, snowt_cry, snowt_agg,                       &
              rainratet, graupratet,                                          &
              lheat_correc_liq, lheat_correc_ice,                             &
              corr, corr2, rocor,                                             &
              tcg, tcgi, tcgc, tcgci, tcgg, tcggi,                            &
              cficekeep,                                                      &
              vm_cry, vm_agg, l_use_agg_vt, vm_used,                          &
              niters_mp )

! Calculate one_over_tsi once here rather than multiple times later
one_over_tsi = one/(timestep_mp*niters_mp)

! Latent heat of sublimation / cp (K).
lsrcp=lcrcp+lfrcp

! ======================================================================
!       Start iterations
! ======================================================================

  DO i = 1, points
         !------------------------------------------------------------
         ! Set fallspeed branch flag for diagnostics
         !------------------------------------------------------------
    IF ( l_use_agg_vt(i) ) THEN
      vtbranch_flag(i) = one   ! Fallspeed setting for aggregates
    ELSE
        vtbranch_flag(i) = zero   ! Fallspeed setting for crystals
    END IF

    ! -------------------------------------------------------------
    ! Check that ice cloud fraction is sensible.
    ! -------------------------------------------------------------

    cfice(i)  = MAX( cficekeep(i), 0.001_real_lsprec )
    cficei(i) = one / cfice(i)
    cf(i)     = cfkeep(i)
    cfliq(i)  = cfliqkeep(i)
    cf(i)     = MIN( MAX(cf(i),cfice(i)) ,(cfice(i)+cfliq(i)) )

         ! -------------------------------------------------------------
         ! Calculate overlaps of liquid, ice and rain fractions
         ! -------------------------------------------------------------

    area_ice(i)   = MAX(cf(i)-cfliq(i),zero)
    area_clear(i) = MAX(one-cf(i),zero)

         ! -------------------------------------------------------------
         ! Update ice cloud top temperature if no ice falling in
         ! -------------------------------------------------------------

    IF (snow_cry(i)+snow_agg(i)  <=  zero) THEN
      cttemp(i) = t(i)
    END IF

  END DO

  ! ======================================================================
  !        Droplet settling
  ! ======================================================================
  CALL lsp_settle(points, one_over_tsi,                                       &
                  q, qcl, t, droplet_flux, bland,                             &
                  cfliq, rho, rhor, corr2, lcrcp,                             &
                  dhi, dhir,                                                  &
                  n_drop_tpr,                                                 &
                  plset, plevpset                                             &
                 )

  ! ======================================================================
  !        Sedimentation of ice, rain and graupel
  ! ======================================================================

  CALL lsp_fall(points, timestep_mp,                                        &
                qcf_cry, qcf_agg, frac_agg, qrain, qgraup, t,               &
                snow_agg, snow_cry, rainrate, grauprate,                    &
                snowt_agg, snowt_cry, rainratet, graupratet,                &
                vf_agg, vf_cry, vf_rain, vf_graup,                          &
                area_clear, area_ice, cfice, cficei,                        &
                frac_ice_above, cfkeep, cfliqkeep, cficekeep,               &
                rho, rhor, tcgi, tcgci,                                     &
                corr, dhi, dhir, rainfrac,                                  &
                pifall, psfall, prfall, pgfall,                             &
                1, one_over_tsi,                                            &
                cf_transfer_diag, cff_transfer_diag,                        &
                uk, vk, ukp1, vkp1,                                         &
                r_theta_levels_c, fv_cos_theta_latitude_c,                  &
                l_use_agg_vt,                                               &
                vm_used                                                     &
               )

  ! ======================================================================
  !        Subgrid-scale set-up calculations and tidy-ups
  ! ======================================================================
  CALL lsp_subgrid(points,                                                  &
                  q, qcf_cry, qcf_agg, qcf_tot, t,                          &
                  qsl, qs,                                                  &
                  q_ice, q_clear, q_ice_1, q_ice_2,                         &
                  area_liq,area_mix,area_ice,area_clear,                    &
                  area_ice_1, area_ice_2,                                   &
                  areamix_over_cfliq,                                       &
                  rain_liq,rain_mix,rain_ice,rain_clear,                    &
                  cf, cfliq, cfice, cficei,                                 &
                  frac_ice_above,                                           &
                  cfkeep, cficekeep, rainfrac,                              &
                  lsrcp, rhcpt                                              &
                 )

  ! ======================================================================
  !        HOMOGENEOUS (PIFRW) AND HETEROGENEOUS NUCLEATION (PIPRM)
  ! ======================================================================
  IF (l_crystals) THEN
         ! Call nucleation with the crystals ice category (qcf_cry)
    CALL lsp_nucleation(points,                                               &
                  q, qcl, qrain, qcf_cry, t,                                  &
                  qs, qsl,                                                    &
                  cfliq, cfice,                                               &
                  area_liq, area_mix,                                         &
                  cfkeep, cfliqkeep, cficekeep, rainfrac,                     &
                  rain_liq, rain_mix,                                         &
                  rhor, lheat_correc_ice,                                     &
                  lfrcp, lsrcp,                                               &
                  piprm, pifrw, pifrr, one_over_tsi,                          &
                  cf_transfer_diag, cfl_transfer_diag,                        &
                  cff_transfer_diag, rf_transfer_diag                         &
                 )
  ELSE
         ! Call nucleation with the only ice category (qcf_agg)
    CALL lsp_nucleation(points,                                               &
                  q, qcl, qrain, qcf_agg, t,                                  &
                  qs, qsl,                                                    &
                  cfliq, cfice,                                               &
                  area_liq, area_mix,                                         &
                  cfkeep, cfliqkeep, cficekeep, rainfrac,                     &
                  rain_liq, rain_mix,                                         &
                  rhor, lheat_correc_ice,                                     &
                  lfrcp, lsrcp,                                               &
                  piprm, pifrw, pifrr, one_over_tsi,                          &
                  cf_transfer_diag, cfl_transfer_diag,                        &
                  cff_transfer_diag, rf_transfer_diag                         &
                 )
  END IF  ! l_crystals

  ! ======================================================================
  !             DEPOSITION/SUBLIMATION OF ICE CRYSTALS (PIDEP)
  ! ======================================================================
  IF (l_crystals) THEN

    DO i = 1, points
         ! Calculate total ice content
      qcft(i) = qcf_cry(i) + qcf_agg(i)
    END DO

    CALL lsp_deposition(points, timestep_mp,                                  &
                    q, qcl, qcf_cry, qcft, t, p, frac_cry_dep,                &
                    q_ice_1, q_ice_2,                                         &
                    area_ice_1, area_ice_2,                                   &
                    esi, qs, qsl,                                             &
                    area_mix, cfliq, cfice, cficei,                           &
                    areamix_over_cfliq,                                       &
                    cfkeep, cfliqkeep, cficekeep,                             &
                    rho, tcgc, tcgci,                                         &
                    corr2, rocor, lheat_correc_ice,                           &
                    lfrcp, lsrcp, 0,                                          &
                    not_generic_size_dist,                                    &
                    pidep, one_over_tsi,                                      &
                    cf_transfer_diag, cfl_transfer_diag,                      &
                    cff_transfer_diag,                                        &
                    l_use_agg_vt                                              &
                   )

  END IF  ! l_crystals
  ! ======================================================================
  !           DEPOSITION/SUBLIMATION OF SNOW AGGREGATES (PSDEP)
  ! ======================================================================
  DO i = 1, points
         ! Calculate total ice content
    qcft(i) = qcf_cry(i) + qcf_agg(i)
  END DO

  CALL lsp_deposition(points, timestep_mp,                                    &
                  q, qcl, qcf_agg, qcft, t, p, frac_agg_dep,                  &
                  q_ice_1, q_ice_2,                                           &
                  area_ice_1, area_ice_2,                                     &
                  esi, qs, qsl,                                               &
                  area_mix, cfliq, cfice, cficei,                             &
                  areamix_over_cfliq,                                         &
                  cfkeep, cfliqkeep, cficekeep,                               &
                  rho, tcg, tcgi,                                             &
                  corr2, rocor, lheat_correc_ice,                             &
                  lfrcp, lsrcp, 1,                                            &
                  l_psd,                                                      &
                  psdep, one_over_tsi,                                        &
                  cf_transfer_diag, cfl_transfer_diag,                        &
                  cff_transfer_diag,                                          &
                  l_use_agg_vt                                                &
                 )

  IF (l_mcr_qcf2) THEN  ! Note that l_mcr_qcf2==true implies
                             ! that l_crystals==true. It is *intended*
                             ! in parametrization development that
                             ! l_mcr_qcf2==true implies l_psd==false.
    ! ======================================================================
    !          AUTOCONVERSION OF ICE CRYSTALS TO AGGREGATES (PSAUT)
    ! ======================================================================
    CALL lsp_snow_autoc(points, timestep_mp,                                  &
                  qcf_cry, qcf_agg, t, cttemp,                                &
                  m0, t_scaling,                                              &
                  qcf0,                                                       &
                  psaut, one_over_tsi                                         &
                 )

    ! ======================================================================
    !           COLLECTION OF CRYSTALS BY AGGREGATES (PSACI)
    ! ======================================================================
    l_use_area=.FALSE.
    l_no_t_check=.TRUE.
    ice_type1=1  ! aggregates
    ice_type2=0  ! crystals
    CALL lsp_collection(points, timestep_mp,                                  &
                qcf_agg, qcf_cry, t,                                          &
                area_mix, area_ice, cficei,                                   &
    !                        cf, cff,
                      rho, rhor, m0, tcg, tcgi, tcgc, tcgci, corr,            &
                      ice_type1,ice_type2,                                    &
                      not_generic_size_dist,                                  &
                      l_use_area, l_no_t_check,                               &
                      psaci, one_over_tsi,                                    &
                      l_use_agg_vt                                            &
    !                        cf_transfer_diag, cff_transfer_diag
                     )

  END IF  ! l_mcr_qcf2

  ! ======================================================================
  !              RIMING OF ICE CRYSTALS BY CLOUD WATER (PIACW)
  ! ======================================================================
  IF (l_crystals) THEN
   IF (l_shape_rime) THEN
    CALL lsp_riming(points, timestep_mp,                                      &
                    qcl, qcf_cry, t,                                          &
                    area_liq, area_mix, cfliq, cficei,                        &
    !    &,                  cfkeep, cfliqkeep, cficekeep
                          rho, m0, tcgc, tcgci, corr, lfrcp, 0,               &
                      not_generic_size_dist,                                  &
                          piacw, one_over_tsi,                                &
                          l_use_agg_vt                                        &
    !    &,                  cf_transfer_diag, cfl_transfer_diag
    !    &,                  cff_transfer_diag
                         )
   ELSE
    CALL lsp_riming_sphere(points, timestep_mp,                               &
                           qcl, qcf_cry, t,                                   &
                           area_liq, area_mix, cfliq, cficei,                 &
    !    &,                  cfkeep, cfliqkeep, cficekeep
                           rho, m0, tcgc, tcgci, corr, lfrcp, 0,              &
                      not_generic_size_dist,                                  &
                           piacw, one_over_tsi,                               &
                           l_use_agg_vt                                       &
    !    &,                  cf_transfer_diag, cfl_transfer_diag
    !    &,                  cff_transfer_diag
                          )
   END IF  ! l_shape_rime
  END IF  ! l_crystals


  ! ======================================================================
  !             RIMING OF SNOW AGGREGATES BY CLOUD WATER (PSACW)
  ! ======================================================================

! Enhance riming using orographic water if required
  IF (l_orograin .AND. l_orogrime) THEN

!   Initialise all arrays used for riming enhancement CHANGE
    DO i = 1, points
      ql_orog(i) = zero
      cf_orog(i) = zero
      qsnow1a(i) = zero
      qsnow1b(i) = zero
      psacw1a(i) = zero
      psacw1b(i) = zero
      t1a(i) = zero
      t1b(i) = zero
      area_mix_orog(i) = zero
      area_liq_orog(i) = zero
      dqsnow(i) = zero
      altitude(i) = zero
      l_enh_rime(i) = .FALSE.
    END DO

!   Calculate orographic water MR and cloud fraction
    CALL lsp_orogwater( points, hmteff, zb,                                   &
      r_theta_levels_c,r_theta_surf_c, fv_cos_theta_latitude_c,               &
      p, t, q, qsl, esw, qcl, cfliq,                                          &
      ql_orog, cf_orog )

!   Loop through points to initialize temp variables for riming calls
    DO i = 1, points

!     Altitude above local surface for debug print statements
      altitude(i) = MAX( zero, r_theta_levels_c(i)-r_theta_surf_c(i) )

!     Save qrain and accr rate into temp variables for orog call
      qsnow1a(i) = qcf_agg(i)
      qsnow1b(i) = qcf_agg(i)
      t1a(i) = t(i)
      t1b(i) = t(i)
      psacw1a(i) = psacw(i)
      psacw1b(i) = psacw(i)
!     Save ql_orog for accretion enhancement test
      ql_orog1(i) = ql_orog(i)
!     Save into diagnostic for subgrid orographic water
      sfwater(i) = sfwater(i) + (ql_orog(i) / niters_mp)

!     MINIMUM overlap of orographic cloud with ICE (as for resolved)
      area_mix_orog(i) = MAX(zero, (cf_orog(i) + cfice(i) - one) )
      area_liq_orog(i) = MAX(zero, (cf_orog(i) - area_mix_orog(i)) )

    END DO ! Loop through all cloudy points


    IF (l_shape_rime) THEN

      CALL lsp_riming(points, timestep_mp,                                    &
                      qcl, qcf_agg, t,                                        &
                      area_liq, area_mix, cfliq, cficei,                      &
  !    &,                  cfkeep, cfliqkeep, cficekeep
                      rho, m0, tcg, tcgi, corr, lfrcp, 1,                     &
                      l_psd,                                                  &
                      psacw, one_over_tsi,                                    &
                      l_use_agg_vt                                            &
  !    &,                  cf_transfer_diag, cfl_transfer_diag
  !    &,                  cff_transfer_diag
                     )

      CALL lsp_riming(points, timestep_mp,                                    &
                      ql_orog, qsnow1b, t1b,                                  &
                      area_liq_orog, area_mix_orog, cf_orog, cficei,          &
  !    &,                  cfkeep, cfliqkeep, cficekeep
                      rho, m0, tcg, tcgi, corr, lfrcp, 1,                     &
                      l_psd,                                                  &
                      psacw1b, one_over_tsi,                                  &
                      l_use_agg_vt                                            &
  !    &,                  cf_transfer_diag, cfl_transfer_diag
  !    &,                  cff_transfer_diag
                     )

    ELSE

      CALL lsp_riming_sphere(points, timestep_mp,                             &
                       qcl, qcf_agg, t,                                       &
                       area_liq, area_mix, cfliq, cficei,                     &
  !    &,                  cfkeep, cfliqkeep, cficekeep
                       rho, m0, tcg, tcgi, corr, lfrcp, 1,                    &
                       l_psd,                                                 &
                       psacw, one_over_tsi,                                   &
                       l_use_agg_vt                                           &
  !    &,                   cf_transfer_diag, cfl_transfer_diag
  !    &,                   cff_transfer_diag
                     )


      CALL lsp_riming_sphere(points, timestep_mp,                             &
                       ql_orog, qsnow1b, t1b,                                 &
                       area_liq_orog, area_mix_orog, cf_orog, cficei,         &
  !    &,                  cfkeep, cfliqkeep, cficekeep
                       rho, m0, tcg, tcgi, corr, lfrcp, 1,                    &
                       l_psd,                                                 &
                       psacw1b, one_over_tsi,                                 &
                       l_use_agg_vt                                           &
  !    &,                   cf_transfer_diag, cfl_transfer_diag
  !    &,                   cff_transfer_diag
                     )

    END IF

!   Apply riming enhancement to actual model variables
    DO i = 1, points

!     Amount of orog water accreted (ensure positive)
      dqsnow(i) = MAX(zero, qsnow1b(i) - qsnow1a(i) )

      IF (dqsnow(i)>zero .AND. ql_orog1(i)>zero) THEN

        l_enh_rime(i) = .TRUE.
        qcf_agg(i) = qcf_agg(i) + dqsnow(i)
        q(i) = q(i) - dqsnow(i)

!       Add LH for cond+freezing of rimed orog water
        t(i) = t(i) + (dqsnow(i) * lsrcp)

!       Add mass transfer
        psacw(i) = psacw(i) +  (psacw1b(i) - psacw1a(i))

        ! Seeder feeder snow production rate / kg kg-1 s-1
        sfsnow(i) = sfsnow(i) + dqsnow(i) * one_over_tsi

      END IF

    END DO ! Loop through all cloudy points


  ELSE  ! seeder feeder off

    IF (l_shape_rime) THEN
      CALL lsp_riming(points, timestep_mp,                                    &
                   qcl, qcf_agg, t,                                           &
                   area_liq, area_mix, cfliq, cficei,                         &
  !    &,                  cfkeep, cfliqkeep, cficekeep
                      rho, m0, tcg, tcgi, corr, lfrcp, 1,                     &
                      l_psd,                                                  &
                      psacw, one_over_tsi,                                    &
                      l_use_agg_vt                                            &
  !    &,                  cf_transfer_diag, cfl_transfer_diag
  !    &,                  cff_transfer_diag
                     )
    ELSE
      CALL lsp_riming_sphere(points, timestep_mp,                             &
                          qcl, qcf_agg, t,                                    &
                          area_liq, area_mix, cfliq, cficei,                  &
  !    &,                  cfkeep, cfliqkeep, cficekeep
                       rho, m0, tcg, tcgi, corr, lfrcp, 1,                    &
                       l_psd,                                                 &
                       psacw, one_over_tsi,                                   &
                       l_use_agg_vt                                           &
  !    &,                   cf_transfer_diag, cfl_transfer_diag
  !    &,                   cff_transfer_diag
                     )

    END IF

  END IF



  IF (l_mcr_qgraup) THEN
    ! ======================================================================
    !             AUTOCONVERSION OF SNOW TO GRAUPEL  (PGAUT)
    ! ======================================================================

    CALL lsp_graup_autoc(points, timestep_mp,                                 &
                  qcf_agg, qgraup, t, rho,                                    &
                  psacw, psdep, pgaut, one_over_tsi                           &
                 )

    ! ======================================================================
    !              RIMING OF GRAUPEL BY CLOUD WATER (PGACW)
    ! ======================================================================
               ! Graupel is not included in ice cloud fraction
    CALL lsp_riming_sphere(points, timestep_mp,                               &
                           qcl, qgraup, t,                                    &
                           area_liq, area_mix, cfliq, cficei,                 &
    !    &,                    cfkeep, cfliqkeep, cficekeep
                           rho, m0, tcgg, tcggi, corr, lfrcp, 3,              &
                            not_generic_size_dist,                            &
                           pgacw, one_over_tsi,                               &
                           l_use_agg_vt                                       &
    !    &,                    cf_transfer_diag, cfl_transfer_diag
    !    &,                    cff_transfer_diag
                          )


    ! ======================================================================
    !              COLLECTION OF SNOW BY GRAUPEL (PGACS)
    ! ======================================================================
            ! For Graupel collecting snow the collision should result in a
            ! shattering of the ice and not a collection. For the moment,
            ! we shall not include this term when there is extra graupel
            ! production from snow-rain collisions. This should maintain
            ! bit-reproducibility with the present graupel scheme.
            !
            ! The break-up of ice aggregates to crystals will be dealt with
            ! at a later release.

    IF ( l_sr2graup ) THEN

      ! Set diagnostic to zero to show there is no transfer due to this process

      pgacs(:) = zero

    ELSE

      ! Perform collection process as the standard UM graupel scheme

      l_use_area=.FALSE.
      l_no_t_check=.FALSE.
      ice_type1=3  ! graupel
      ice_type2=1  ! aggregates

      CALL lsp_collection(points, timestep_mp,                                &
                  qgraup, qcf_agg, t,                                         &
                  area_mix, area_ice, cficei,                                 &
      !                          cf, cff,
                            rho, rhor, m0, tcgg, tcggi,                       &
                            tcg, tcgi, corr,                                  &
                            ice_type1, ice_type2,                             &
                            l_psd,                                            &
                            l_use_area, l_no_t_check,                         &
                            pgacs, one_over_tsi,                              &
                            l_use_agg_vt                                      &
      !                          cf_transfer_diag, cff_transfer_diag
                                  )

    END IF ! l_sr2graup

  END IF  ! l_mcr_graup

  ! ======================================================================
  !               COLLECTION OF RAIN BY ICE CRYSTALS (PIACR)
  ! ======================================================================
  IF (l_crystals) THEN
    CALL lsp_capture(points, timestep_mp,                                     &
                     qcf_cry, qrain, qgraup, t,                               &
                     area_liq, area_mix, area_ice, cficei,                    &
                     rainfrac,rain_liq,rain_mix,rain_ice,rain_clear,          &
    !    &,                   cfkeep, cficekeep
                           rho, rhor, m0, tcgc, tcgci, corr, dhir,            &
                           lfrcp, 0,                                          &
                           not_generic_size_dist,                             &
                           piacr, one_over_tsi,                               &
    !    &,                   cf_transfer_diag, cff_transfer_diag
                           rf_transfer_diag,                                  &
                           l_use_agg_vt                                       &
                          )
  END IF

  ! ======================================================================
  !               COLLECTION OF RAIN BY SNOW AGGREGATES (PSACR)
  ! ======================================================================
  CALL lsp_capture(points, timestep_mp,                                       &
                   qcf_agg, qrain, qgraup, t,                                 &
                   area_liq, area_mix, area_ice, cficei,                      &
                   rainfrac,rain_liq,rain_mix,rain_ice,rain_clear,            &
  !    &,                   cfkeep, cficekeep
                       rho, rhor, m0, tcg, tcgi, corr,                        &
                       dhir,lfrcp, 1,                                         &
                       l_psd,                                                 &
                       psacr, one_over_tsi,                                   &
  !    &,                   cf_transfer_diag, cff_transfer_diag
                       rf_transfer_diag,                                      &
                       l_use_agg_vt                                           &
                      )

  ! ======================================================================
  !                 EVAPORATION OF MELTING ICE CRYSTALS (PIMLTEVP)
  ! ======================================================================
  IF (l_crystals) THEN

    DO i = 1, points
         ! Calculate total ice content
      qcft(i) = qcf_cry(i) + qcf_agg(i)
    END DO

    CALL lsp_evap_snow(points, timestep_mp,                                   &
                       q, q_ice, qcf_cry, qcft, t, p, esw, qsl,               &
                       area_ice, cficei, cfkeep, cficekeep,                   &
                       rho, tcgc, tcgci,                                      &
                       corr2, rocor, lheat_correc_liq,                        &
                       lsrcp, 0,                                              &
                       not_generic_size_dist, pimltevp, one_over_tsi,         &
                       cf_transfer_diag, cff_transfer_diag,                   &
                       l_use_agg_vt                                           &
                      )
  END IF  ! l_crystals

  ! ======================================================================
  !                 EVAPORATION OF MELTING SNOW AGGREGATES (PSMLTEVP)
  ! ======================================================================
  DO i = 1, points
         ! Calculate total ice content
    qcft(i) = qcf_cry(i) + qcf_agg(i)
  END DO

  CALL lsp_evap_snow(points, timestep_mp,                                     &
                     q, q_ice, qcf_agg, qcft, t, p, esw, qsl,                 &
                     area_ice, cficei, cfkeep, cficekeep,                     &
                     rho, tcg, tcgi,                                          &
                     corr2, rocor, lheat_correc_liq, lsrcp, 1,                &
                     l_psd, psmltevp, one_over_tsi,                           &
                     cf_transfer_diag,cff_transfer_diag,                      &
                     l_use_agg_vt                                             &
                    )

       ! Perform the melting steps here

  ! ======================================================================
  !                    MELTING OF ICE CRYSTALS (PIMLT)
  ! ======================================================================
  IF (l_crystals) THEN

    DO i = 1, points
       ! Calculate total ice content
      qcft(i) = qcf_cry(i) + qcf_agg(i)
    END DO

    CALL lsp_melting(points, timestep_mp,                                   &
                     q, q_ice, qcf_cry, qcft, qrain, qsl, t, p,             &
                     area_liq, area_mix, area_ice, area_clear,              &
                     cfice, cficei, frac_ice_above,                         &
                     rainfrac, rain_liq, rain_mix,                          &
                     rain_ice, rain_clear, cfkeep, cficekeep,               &
                     rho, rhor, m0, tcg, tcgi, corr2, rocor,                &
                     lfrcp, 0,                                              &
                     not_generic_size_dist,                                 &
                     pimlt, one_over_tsi,                                   &
                     cf_transfer_diag, cff_transfer_diag,                   &
                     rf_transfer_diag,                                      &
                     l_use_agg_vt                                           &
                    )
  END IF  ! l_crystals

  ! ======================================================================
  !                    MELTING OF SNOW AGGREGATES (PSMLT)
  ! ======================================================================
  DO i = 1, points
       ! Calculate total ice content
    qcft(i) = qcf_cry(i) + qcf_agg(i)
  END DO

  CALL lsp_melting(points, timestep_mp,                                     &
                   q, q_ice, qcf_agg, qcft, qrain, qsl, t, p,               &
                   area_liq, area_mix, area_ice, area_clear,                &
                   cfice, cficei, frac_ice_above,                           &
                   rainfrac, rain_liq, rain_mix,                            &
                   rain_ice, rain_clear, cfkeep, cficekeep,                 &
                   rho, rhor, m0, tcgc, tcgci,                              &
                   corr2, rocor,                                            &
                   lfrcp, 1,                                                &
                   l_psd,                                                   &
                   psmlt, one_over_tsi,                                     &
                   cf_transfer_diag, cff_transfer_diag,                     &
                   rf_transfer_diag,                                        &
                   l_use_agg_vt                                             &
                  )

  ! ======================================================================
  !                    MELTING OF GRAUPEL (PGMLT)
  ! ======================================================================
  IF (l_mcr_qgraup) THEN
     ! Graupel does not update cloud fractions so there is no need
     ! to update qcft (it is not used)
    CALL lsp_melting(points, timestep_mp,                                   &
                     q, q_ice, qgraup, qcft, qrain, qsl, t, p,              &
                     area_liq, area_mix, area_ice, area_clear,              &
                     cfice, cficei, frac_ice_above,                         &
                     rainfrac, rain_liq, rain_mix,                          &
                     rain_ice, rain_clear, cfkeep, cficekeep,               &
                     rho, rhor, m0, tcgg, tcggi, corr2, rocor,              &
                     lfrcp, 3,                                              &
                     not_generic_size_dist,                                 &
                     pgmlt, one_over_tsi,                                   &
                     cf_transfer_diag, cff_transfer_diag,                   &
                     rf_transfer_diag,                                      &
                     l_use_agg_vt                                           &
                    )
  END IF  ! l_mcr_qgraup

  ! ======================================================================
  !                   EVAPORATION OF RAINDROPS (PREVP)
  ! ======================================================================
  CALL lsp_evap(points, timestep_mp, p, q, qrain, t,                          &
                q_ice, q_clear,                                               &
                area_liq, area_mix, area_ice, area_clear,                     &
                rainfrac, rain_liq, rain_mix,                                 &
                rain_ice, rain_clear,                                         &
                rho, corr, corr2, rocor,                                      &
                dhir, lcrcp, lheat_correc_liq,                                &
                qsl, esw,                                                     &
                prevp, rf_transfer_diag, one_over_tsi                         &
               )

  ! ======================================================================
  !             ACCRETION OF CLOUD DROPLETS ON RAINDROPS (PRACW)
  ! ======================================================================
  !  Enhanced by subgrid orographic motions
  ! ================================================================

  IF (l_orograin) THEN

!   Initialise all arrays used for accretion enhancement
    DO i = 1, points
      ql_orog(i) = zero
      cf_orog(i) = zero
      ql_orog1(i) = zero
      qrain1a(i) = zero
      qrain1b(i) = zero
      pracw1a(i) = zero
      pracw1b(i) = zero
      rain_liq_orog(i) = zero
      rain_mix_orog(i) = zero
      dqrain(i) = zero
    END DO

!   Calculate orographic water MR and cloud fraction
    CALL lsp_orogwater( points, hmteff, zb,                                   &
      r_theta_levels_c,r_theta_surf_c, fv_cos_theta_latitude_c,               &
      p, t, q, qsl, esw, qcl, cfliq,                                          &
      ql_orog, cf_orog )

!   Loop through points to initialize temp variables for accr calls
    DO i = 1, points

!     Save qrain into temp variables for orog call
      qrain1a(i) = qrain(i)
      qrain1b(i) = qrain(i)
      pracw1a(i) = pracw(i)
      pracw1b(i) = pracw(i)
!     Save ql_orog for outputting afterwards
      ql_orog1(i) = ql_orog(i)

!     Overlap with RAIN - assume max value (only use sum)
      rain_liq_orog(i) = MAX(zero, MIN( cf_orog(i), rainfrac(i) ) )
      rain_mix_orog(i) =  zero

    END DO

    CALL lsp_accretion(points, timestep_mp, qcl, qrain,                       &
                      cfliq, rainfrac, rain_liq, rain_mix,                    &
                      rho, corr, dhir,                                        &
                      pracw, one_over_tsi,                                    &
                      r_theta_levels_c, fv_cos_theta_latitude_c,              &
                      f_arr1, f_arr2, f_arr3                                  &
                          )

    CALL lsp_accretion(points, timestep_mp, ql_orog, qrain1b,                 &
                     cf_orog, rainfrac, rain_liq_orog, rain_mix_orog,         &
                     rho, corr, dhir,                                         &
                     pracw1b, one_over_tsi,                                   &
                     r_theta_levels_c, fv_cos_theta_latitude_c,               &
                     f_arr1, f_arr2, f_arr3                                   &
                          )


!   Apply accretion enhancement to actual model variables
!   ...(but only if riming wasnt enhanced at that point).
    DO i = 1, points

!     Subgrid orog water MR actually accreted by rain
      dqrain(i) = MAX(zero, qrain1b(i) - qrain1a(i) )

      IF (dqrain(i)>zero .AND. ql_orog1(i)>zero .AND.                         &
            .NOT. l_enh_rime(i) ) THEN

!       Add orographic accretion to that from resolved cloud
        qrain(i) = qrain(i) + dqrain(i)
        q(i) = q(i) - dqrain(i)
        t(i) = t(i) + (dqrain(i) * lcrcp)

!       Add mass transfer
        pracw(i) = pracw(i) +  (pracw1b(i) - pracw1a(i))

        ! Seeder feeder rain production rate / kg kg-1 s-1
        sfrain(i) = sfrain(i) + dqrain(i) * one_over_tsi

      END IF

    END DO

  ELSE  ! seeder feeder off

    CALL lsp_accretion(points, timestep_mp, qcl, qrain,                       &
                       cfliq, rainfrac, rain_liq, rain_mix,                   &
                       rho, corr, dhir,                                       &
                       pracw, one_over_tsi,                                   &
                      r_theta_levels_c, fv_cos_theta_latitude_c,              &
                      f_arr1, f_arr2, f_arr3                                  &
                          )

  END IF


  ! ======================================================================
  !              AUTOCONVERSION OF CLOUD LIQUID TO RAIN (PRAUT)
  ! ======================================================================
  CALL lsp_autoc(points, timestep_mp, qcl, qrain, t, p,                       &
                cfliq, rhcpt,                                                 &
  !    &,                  cfkeep, cfliqkeep
                    area_liq, area_mix, area_ice, rainfrac,                   &
                    rain_liq, rain_mix, rain_ice, rain_clear,                 &
                    rho, rhor, corr2, lcrcp,                                  &
                    one_over_tsi, praut, rf_transfer_diag,                    &
  !    &,                  cf_transfer_diag, cfl_transfer_diag
                    land_fract,                                               &
                    n_drop_tpr, n_drop_out,                                   &
                    r_theta_levels_c, fv_cos_theta_latitude_c,                &
                    f_arr1, f_arr2, f_arr3                                    &
                   )

! ======================================================================


!                  COPY ICE/SNOW VARIABLES AND FLUXES
!                        TO OUTPUT VARIABLES


! ======================================================================
! Copy contents of ice/snow variables and fluxes to output variables
! to fall into next layer down
! ----------------------------------------------------------------------

DO i = 1, points

  IF (l_mcr_qcf2) THEN ! two ice prognostics
    qcf(i)  = qcf_agg(i)
    qcf2(i) = qcf_cry(i)
    snow_cry(i) = snowt_cry(i)
    snow_agg(i) = snowt_agg(i)
  ELSE ! only one ice prognostic, put all snow in to snow_agg
    qcf(i) = qcf_cry(i) + qcf_agg(i)
    snow_cry(i) = zero  ! Redundant variable
    snow_agg(i) = snowt_cry(i) + snowt_agg(i)
  END IF ! on l_mcr_qcf2

  IF (l_mcr_qrain)  rainrate(i) = rainratet(i)
  IF (l_mcr_qgraup) grauprate(i) = graupratet(i)

END DO ! Points

! ======================================================================
!              NUMERICAL TIDYING UP OF SMALL VALUES
! ======================================================================
CALL lsp_tidy(points, one_over_tsi,                                           &
             q, qcl, qcf, qcf2, qrain, t,                                     &
             area_liq, area_mix, area_ice,                                    &
             cfice, cficei, cfkeep, cfliqkeep, cficekeep,                     &
             rainfrac, rain_liq, rain_mix, rain_ice, rain_clear,              &
             q_ice, qs, qsl, snow_agg, snow_cry,                              &
             rho, rhor, p,                                                    &
             cttemp,dhi,dhir,frac_ice_above,                                  &
             lcrcp, lfrcp, lsrcp,                                             &
             psdep, pidep, psmlt, pimlt, prevp,                               &
             cf_transfer_diag,cfl_transfer_diag,                              &
             cff_transfer_diag, rf_transfer_diag                              &
             )

DO i = 1, points
      !------------------------------------------------
      ! Now update fraction of ice in layer above
      ! for next layer down
      !------------------------------------------------
  frac_ice_above(i)=cficekeep(i)

END DO ! Points


!========================================================================
!                   RADAR REFLECTIVITY DIAGNOSTICS
!========================================================================

IF (l_ref_diag) THEN
  CALL mphys_reflec(points, rho,t, qgraup, qcf, qcf_cry, qrain, qcl,          &
                    n_drop_tpr, cficekeep, cfliqkeep, rainfrac_impr, tcg,     &
                    tcgc, dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_l, dbz_r )

END IF

! ======================================================================
!       IF DIAGNOSTIC RAIN, CONVERT MASS (kg/kg) TO FLUX (kg/m2/s)
! ======================================================================
IF (.NOT. l_mcr_qrain) THEN

  IF (l_mr_physics) THEN

    DO i = 1, points

      ! Use mixing ratio formulation
      rainrate(i) = qrain(i) * rhodz_dry(i) / timestep

    END DO ! points

  ELSE ! l_mr_physics

    DO i = 1, points

      ! Use specific humidity formulation
      rainrate(i) = qrain(i) * rhodz_moist(i) / timestep

    END DO ! points

  END IF  ! l_mr_physics

END IF ! on prognostic rain mixing ratio

! ----------------------------------------------------------------------
!   End of the LSP_ICE subroutine
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE lsp_ice
END MODULE lsp_ice_mod
