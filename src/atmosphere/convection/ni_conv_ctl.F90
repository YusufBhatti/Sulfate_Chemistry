! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine NI_conv_ctl

MODULE ni_conv_ctl_mod

IMPLICIT NONE

!---------------------------------------------------------------------------
! Description: Interface to Atmospheric Physics convection code.

! method:
! Note all working arrays other than prognostics are defined without surface
! level values for the ENDGame case. This matches what currently happens for
! new dynamics, i.e. the convection scheme only operates on model levels
! above the surface. All loops over model levels should always work from 1
! to model levels not 0.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 8.3.
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='NI_CONV_CTL_MOD'

CONTAINS

SUBROUTINE ni_conv_ctl (                                          &

! Parallel variables
  at_extremity, n_proc                                            &
, delta_lambda,delta_phi                                          &

! Parameters for iterative SISL
, numcycles, cycleno                                              &

! Model dimensions.
, row_length, rows                                                &
, rowsrowlength                                                   &
, bl_levels, n_cca_lev                                            &
, tr_vars, tr_ukca                                                &

! Model switches
, l_calc_dxek, l_q_interact                                       &
, l_spec_z0                                                       &

! in coordinate information
, z_rho, z_theta                                                  &

! SCM/Idealised UM
, flux_e,flux_h,ustar_in,z0m_scm,z0h_scm                          &
! require for conv_diag call
, t_surf, zh, u_0_p, v_0_p, fb_surf, ustar                        &
, ls_rain, ls_snow                                                &

! SCM diagnostics (dummy in full UM)
, nscmdpkgs, l_scmdiags, conv_mode                                &

! in data fields.
, rho,  rho_only, rho_theta, rho_dry, rho_dry_theta               &
, u_p, v_p, ustar_p, vstar_p, w, p, p_star, exner_rho_levels      &
, land_sea_mask, flandg, ice_fract                                &
, tstar_land, tstar_sea, tstar_sice, z0msea                       &
, p_layer_boundaries, p_layer_centres                             &
, exner_layer_boundaries, exner_layer_centres                     &
, t1_sd, q1_sd, exner_theta_levels                                &
, uw0_p, vw0_p, w_max, zlcl, zlcl_uv, ztop, dzh, qcl_inv_top      &
, entrain_coef, conv_type                                         &
, cumulus, l_shallow, l_congestus, l_congestus2, it_mid_level     &
, l_pc2_diag_sh_pts                                               &
, no_cumulus                                                      &
, ntml, ntpar                                                     &
, wstar, wthvs, delthvu, ql_ad, qsat_lcl, ftl, fqt                &
, shallowc, cu_over_orog, cape_undilute, cin_undilute             &
, deep_flag, past_precip, past_conv_ht                            &

! in Prognostics at beginning of time step (winds passed in higher up)
, theta_n, q_n, qcl_n, qcf_n, qrain_n, qgraup_n, qcf2_n           &
, cf_liquid_n, cf_frozen_n, bulk_cf_n                             &
! in Cold-pool fields
, g_ccp, h_ccp                                                    &
! in/out primary fields with all increments so far added on
, theta_star, q_star, qcl_star, qcf_star                          &
, qcf2_star, qrain_star, qgraup_star                              &
, cf_liquid_star, cf_frozen_star, bulk_cf_star                    &
! out convective momentum transport tendencies on p-grid
, dubydt_pout, dvbydt_pout                                        &
! in/out values with all increments added
, conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip         &
!  In  - total wind increments before convection on p grid
, r_u_p, r_v_p                                                    &
! in/out
, aerosol                                                         &
, dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6     &
, so2, so4_aitken, so4_accu, so4_diss                             &
, dms, nh3, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd    &
, bmass_cld, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss    &
, co2,  free_tracers, tracer_ukca                                 &
, ozone_tracer                                                    &

! out fields
, cca0_dp, cca0_md, cca0_sh                                       &
, cca0, ccw0, ccb0, cct0, cclwp0, lcbase0, cca0_2d, lctop, lcca   &
, cca,  ccw,  ccb,  cct,  cclwp,  lcbase,  cca_2d                 &
, conv_rain, conv_snow,  ddmfx                                    &

! error information
, error_code  )


!$ USE omp_lib

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                          &
   wdims, tdims, pdims, tdims_s, pdims_s,                                 &
   tdims_l, ScmRowLen, ScmRow

USE nlsizes_namelist_mod,   ONLY: model_levels

! Model time stepping information
USE timestep_mod, ONLY: timestep, timestep_number, recip_timestep

USE segments_mod, ONLY:                                                   &
  segment_type, meta_segment_type,                                        &
  segments_mod_seg_meta, segments_mod_segments

USE autotune_mod, ONLY:       &
    autotune_type,            &
    autotune_init,            &
    autotune_get_trial_size,  &
    autotune_start_region,    &
    autotune_stop_region,     &
    autotune_advance,         &
    autotune_report

! Model level heights from centre of Earth
USE level_heights_mod, ONLY: &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

USE gen_phys_inputs_mod, ONLY: l_mr_physics

USE cv_run_mod,  ONLY:                                                   &
    i_convection_vn,                                                     &
    i_convection_vn_5a,                                                  &
    i_convection_vn_6a,                                                  &
    l_fix_udfactor,          l_mom,                                      &
    l_safe_conv,             l_rediagnosis,         ud_factor,           &
    iconv_shallow,           iconv_congestus,       iconv_deep,          &
                             n_conv_calls,                               &
    rad_cloud_decay_opt,     cld_life_opt,          cca_min,             &
    fixed_cld_life,                                 l_murk_conv,         &
    l_dcpl_cld4pc2,          qmin_conv,             l_ccrad,             &
    l_3d_cca,                l_pc2_diag_sh,         l_conv_hist,         &
    l_conv_prog_group_1,     l_conv_prog_group_2,   l_conv_prog_group_3, &
    l_conv_prog_precip,      tau_conv_prog_precip

USE tuning_segments_mod, ONLY: l_autotune_segments,                      &
                               a_convect_segments,a_convect_seg_size

USE model_domain_mod, ONLY: model_type, mt_global, mt_single_column
USE dynamics_input_mod,  ONLY:                                            &
    l_check_moist_inc

! Extra moist variable switches
USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

USE cv_hist_constants_mod, ONLY:                                  &
    decay_period

USE cv_param_mod,  ONLY:                                          &
    cld_life_constant,       cld_life_func_hgt,     rad_decay_off,&
    rad_decay_full_timestep, rad_decay_conv_substep,              &
    dthetadt_conv_active_threshold, conv_prog_precip_min_threshold

! Convective diagnostic output arrays
USE cv_diagnostic_array_mod, ONLY:                                      &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,dwn_flux    &
       ,ep1, ep2, ep3

! Flags controlling convective diagnostic output
USE cv_stash_flg_mod, ONLY: l_apply_diag, flg_w_eqn 

USE s_scmop_mod, ONLY: scmdiag_convss

! Tracer switches
USE run_aerosol_mod, ONLY: l_sulpc_so2, l_sulpc_dms, l_sulpc_nh3,       &
                           l_soot, l_ocff, l_biomass, l_nitrate
USE carbon_options_mod, ONLY: l_co2_interactive
USE dust_parameters_mod, ONLY: l_dust


USE turb_diff_ctl_mod, ONLY: delta_smag

! Flags controlling the saving of fields
USE jules_surface_mod, ONLY: ISrfExCnvGust, IP_SrfExWithCnv

USE cloud_inputs_mod, ONLY: i_pc2_conv_coupling,                  &
                            l_fixbug_pc2_mixph

USE rad_input_mod, ONLY: l_cca_dp_prog, l_cca_md_prog, l_cca_sh_prog, &
                         l_use_cariolle

! Add convection tendencies to physics_tendencies module SPT scheme
USE physics_tendencies_mod,  ONLY:                                     &
    l_retain_conv_tendencies, dt_conv, dq_conv

USE UM_ParCore, ONLY: mype
USE UM_ParParams, ONLY: pnorth, psouth

! Subroutines
USE conv_diag_5a_mod, ONLY: conv_diag_5a
USE conv_diag_6a_mod, ONLY: conv_diag_6a
USE glue_conv_5a_mod, ONLY: glue_conv_5a
USE glue_conv_6a_mod, ONLY: glue_conv_6a
USE tracer_total_var_mod, ONLY: tracer_total_var
USE tracer_copy_mod,      ONLY: tracer_copy
USE tracer_restore_mod,   ONLY: tracer_restore
USE all_scav_calls_mod,   ONLY: all_scav_calls
USE pc2_from_conv_ctl_mod, ONLY: pc2_from_conv_ctl
USE scm_convss_dg_mod, ONLY: scm_convss_dg_type,                       &
                             scm_convss_dg_allocate,                   &
                             scm_convss_dg_deallocate,                 &
                             scm_convss_dg_initzero,                   &
                             scm_convss_dg_output
USE update_conv_diags_mod, ONLY: update_conv_diags
USE update_conv_cloud_mod, ONLY: update_conv_cloud
USE check_dmoist_inc_mod,  ONLY: check_dmoist_inc
USE save_conv_diags_mod,   ONLY: save_conv_diags

! Required by most routines - message printing and error reporting
USE umPrintMgr
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE conv_pc2_init_mod, ONLY: conv_pc2_init
IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments
! arguments with intent in. ie: input variables.
! Parallel setup variables
INTEGER, INTENT(IN)  ::      &
  n_proc                     & ! Total number of processors
, numcycles                  & ! Number of cycles
, cycleno                      ! cycle number

REAL , INTENT(IN)  :: delta_lambda,delta_phi ! model grid spacing in radians

LOGICAL, INTENT(IN) ::  &
  at_extremity(4)         ! Indicates if this processor is at north, south
                          ! east or west of the processor grid

! Model dimensions
INTEGER, INTENT(IN)  ::      &
  row_length                 & ! Row length
, rows                       & ! Number of rows
, rowsrowlength              & ! rows*row_length
, bl_levels                  & ! Number of boundary layer levels
, n_cca_lev                  & ! Number of levels for conv cloud
                               ! amount: 1 for 2D, nlevs for 3D.
, tr_vars                    & ! number of free tracers
, tr_ukca                      ! number of ukca tracers

! Model switches
LOGICAL, INTENT(IN)  ::  &
  l_calc_dxek            & ! Switch for calculation of condensate increment
, l_q_interact           & ! Switch allows overwriting of parcel variables
                           ! when calculating condensate increments.
, l_spec_z0                ! true if roughness length has been specified

! Co-ordinate arrays  ! local vertical co-ordinate information
REAL, INTENT(IN) ::                               &
  z_rho(pdims%i_end,pdims%j_end,pdims%k_end)      & ! rho levels
, z_theta(tdims%i_end,tdims%j_end,tdims%k_end)      ! theta levels

! Used in SCM for prescribed surface flux forcing (required for call to
!  conv_diag)
REAL, INTENT(IN) ::           &
  flux_e(row_length, rows)    & ! Surface latent heat flux (W/m^2)
 ,flux_h(row_length, rows)    & ! Surface sensible heat flux (W/m^2)
 ,ustar_in(row_length, rows)  & ! Specified surface friction velocity (m/s)
 ,z0m_scm(row_length, rows)   & ! Fixed sea surface roughness lengths
 ,z0h_scm(row_length, rows)     ! for momentum and heat, used in SCM
! Required for call to conv_diag
REAL, INTENT(INOUT) ::        &
  t_surf(row_length, rows)      ! surface temperature (K)

REAL, INTENT(IN) ::                                                      &
   fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   ustar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! GBM surface friction velocity

REAL, INTENT(IN) ::           &
  zh(row_length, rows)        & ! boundary layer height (m)
, u_0_p(row_length, rows)     & ! set to zero (surface current m/s)
, v_0_p(row_length, rows)       ! set to zero (surface current m/s)

REAL, INTENT(INOUT) ::       &
  ls_rain(row_length, rows)  &  ! Large-scale rainfall from section 4
, ls_snow(row_length, rows)     ! Large-scale snowfall from section 4

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN)  ::  &
  nscmdpkgs                ! No of SCM diagnostics packages

LOGICAL, INTENT(IN)  ::  &
  l_scmdiags(nscmdpkgs)    ! Logicals for SCM diagnostics packages

INTEGER, INTENT(IN) :: conv_mode ! SCM switch for whether convection running
                                 ! in fully-coupled mode or diagnostics-only

! Data fields
REAL, INTENT(IN) ::                                    &
  rho(pdims_s%i_start:pdims_s%i_end,                   & ! density *r*r (kg/m)
      pdims_s%j_start:pdims_s%j_end,                   & !
      pdims_s%k_start:pdims_s%k_end)                   & !
, rho_only(pdims%i_end,pdims%j_end,pdims%k_end)        & ! wet density (kg/m3)
, rho_theta(tdims%i_end,tdims%j_end,tdims%k_end-1)     & ! wet on th lev (kg/m3)
, rho_dry(pdims%i_end,pdims%j_end,pdims%k_end)         & ! dry density (kg/m3)
, rho_dry_theta(tdims%i_end,tdims%j_end,tdims%k_end-1) & ! dry on th lev (kg/m3)
, u_p(pdims%i_end,pdims%j_end,pdims%k_end)             & ! U (m/s) P-grid
, v_p(pdims%i_end,pdims%j_end,pdims%k_end)             & ! V (m/s) P-grid
, ustar_p(pdims%i_end,pdims%j_end,pdims%k_end)         & ! U+du (m/s) P-grid
, vstar_p(pdims%i_end,pdims%j_end,pdims%k_end)         & ! V+dv (m/s) P-grid
, w(wdims%i_start:wdims%i_end,        & ! W (m/s) (copy of prognostic array)
    wdims%j_start:wdims%j_end,        &
                0:wdims%k_end)        & ! Not exact match to module values

, p(pdims_s%i_start:pdims_s%i_end,    & ! pressure (Pa)
    pdims_s%j_start:pdims_s%j_end,    &
    pdims_s%k_start:pdims_s%k_end)    &
, p_star(row_length, rows)                         & ! surface pressure
, exner_rho_levels(pdims_s%i_start:pdims_s%i_end,  & ! Exner on rho level
                   pdims_s%j_start:pdims_s%j_end,  & !
                   pdims_s%k_start:pdims_s%k_end)  &
, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)

LOGICAL, INTENT(IN) ::             &
  land_sea_mask(row_length, rows)    ! land sea land

REAL, INTENT(IN) ::                &
  flandg(pdims_s%i_start:pdims_s%i_end, & ! Land fraction of gridbox
       pdims_s%j_start:pdims_s%j_end)   & ! on all points
 ,ice_fract(row_length,rows)              ! fraction of sea that has ice

REAL, INTENT(IN) :: &
  tstar_land(row_length, rows) &
 ,tstar_sea(row_length, rows)  & ! Surface T on sea
 ,tstar_sice(row_length, rows) & ! Surface T on sea-ice
 ,z0msea(row_length, rows)       ! Sea surface roughness

REAL, INTENT(IN) ::                                                    &
  p_layer_boundaries(pdims%i_end,pdims%j_end,0:pdims%k_end)            &
      ! pressure at layer boundaries. Same as p except at
      ! bottom level = pstar, and at top = 0.
, p_layer_centres(tdims%i_end,tdims%j_end,0:tdims%k_end)               &
      ! pressure at layer centres. Same as p_theta_levels
      !except bottom level = pstar, and at top = 0.
, exner_layer_boundaries(pdims%i_end,pdims%j_end,0:pdims%k_end)        &
, exner_layer_centres(tdims%i_end,tdims%j_end,0:tdims%k_end)

REAL,INTENT(IN) ::           &
  t1_sd(row_length, rows)    & ! set to zero initially
, q1_sd(row_length, rows)      ! set to zero initially

REAL,INTENT(IN) ::               &
  uw0_p(row_length, rows)        & ! X-COMP SURFACE STRESS ON P GRID
, vw0_p(row_length, rows)        & ! Y-COMP SURFACE STRESS ON P GRID
, w_max(row_length, rows)          ! max w in column

REAL, INTENT (IN) ::         g_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)          &
!            gridbox c.c.p. reduced gravity (m/s^2)
                        ,    h_ccp( pdims%i_start:pdims%i_end,          &
                                    pdims%j_start:pdims%j_end)
!            gridbox c.c.p. height (m)

! Data arrays coming in from conv_diag & BL

REAL,INTENT(INOUT) ::            &
  zlcl(row_length, rows)         & ! accurate lifting condensation level(m)
, zlcl_uv(row_length, rows)      & ! LCL at level for uv calculation (m)
, ztop(row_length, rows)         & ! Top of cloud layer (m).
, dzh(row_length, rows)          & ! Inversion thickness (m).
, qcl_inv_top(row_length,rows)   & ! Parcel water content at inversion top
, entrain_coef(row_length, rows)   ! Entrainment coefficient

! Convective type array ::
INTEGER, INTENT(INOUT) :: conv_type(row_length, rows)
              ! Integer index describing convective type:
              !    0=no convection
              !    1=non-precipitating shallow
              !    2=drizzling shallow
              !    3=warm congestus
              !    ...
              !    8=deep convection

LOGICAL, INTENT(INOUT) ::          &
  cumulus(row_length, rows)        & ! Logical switch from boundary
                                     !    layer for presence of Cu
, l_shallow(row_length, rows)      & ! Logical switch for shallow Cu
, l_congestus(row_length, rows)    & ! Logical switch for congestus
, l_congestus2(row_length, rows)   & ! congestus in descending air
, it_mid_level(row_length, rows)   & ! mid-level convection occurred
, l_pc2_diag_sh_pts(row_length,rows) ! Carry diagnostic shallow convective
                                     ! information for PC2
LOGICAL, INTENT(IN) ::             &
  no_cumulus(row_length, rows)       ! Points which BL says cannot be
                                     ! cumulus

INTEGER, INTENT(INOUT) ::   &
  ntml(row_length, rows)    & ! Top level of surface mixed layer
, ntpar(row_length, rows)     ! Top level of initial parcel ascent

REAL,INTENT(INOUT) ::            &
  wstar(row_length, rows)        & ! Mixed layer convective velocity scale
, wthvs(row_length, rows)        & ! Surface flux of THV
, delthvu(row_length, rows)      & ! Integral of undilute parcel
                                   ! buoyancy over convective cloud layer
, ql_ad(row_length, rows)        & ! adiabatic liquid water content
                                   ! at inversion (kg/kg)
, qsat_lcl(row_length, rows)       ! qsat at LCL (kg/kg)

REAL,INTENT(IN) ::               &
  ftl(row_length, rows)          & ! Surface sensible heat flux from BL
                                   ! (W/m2) i.e. cp*rho*w'tl'
, fqt(row_length, rows)          & ! Total water flux from surface
                                   ! (kg/m2/s) i.e. rho*w'qT'
, shallowc(row_length,rows)      & ! Indicator set to 1.0 if shallow,
                                   !  0.0 if not shallow or not cumulus
, cu_over_orog(row_length,rows)    ! Indicator for cumulus over steep
                                   ! orography. Indicator set to 1.0 if
                                   ! true, 0.0 if false. Exclusive.

REAL,INTENT(INOUT) ::            &
  cape_undilute(row_length,rows) & ! Undilute CAPE from parcel ascent m2/s2
, cin_undilute(row_length,rows)    ! Undilute CIN from parcel ascent m2/s2

! Past history of convection - only used if prognostics requested otherwise
! these fields should remain unset as no space is allocated at the top
! level of the UM.

REAL,INTENT(INOUT) ::          &
  deep_flag(row_length,rows)   &  ! Value between 0.0 and 1.0
                                  ! 1 if deep convection last timestep
 ,past_precip(row_length,rows) &  ! Rate of convective precip (kg/m2/s)
                                  ! decayed over time.
 ,past_conv_ht(row_length,rows)   ! Height of previous convection (m)
                                

! Prognostics at start of time step
REAL,INTENT(IN) ::                           &
  theta_n(tdims_s%i_start:tdims_s%i_end,     & ! theta (K)
          tdims_s%j_start:tdims_s%j_end,     &
          tdims_s%k_start:tdims_s%k_end)     &
, q_n(tdims_l%i_start:tdims_l%i_end,         & ! q (kg/kg)
      tdims_l%j_start:tdims_l%j_end,         &
      tdims_l%k_start:tdims_l%k_end)         &
, qcl_n(tdims_l%i_start:tdims_l%i_end,       & ! qcl (kg/kg)
        tdims_l%j_start:tdims_l%j_end,       &
        tdims_l%k_start:tdims_l%k_end)       &
, qcf_n(tdims_l%i_start:tdims_l%i_end,       & ! qcf (kg/kg)
        tdims_l%j_start:tdims_l%j_end,       &
        tdims_l%k_start:tdims_l%k_end)       &
, qrain_n(tdims_l%i_start:tdims_l%i_end,     & ! qrain (kg/kg)
          tdims_l%j_start:tdims_l%j_end,     &
          tdims_l%k_start:tdims_l%k_end)     &
, qgraup_n(tdims_l%i_start:tdims_l%i_end,    & ! qraup (kg/kg)
           tdims_l%j_start:tdims_l%j_end,    &
           tdims_l%k_start:tdims_l%k_end)    &
, qcf2_n(tdims_l%i_start:tdims_l%i_end,      & ! qcf2 (kg/kg)
         tdims_l%j_start:tdims_l%j_end,      &
         tdims_l%k_start:tdims_l%k_end)      &
, cf_liquid_n(tdims_l%i_start:tdims_l%i_end, & ! liquid cloud fraction
              tdims_l%j_start:tdims_l%j_end, &
              tdims_l%k_start:tdims_l%k_end) &
, cf_frozen_n(tdims_l%i_start:tdims_l%i_end, & ! ice cloud fraction
              tdims_l%j_start:tdims_l%j_end, &
              tdims_l%k_start:tdims_l%k_end) &
, bulk_cf_n(tdims_l%i_start:tdims_l%i_end,   & ! total cloud fraction
            tdims_l%j_start:tdims_l%j_end,   &
            tdims_l%k_start:tdims_l%k_end)

! '_star' = '_n' + all increments up to now in time step
! _star variables held on same levels as prognostics in the vertical for
! ENDGame. Note Convection will not update any variables stored on level 0.
REAL,INTENT(INOUT) ::                                                    &
  theta_star(tdims%i_start:tdims%i_end,                                  &
             tdims%j_start:tdims%j_end,                                  &
             tdims%k_start:tdims%k_end)                                  &
, q_star(tdims%i_start:tdims%i_end,                                      &
         tdims%j_start:tdims%j_end,                                      &
         tdims%k_start:tdims%k_end)                                      &
, qcl_star(tdims%i_start:tdims%i_end,                                    &
           tdims%j_start:tdims%j_end,                                    &
           tdims%k_start:tdims%k_end)                                    &
, qcf_star(tdims%i_start:tdims%i_end,                                    &
           tdims%j_start:tdims%j_end,                                    &
           tdims%k_start:tdims%k_end)
REAL,INTENT(IN) ::                                                       &
  qcf2_star(tdims%i_start:tdims%i_end,                                   &
            tdims%j_start:tdims%j_end,                                   &
            tdims%k_start:tdims%k_end)                                   &
, qrain_star(tdims%i_start:tdims%i_end,                                  &
             tdims%j_start:tdims%j_end,                                  &
             tdims%k_start:tdims%k_end)                                  &
, qgraup_star(tdims%i_start:tdims%i_end,                                 &
              tdims%j_start:tdims%j_end,                                 &
              tdims%k_start:tdims%k_end)
REAL,INTENT(INOUT) ::                                                    &
  cf_liquid_star(tdims%i_start:tdims%i_end,                              &
                 tdims%j_start:tdims%j_end,                              &
                 tdims%k_start:tdims%k_end)                              &
, cf_frozen_star(tdims%i_start:tdims%i_end,                              &
                 tdims%j_start:tdims%j_end,                              &
                 tdims%k_start:tdims%k_end)                              &
, bulk_cf_star(tdims%i_start:tdims%i_end,                                &
               tdims%j_start:tdims%j_end,                                &
               tdims%k_start:tdims%k_end)
! OUT convective momentum transport tendencies on p-grid
REAL, INTENT(OUT) :: dubydt_pout ( pdims_s%i_start:pdims_s%i_end,        &
                                   pdims_s%j_start:pdims_s%j_end,        &
                                   pdims_s%k_start:pdims_s%k_end )
REAL, INTENT(OUT) :: dvbydt_pout ( pdims_s%i_start:pdims_s%i_end,        &
                                   pdims_s%j_start:pdims_s%j_end,        &
                                   pdims_s%k_start:pdims_s%k_end )

! arguments with intent in/out. ie: input variables changed on output.
! Convection prognostics
REAL, INTENT (INOUT) :: conv_prog_1(      tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_2(      tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )
REAL, INTENT (INOUT) :: conv_prog_3(      tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )
! Surface precipitation based 3d convective prognostic in kg/m2/s
REAL, INTENT (INOUT) :: conv_prog_precip( tdims_s%i_start:tdims_s%i_end,  &
                                          tdims_s%j_start:tdims_s%j_end,  &
                                          tdims_s%k_start:tdims_s%k_end )

REAL,INTENT(IN)  ::                  &
r_u_p(pdims%i_start:pdims%i_end,     &  ! r_u (m/s)on P-grid
      pdims%j_start:pdims%j_end,     &
      pdims%k_start:pdims%k_end)     &
,r_v_p(pdims%i_start:pdims%i_end,    &  ! r_v (m/s) on P-grid
      pdims%j_start:pdims%j_end,     &
      pdims%k_start:pdims%k_end)

! Tracers
REAL,INTENT(INOUT) ::                         &
  aerosol(tdims_s%i_start:tdims_s%i_end,      & ! aerosol tracer
          tdims_s%j_start:tdims_s%j_end,      &
          tdims_s%k_start:tdims_s%k_end)      &
, dust_div1(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div1
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div2(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div2
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div3(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div3
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div4(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div4
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div5(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div5
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div6(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div6
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)
                                                !DIV1,DIV2,...,DIV6 DUST
REAL,INTENT(INOUT) ::                         &
  so2(tdims_s%i_start:tdims_s%i_end,          & ! SO2
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, so4_aitken(tdims_s%i_start:tdims_s%i_end,   & ! SO4 aitken mode
             tdims_s%j_start:tdims_s%j_end,   &
             tdims_s%k_start:tdims_s%k_end)   &
, so4_accu(  tdims_s%i_start:tdims_s%i_end,   & ! SO4 accumulation mode
             tdims_s%j_start:tdims_s%j_end,   &
             tdims_s%k_start:tdims_s%k_end)   &
, so4_diss(tdims_s%i_start:tdims_s%i_end,     & ! SO4 dissipation mode
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, dms(tdims_s%i_start:tdims_s%i_end,          & ! DMS
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, nh3(tdims_s%i_start:tdims_s%i_end,          & ! NH3
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, soot_new(tdims_s%i_start:tdims_s%i_end,     & ! New soot
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, soot_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged soot
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, soot_cld(tdims_s%i_start:tdims_s%i_end,     & ! soot in cloud
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, bmass_new(tdims_s%i_start:tdims_s%i_end,    & ! New biomass
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, bmass_agd(tdims_s%i_start:tdims_s%i_end,    & ! Aged biomass
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, bmass_cld(tdims_s%i_start:tdims_s%i_end,    & ! Biomass in cloud
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                         &
  ocff_new(tdims_s%i_start:tdims_s%i_end,     & ! New ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, ocff_agd(tdims_s%i_start:tdims_s%i_end,     & ! Aged ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, ocff_cld(tdims_s%i_start:tdims_s%i_end,     & ! Cloud ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_acc(tdims_s%i_start:tdims_s%i_end,     & ! nitrate accumulation mode
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_diss(tdims_s%i_start:tdims_s%i_end,    & ! nitrate dissipation mode
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, co2(tdims_s%i_start:tdims_s%i_end,          & ! CO2
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                                                  &
  free_tracers(tdims_s%i_start:tdims_s%i_end,                          &
               tdims_s%j_start:tdims_s%j_end,                          &
               tdims_s%k_start:tdims_s%k_end, tr_vars)                 &
, tracer_ukca(tdims_s%i_start:tdims_s%i_end,                           &
              tdims_s%j_start:tdims_s%j_end,                           &
              tdims_s%k_start:tdims_s%k_end, tr_ukca)                  &
, ozone_tracer( tdims_s%i_start:tdims_s%i_end,          & ! Ozone tracer
                tdims_s%j_start:tdims_s%j_end,          &
                tdims_s%k_start:tdims_s%k_end)

! arguments with intent out. ie: output variables.
REAL, INTENT(OUT) ::                 &
  conv_rain(row_length, rows)        & ! Convective surface rainfall (kg/m2/s)
, conv_snow(row_length, rows)          ! Convective surface snowfall (kg/m2/s)

REAL, INTENT(OUT) :: ddmfx(row_length, rows)
!           ! Convective downdraught mass-flux at cloud-base

! Section 0 cloud properties:
!   These are passed to radiation and are subject to PC2/CCRad
!   and/or Cloud Property decay.

REAL, INTENT(INOUT) ::                  &
  cca0    (row_length,rows,n_cca_lev)   &! Cnv.Cld Amount (0-1)
, cca0_dp (row_length,rows,n_cca_lev)   &! Deep cnv.Cld Amount (0-1)
, cca0_md (row_length,rows,n_cca_lev)   &! Mid-level cnv.Cld Amount (0-1)
, cca0_sh (row_length,rows,n_cca_lev)   &! Shallow cnv.Cld Amount (0-1)
, ccw0    (tdims%i_end,tdims%j_end,tdims%k_end)  &! Cnv.Cld Water (kg/kg)
, cclwp0  (row_length,rows)              ! Cloud Condensed water path
                                         ! (kg/m^2)

REAL, INTENT(INOUT) ::                  &! Cnv.Cld Amount (2d) with
  cca0_2d (row_length,rows)              ! no anvil
          ! Only used if: l_3d_cca = .False. .and. i_cld_pc2 .and.
          ! l_pc2_diag_sh

INTEGER, INTENT(INOUT) ::    &
  ccb0    (row_length,rows)  &! Cnv Cld Base of highest layer on gridpoint
, cct0    (row_length,rows)  &! Cnv Cld Top  of highest layer on gridpoint
, lcbase0 (row_length,rows)   ! Cnv Cld Top  of lowest  layer on gridpoint

! Section 5 cloud properties:
!   Refer to diagnostics from calls to convection only.
!   Intent InOut because they are passed upto atmos_physics2
!   so they can be held until the last convection call before
!   exiting atmos_physics2


REAL, INTENT(INOUT) ::                &
  cca    (row_length,rows,n_cca_lev)  &! Cnv.Cld Amount (0-1)
, ccw    (tdims%i_end,tdims%j_end,tdims%k_end) &! Cnv.Cld Water (kg/kg)
, cclwp  (row_length,rows)            &! Cloud Condensed water path
                                       ! (kg/m^2)
, cca_2d(row_length,rows)             &! Cnv.Cld Amount (2d)
                                       ! with no anvil
, lcca   (row_length,rows)             ! Conv. Cloud Amount of
                                       ! lowest Convective Cloud Layer
                                       ! on gridpoint. With no anvil
                                       ! scheme modification (0-1)

INTEGER, INTENT(INOUT) ::  &!
  ccb    (row_length,rows) &! Cnv Cld Base of highest layer on gridpoint
, cct    (row_length,rows) &! Cnv Cld Top  of highest layer on gridpoint
, lcbase (row_length,rows) &! Cnv Cld Base of lowest  layer on gridpoint
, lctop  (row_length,rows)  ! Cnv Cld Top  of lowest  layer on gridpoint

INTEGER, INTENT(INOUT) ::                                         &
  error_code

!------------------------------------------------------------------------------
! End of subroutine argument definitions
!------------------------------------------------------------------------------
!---------------------------------------------------------------------------
! local variables.
!---------------------------------------------------------------------------

! loop counters
INTEGER ::                                                        &
  i, j, k, ii, jj                                                 &
, call_number

INTEGER ::                                                        &
  n_conv_levels                                                   &
, n_conv_points                                                   &
, this_point                                                      &
, info                                                            &
, step                                                            &
, seg_points                                                      &
, first_point

! CCRad - Local variables and parameters
REAL    :: cld_life_3d(row_length, rows, n_cca_lev)
REAL    :: decay_time
INTEGER :: ndecay_steps

! Local variables for convection prognostic calculations
REAL :: conv_active  ! = 1 on levels / points where the magnitude of the
                     ! convective theta tendency exceeds the threshold
                     ! dthetadt_conv_active_threshold,
                     ! = 0 elsewhere.
REAL :: theta_inc_threshold  ! dthetadt_conv_active_threshold converted to
                             ! increment over the convection timestep

! '_inc' - convective increments
REAL :: theta_inc      ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: q_inc          ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: qcl_inc        ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: qcf_inc        ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: cf_liquid_inc  ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: cf_frozen_inc  ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: bulk_cf_inc    ( tdims%i_end, tdims%j_end, tdims%k_end )

! Local copies of primary fields to update during convective substepping
REAL :: theta_conv     ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: q_conv         ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: qcl_conv       ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: qcf_conv       ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: qcf2_conv      ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: qrain_conv     ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: qgraup_conv    ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: cf_liquid_conv ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: cf_frozen_conv ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: bulk_cf_conv   ( tdims%i_end, tdims%j_end, tdims%k_end )
REAL :: u_conv         ( pdims%i_end, pdims%j_end, pdims%k_end )
REAL :: v_conv         ( pdims%i_end, pdims%j_end, pdims%k_end )

! To hold mean profiles over convection substeps until they
! are decayed on either each substep or at the end of a set of
! convection substeps(n_conv_calls)
INTEGER ::                         &
  ccb0_local    (row_length,rows)  &
, cct0_local    (row_length,rows)  &
, lcbase0_local (row_length,rows)

REAL ::                                              &
  cca0_local   (row_length,rows,n_cca_lev)           &
, cca0_local_dp   (row_length,rows,n_cca_lev)        &
, cca0_local_md   (row_length,rows,n_cca_lev)        &
, cca0_local_sh   (row_length,rows,n_cca_lev)        &
, ccw0_local   (tdims%i_end,tdims%j_end,tdims%k_end) &
, cclwp0_local (row_length,rows)


INTEGER ::                       &
  freeze_lev(row_length,rows)    & ! freezing level no
, it_kterm_deep(row_length,rows) & !lev no for terminating deep con
, it_kterm_shall(row_length,rows) & !lev no for terminating shallow con
, it_cg_term(row_length,rows)      !lev no for terminating congestus


LOGICAL ::                       &
  cumulus_1d(rowsrowlength)      & ! 1-d version of CUMULUS
 ,l_shallow_1d(rowsrowlength)    & ! 1-d version of SHALLOW_BL
 ,l_congestus_1d(rowsrowlength)  & ! 1-d version of congestus
 ,l_mid_1d(rowsrowlength)        & ! 1-d version of mid
 ,l_no_cumulus_1d(rowsrowlength)   ! 1-d version of no cumulus


! Total increments applied up to now across model timestep
! ie.  d(slow physics) + d(dyn),
REAL ::                                                                &
  u_inc_step(tdims%i_end,tdims%j_end,tdims%k_end)                      &
 ,v_inc_step(tdims%i_end,tdims%j_end,tdims%k_end)

LOGICAL ::       &
  l_tracer         ! Switch for tracer variables used


! Convection
REAL ::                                           &
  dthbydt(tdims%i_end,tdims%j_end,tdims%k_end)    & ! theta increment
, dqbydt(tdims%i_end,tdims%j_end,tdims%k_end)     & ! q increment
, dqclbydt(tdims%i_end,tdims%j_end,tdims%k_end)   & ! Q4 Increment qCL
, dqcfbydt(tdims%i_end,tdims%j_end,tdims%k_end)   & ! Q4 Increment qCF
, dcflbydt(tdims%i_end,tdims%j_end,tdims%k_end)   & ! Cloud Increment
, dcffbydt(tdims%i_end,tdims%j_end,tdims%k_end)   & ! Cloud Increment
, dbcfbydt(tdims%i_end,tdims%j_end,tdims%k_end)   & ! Cloud Increment
, dubydt_p(tdims%i_end,tdims%j_end,tdims%k_end)   & ! du/dt on p grid
, dvbydt_p(tdims%i_end,tdims%j_end,tdims%k_end)     ! dv/dt on p grid

REAL ::                                           &
  dq_add(tdims%i_end,tdims%j_end,tdims%k_end)      ! q added for safe profile

!  allocatable array holding copies of all tracers
REAL, ALLOCATABLE :: tot_tracer(:,:,:,:)

REAL, ALLOCATABLE ::    &
  dqt(:,:,:)              ! Total moisture increment for moisture checking


! Work arrays SCM diagnostic output
REAL :: TmpScm2d(ScmRowLen,ScmRow)
REAL :: TmpScm3d(ScmRowLen,ScmRow,tdims%k_end)

! Declare structure to store SCM convection sub-step diagnostics 
! passed up from the convection scheme
TYPE(scm_convss_dg_type), ALLOCATABLE :: scm_convss_dg(:,:)
! Flag for whether to calculate the SCM diagnostics contained in the
! above structure (and whether to allocate it at all).  This flag is
! passed down into the convection routines; we can set it in this
! routine to control which substeps the diagnostics are calculated on.
LOGICAL :: l_scm_convss_dg

LOGICAL :: l_extra_call = .TRUE.   ! true this is an additional
                                    ! call to conv_diag within a timestep

!  sizes for temp. tracer calculations
INTEGER ::                                                        &
  ntra_fld                                                        &
, ntra_lev

! Holding arrays for Section 0 cloud properties around convection substeps
REAL ::                                  &
  it_cca0   (row_length,rows,n_cca_lev)  &! Cnv.Cld Amount (0-1)
, it_cca0_dp (row_length,rows,n_cca_lev) &! Deep cnv.Cld Amount (0-1)
, it_cca0_md (row_length,rows,n_cca_lev) &! Mid-level cnv.Cld Amount (0-1)
, it_cca0_sh (row_length,rows,n_cca_lev) &! Shallow cnv.Cld Amount (0-1)
, it_ccw0   (tdims%i_end,tdims%j_end,tdims%k_end) &! Cnv.Cld Water (kg/kg)
, it_cclwp0 (row_length,rows)            &! Cloud Condensed water path
                                          ! (kg/m^2)
, it_cca0_2d(row_length,rows)             ! Cnv.Cld Amount (2d) with
                                          ! no anvil

INTEGER ::                     &
  it_ccb0   (row_length,rows)  &! Cnv Cld Base of highest layer on gridpoint
, it_cct0   (row_length,rows)  &! Cnv Cld Top  of highest layer on gridpoint
, it_lcbase0(row_length,rows)   ! Cnv Cld Base of lowest  layer on gridpoint


! Holding arrays for Section 5 cloud properties around convection substeps
REAL ::                                  &
  it_cca    (row_length,rows,n_cca_lev)  &! Cnv.Cld Amount (0-1)
, it_ccw    (tdims%i_end,tdims%j_end,tdims%k_end) &! Cnv.Cld Water (kg/kg)
, it_cclwp  (row_length,rows)            &! Cloud Condensed water path
                                          ! (kg/m^2)
, it_cca_2d (row_length,rows)            &! Cnv.Cld Amount (2d)
                                          ! with no anvil
, it_lcca   (row_length,rows)             ! Conv. Cloud Amount of
                                          ! lowest Convective Cloud Layer
                                          ! on gridpoint. With no anvil
                                          ! scheme modification (0-1)

INTEGER ::                     &
  it_ccb    (row_length,rows)  &! Cnv Cld Base of highest layer on gridpoint
, it_cct    (row_length,rows)  &! Cnv Cld Top  of highest layer on gridpoint
, it_lcbase (row_length,rows)  &! Cnv Cld Base of lowest  layer on gridpoint
, it_lctop  (row_length,rows)   ! Cnv Cld Top  of lowest  layer on gridpoin

! Arrays holding diagnostics around substeps
REAL  ::                                                          &
  it_conv_rain    (row_length, rows)                              &
, it_conv_snow    (row_length, rows)                              &
, it_conv_rain_3d (tdims%i_end,tdims%j_end,tdims%k_end)           &
, it_conv_snow_3d (tdims%i_end,tdims%j_end,tdims%k_end)           &
, it_cape_out(row_length, rows)                                   &
                                ! CAPE
, it_precip_dp(row_length, rows)                                  &
                                 ! deep precip
, it_precip_sh(row_length, rows)                                  &
                                 ! shallow precip
, it_precip_md(row_length, rows)                                  &
                                 ! mid-level precip
, it_precip_cg(row_length, rows)                                  &
                                 ! congestus precip
, it_wstar_dn(row_length, rows)                                   &
, it_wstar_up(row_length, rows)                                   &
, it_mb1(row_length, rows)                                        &
, it_mb2(row_length, rows)                                        &
, it_up_flux_half(pdims%i_end,pdims%j_end,pdims%k_end)            &
                                !up flux on half levs.
, ind_cape_reduced(row_length, rows)                              &
                             ! indicator of reduced cape timescale
, cape_ts_used(row_length, rows)                                  &
                             ! Actual cape timescale used for deep (s)
, it_ind_deep(row_length, rows)        & ! indicator deep actually occurs
, it_ind_shall(row_length, rows)       & ! indicator shallow actually occurs
, it_dp_cfl_limited(row_length, rows)  & ! Indicator of CFL limited deep
, it_md_cfl_limited(row_length, rows)    ! Indicator of CFL limited mid


LOGICAL :: l_mid(row_length, rows)
                          ! true if mid level convection is
                          ! possible on a convection substep

! u'w' and v'w' fluxes are on theta levels
REAL ::                                             &
  it_uw_dp(tdims%i_end,tdims%j_end,tdims%k_end)     &
, it_vw_dp(tdims%i_end,tdims%j_end,tdims%k_end)     &
, it_uw_shall(tdims%i_end,tdims%j_end,tdims%k_end)  &
, it_vw_shall(tdims%i_end,tdims%j_end,tdims%k_end)  &
, it_uw_mid(tdims%i_end,tdims%j_end,tdims%k_end)    &
, it_vw_mid(tdims%i_end,tdims%j_end,tdims%k_end)

! Mass flux, entrain etc in our model is calculated on theta levels
REAL ::                                                                &
  it_up_flux(tdims%i_end,tdims%j_end,tdims%k_end)                      &
, it_dwn_flux(tdims%i_end,tdims%j_end,tdims%k_end)                     &
, it_entrain_up(tdims%i_end,tdims%j_end,tdims%k_end)                   &
, it_detrain_up(tdims%i_end,tdims%j_end,tdims%k_end)                   &
, it_entrain_dwn(tdims%i_end,tdims%j_end,tdims%k_end)                  &
, it_detrain_dwn(tdims%i_end,tdims%j_end,tdims%k_end)                  &
, r_rho_lev(pdims%i_end,pdims%j_end,pdims%k_end)                       &
, r_theta_lev(tdims%i_end,tdims%j_end,0:tdims%k_end)                   &
, exner_rho_lev(pdims%i_end,pdims%j_end,pdims%k_end)                   &
, zh_copy(row_length, rows)                                            &
, one_over_conv_calls                                                  &
, fraction_step                                                        &
, timestep_conv

REAL ::                               &! (Parcel vertical velocity)^2 [(m/s)^2
  it_w2p  (tdims%i_start:tdims%i_end, &! calculated on convection sub-step
           tdims%j_start:tdims%j_end, &
                       1:tdims%k_end)

! Compressed Surface precipitation based 3d convective prognostic in kg/m2/s
REAL :: conv_prog_precip_c(tdims%i_start:tdims%i_end,                  &
                           tdims%j_start:tdims%j_end,                  &
                           1:tdims%k_end)

! Downdraught increments as rates
REAL ::                                                                &
  it_dt_dd(tdims%i_end,tdims%j_end,tdims%k_end)                        &
, it_dq_dd(tdims%i_end,tdims%j_end,tdims%k_end)                        &
, it_du_dd(tdims%i_end,tdims%j_end,tdims%k_end)                        &
, it_dv_dd(tdims%i_end,tdims%j_end,tdims%k_end)

! Fractional areas updraught and downdraught
REAL ::                                                                &
  it_area_ud(tdims%i_end,tdims%j_end,tdims%k_end)                      &
, it_area_dd(tdims%i_end,tdims%j_end,tdims%k_end)

INTEGER ::                     &
  nlcl(row_length, rows)         ! dummy array for conv_diag call

!-----------------------------------------------------------------------
! More diagnostic arrays holding values for a substep
REAL ::                                                            &
  it_mf_deep(tdims%i_end,tdims%j_end,tdims%k_end)                  &
, it_mf_congest (tdims%i_end,tdims%j_end,tdims%k_end)              &
, it_mf_shall(tdims%i_end,tdims%j_end,tdims%k_end)                 &
, it_mf_midlev(tdims%i_end,tdims%j_end,tdims%k_end)                &
, it_dt_deep(tdims%i_end,tdims%j_end,tdims%k_end)                  &
, it_dt_congest(tdims%i_end,tdims%j_end,tdims%k_end)               &
, it_dt_shall(tdims%i_end,tdims%j_end,tdims%k_end)                 &
, it_dt_midlev(tdims%i_end,tdims%j_end,tdims%k_end)                &
, it_dq_deep(tdims%i_end,tdims%j_end,tdims%k_end)                  &
, it_dq_congest(tdims%i_end,tdims%j_end,tdims%k_end)               &
, it_dq_shall(tdims%i_end,tdims%j_end,tdims%k_end)                 &
, it_dq_midlev (tdims%i_end,tdims%j_end,tdims%k_end)               &
, it_du_deep (tdims%i_end,tdims%j_end,tdims%k_end)                 &
, it_du_congest(tdims%i_end,tdims%j_end,tdims%k_end)               &
, it_du_shall(tdims%i_end,tdims%j_end,tdims%k_end)                 &
, it_du_midlev (tdims%i_end,tdims%j_end,tdims%k_end)               &
, it_dv_deep (tdims%i_end,tdims%j_end,tdims%k_end)                 &
, it_dv_congest(tdims%i_end,tdims%j_end,tdims%k_end)               &
, it_dv_shall (tdims%i_end,tdims%j_end,tdims%k_end)                &
, it_dv_midlev (tdims%i_end,tdims%j_end,tdims%k_end)

! These fluxes are on rho levels
REAL ::                                                            &
 it_wqt_flux(pdims%i_end,pdims%j_end,tdims%k_end)                  &
,it_wql_flux(pdims%i_end,pdims%j_end,tdims%k_end)                  &
,it_wthetal_flux(pdims%i_end,pdims%j_end,tdims%k_end)              &
,it_wthetav_flux(pdims%i_end,pdims%j_end,tdims%k_end)

! The '_conv' arrays hold the prognostic values input to the convection scheme.
! The convection scheme requires all the fields on the p_grid in the
! horizontal. Values held in _conv arrays depend on the choice of sub-stepping
! method and are never the start of time step prognostic values held in the
! '_n' arrays.

! If l_rediagnosis = T: '_conv' = '_n' + (1/n_conv_calls)*(slow phys+dyn incs)
! If l_rediagnosis = F: '_conv' = '_n' + (slow phys+dyn incs) = '_star'

REAL ::                                                                 &
  ep4(row_length, rows)             ! diagnostic Evap/Precip

!----------------------------------------------------------------------
! Arrays to hold total increments to T, q (start of time step values) etc
! before calling convection. These increments are from the slow physics
! (atmos_physics1) and from semi-lagrangian dynamics.
! Note total wind increments are already held in r_u_p and r_v_p arrays.
! Want to know all these increments so they can be applied gradually
! over sub-steps in convection. Assuming l_rediagnosis = .true.
!----------------------------------------------------------------------
REAL, ALLOCATABLE ::         &
  theta_inc_step(:,:,:)      &
 ,q_inc_step(:,:,:)          &
 ,qcl_inc_step(:,:,:)        &
 ,qcf_inc_step(:,:,:)        &
 ,qrain_inc_step(:,:,:)      &
 ,qgraup_inc_step(:,:,:)     &
 ,qcf2_inc_step(:,:,:)       &
 ,cf_liquid_inc_step(:,:,:)  &
 ,cf_frozen_inc_step(:,:,:)  &
 ,bulk_cf_inc_step(:,:,:)

! Segmentation variables
TYPE(segment_type),      ALLOCATABLE  :: segments(:)
TYPE(meta_segment_type)               :: meta_segments

INTEGER :: n_cumulus   ! Number of CUMULUS points
INTEGER :: n_deep      ! Number of DEEP points
INTEGER :: n_shallow   ! Number of SHALLOW points
INTEGER :: n_congestus ! Number of congestus points
INTEGER :: n_mid       ! Number of mid-level points

! convective history
REAL ::              &
  decay_amount       &  ! decay fraction
, tot_conv_precip       ! total convective precip rate

REAL ::              &
  tot_conv_precip_2d(row_length, rows) ! 2D total convective precip rate

! Several options are available:
INTEGER, PARAMETER :: pc2_conv_original = 1
! As used in Wilson et al (2008). Condensate and cloud fraction
! increments are calculated independently and there is no checking that
! the implied in-cloud condensate amounts are sensible.
!
! i_pc2_conv_coupling /= pc2_conv_original
! Protect against generation of inconsistently low cloud
! fraction implying very high in-cloud condensate amounts.
! If the in-cloud condensate amount is about to be made greater than
! 2.0e-3 then the cloud fraction is increased (up to a value of 1.0)

REAL    :: tmpb4 ! Temporary storing space for cloud fraction before
                 ! making consistency checks.

!-----------------------------------------------------------------------

! Error reporting variables
INTEGER ::                                                        &
  errorstatus

CHARACTER (LEN=errormessagelength) ::                             &
  cmessage

CHARACTER (LEN=*), PARAMETER ::                                   &
  RoutineName='NI_CONV_CTL'

CHARACTER (LEN=80) ::                                             &
  scheme_name

!Automatic segment size tuning
TYPE(autotune_type), ALLOCATABLE, SAVE :: autotune_state
INTEGER :: segment_size

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!Segmentation
INTEGER :: num_parallel
INTEGER :: ipar

!=============================================================================
! Start of convection control routine
!=============================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------------
IF (error_code  ==  0 ) THEN     ! continue with convection

  !-----------------------------------------------------------------------------
  ! Section 1 - initialisation and setup
  !-----------------------------------------------------------------------------

  !-------------------------------------
  ! 1.1 Initialise local cloud arrays
  !-------------------------------------
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(k, i, j)                                &
!$OMP& SHARED(rows,row_length,cclwp0_local,ccb0_local,cct0_local,n_cca_lev,   &
!$OMP&    ccw0_local,dqbydt,dqclbydt,dqcfbydt,cca0_local,tdims,l_cca_dp_prog, &
!$OMP&    cca0_local_dp,l_cca_md_prog,cca0_local_md,l_cca_sh_prog,            &
!$OMP&    cca0_local_sh)

!$OMP DO SCHEDULE(STATIC)  
  DO j=1, rows
    DO i=1, row_length
      cclwp0_local(i,j) = 0.0
      ccb0_local(i,j)   = 0
      cct0_local(i,j)   = 0
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)  
  DO k=1, n_cca_lev
    DO j=1, rows
      DO i=1, row_length
        cca0_local(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)  
  DO k=1, tdims%k_end
    DO j=1, rows
      DO i=1, row_length
        ccw0_local(i,j,k) = 0.0
        ! for q checking
        dqbydt(i,j,k)    = 0.0
        dqclbydt(i,j,k)  = 0.0
        dqcfbydt(i,j,k)  = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  IF (l_cca_dp_prog) THEN
!$OMP DO SCHEDULE(STATIC)  
    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          cca0_local_dp(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_cca_md_prog) THEN
!$OMP DO SCHEDULE(STATIC)  
    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          cca0_local_md(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  IF (l_cca_sh_prog) THEN
!$OMP DO SCHEDULE(STATIC)  
    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          cca0_local_sh(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

!$OMP END PARALLEL

  !-----------------------------------------------------------------------
  ! 1.3 Convection time step information and level number
  !-----------------------------------------------------------------------

  timestep_conv=timestep/(n_conv_calls*1.0)

  ! Sub-timestep scheme

  fraction_step = 1.0/(REAL(n_conv_calls))

  one_over_conv_calls = 1.0/(n_conv_calls*1.0)

  ! Set number of levels to call convection for

  n_conv_levels = tdims%k_end

  IF (l_mom) THEN
    ! Limit convection calling levels to maximum of model_levels - 1
    ! This is because CMT increments to u and v exist on n_levels + 1
    IF (n_conv_levels  >   model_levels - 1 ) THEN
      n_conv_levels = model_levels - 1
    END IF
  END IF

  !-----------------------------------------------------------------------
  ! 1.4  Initialise arrays
  !-----------------------------------------------------------------------

  ! Initialise arrays for convection - without halos

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(k, i, j)                        &
!$OMP& SHARED(tdims,pdims_s,rows,row_length,                          &
!$OMP&     r_rho_lev,r_rho_levels,                                    &
!$OMP&     r_theta_lev,r_theta_levels,pdims,exner_rho_lev,            &
!$OMP&     exner_rho_levels, theta_inc, q_inc, qcl_inc, qcf_inc,      &
!$OMP&     cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,                 &
!$OMP&     l_mom, dubydt_pout, dvbydt_pout,                           &
!$OMP&     conv_rain, conv_snow, l_mid)

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        r_rho_lev(i,j,k) = r_rho_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = 0, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        r_theta_lev(i,j,k) = r_theta_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = pdims%k_start, pdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        exner_rho_lev(i,j,k) = exner_rho_levels(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

! Initialise convection increments to zero

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end
    DO j = 1, rows
      DO i = 1, row_length
        theta_inc(i,j,k) = 0.0
        q_inc(i,j,k) = 0.0
        qcl_inc(i,j,k) = 0.0
        qcf_inc(i,j,k) = 0.0
        cf_liquid_inc(i,j,k) = 0.0
        cf_frozen_inc(i,j,k) = 0.0
        bulk_cf_inc(i,j,k)   = 0.0
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  ! If convective momentum transport is on, initialise wind tendencies to zero
  IF (l_mom) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = pdims_s%k_start, pdims_s%k_end
      DO j = pdims_s%j_start, pdims_s%j_end
        DO i = pdims_s%i_start, pdims_s%i_end
          dubydt_pout(i,j,k) = 0.0
          dvbydt_pout(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

  ! Initialise output convective precip fields to zero

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      conv_rain(i,j) = 0.0
      conv_snow(i,j) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

  !---------------------------------------------------------------------
  ! Initialise flag to indicate whether mid-level convection is possible
  ! On first convection substep it is possible everywhere.
  !---------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO j=1, rows
    DO i=1, row_length
      l_mid(i,j) = .TRUE.
    END DO  !i
  END DO    !j
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  !---------------------------------------------------------------------
  ! Initialise 3d precipitation based prognostic without halos
  !---------------------------------------------------------------------
  IF (l_conv_prog_precip) THEN
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, rows, row_length, conv_prog_precip_c,             &
!$OMP         conv_prog_precip )                                       &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO  j = 1, rows
        DO i = 1, row_length
          conv_prog_precip_c(i,j,k) = conv_prog_precip(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  !---------------------------------------------------------------------
  ! 1.5 Set values of theta_conv, q_conv etc for input to convection
  !-----------------------------------------------------------------------
  IF (l_rediagnosis .AND. n_conv_calls > 1) THEN

    ! Allocate arrays to hold increments
    ALLOCATE ( theta_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    ALLOCATE ( q_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    ALLOCATE ( qcl_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    ALLOCATE ( qcf_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    ALLOCATE ( cf_liquid_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    ALLOCATE ( cf_frozen_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    ALLOCATE ( bulk_cf_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    IF (l_mcr_qrain) THEN
      ALLOCATE ( qrain_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    END IF
    IF (l_mcr_qgraup) THEN
      ALLOCATE ( qgraup_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    END IF
    IF (l_mcr_qcf2) THEN
      ALLOCATE ( qcf2_inc_step(tdims%i_end,tdims%j_end,tdims%k_end) )
    END IF


    !-----------------------------------------------------------------------
    ! Fill arrays
    !-----------------------------------------------------------------------
    ! Work out timestep increments so far i.e. from slow physics and
    ! semi-lagrangian dynamics. Note may be possible in future to pass
    ! these in directly from atm_step rather than pass in theta_star etc.
    !-----------------------------------------------------------------------
    DO k=1,tdims%k_end
      DO j=1,rows
        DO i=1,row_length
          theta_inc_step(i,j,k) = theta_star(i,j,k)-theta_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    DO k=1,tdims%k_end
      DO j=1,rows
        DO i=1,row_length
          q_inc_step(i,j,k)   = q_star(i,j,k)   - q_n(i,j,k)
          qcl_inc_step(i,j,k) = qcl_star(i,j,k) - qcl_n(i,j,k)
          qcf_inc_step(i,j,k) = qcf_star(i,j,k) - qcf_n(i,j,k)
          cf_liquid_inc_step(i,j,k) = cf_liquid_star(i,j,k)                  &
                                                 -cf_liquid_n(i,j,k)
          cf_frozen_inc_step(i,j,k) = cf_frozen_star(i,j,k)                  &
                                                 -cf_frozen_n(i,j,k)
          bulk_cf_inc_step(i,j,k)   = bulk_cf_star(i,j,k)                    &
                                                 -bulk_cf_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    ! Extra prognostics if present
    IF (l_mcr_qrain) THEN
      DO k=1,tdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qrain_inc_step(i,j,k)   = qrain_star(i,j,k)   - qrain_n(i,j,k)
            qrain_conv(i,j,k)   = qrain_n(i,j,k)                            &
                                   + fraction_step*qrain_inc_step(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
    IF (l_mcr_qgraup) THEN
      DO k=1,tdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qgraup_inc_step(i,j,k)   = qgraup_star(i,j,k) - qgraup_n(i,j,k)
            qgraup_conv(i,j,k)   = qgraup_n(i,j,k)                            &
                                   + fraction_step*qgraup_inc_step(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF
    IF (l_mcr_qcf2) THEN
      DO k=1,tdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qcf2_inc_step(i,j,k)   = qcf2_star(i,j,k) - qcf2_n(i,j,k)
            qcf2_conv(i,j,k)   = qcf2_n(i,j,k)                            &
                                   + fraction_step*qcf2_inc_step(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
    END IF


    !-----------------------------------------------------------------------
    ! Set up input prognostic values for convection - sub-step 1
    !-----------------------------------------------------------------------

    ! We want theta_conv = theta_n(start_step)+fraction_step*theta_inc_step(i,j,k)

    DO k=1,tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          theta_conv(i,j,k) = theta_n(i,j,k)                                &
                            +fraction_step*theta_inc_step(i,j,k)
        END DO ! i
      END DO ! j
    END DO

    DO k=1,tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          q_conv(i,j,k)   = q_n(i,j,k)  + fraction_step*q_inc_step(i,j,k)

          qcl_conv(i,j,k) = qcl_n(i,j,k)                                    &
                                        + fraction_step*qcl_inc_step(i,j,k)
          qcf_conv(i,j,k) = qcf_n(i,j,k)                                    &
                                        + fraction_step*qcf_inc_step(i,j,k)
          cf_liquid_conv(i,j,k) = cf_liquid_n(i,j,k)                        &
                                      + fraction_step*cf_liquid_inc_step(i,j,k)
          cf_frozen_conv(i,j,k) = cf_frozen_n(i,j,k)                        &
                                      + fraction_step*cf_frozen_inc_step(i,j,k)
          bulk_cf_conv(i,j,k) = bulk_cf_n(i,j,k)                            &
                                      + fraction_step*bulk_cf_inc_step(i,j,k)
        END DO ! i
      END DO ! j
    END DO

  ELSE     ! non rediagnosis option
    !-----------------------------------------------------------------------
    ! Original code uses all increments added to start of time step for
    ! sweep one
    !-----------------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(k, i, j)                       &
!$OMP& SHARED(tdims,rows,row_length,theta_conv,theta_star,           &
!$OMP&     q_conv,q_star,qcl_conv,qcl_star,qcf_conv,qcf_star,        &
!$OMP&     cf_liquid_conv,cf_liquid_star,cf_frozen_conv,             &
!$OMP&     cf_frozen_star,bulk_cf_conv,bulk_cf_star,l_mcr_qrain,     &
!$OMP&     qrain_conv,qrain_star,l_mcr_qgraup,qgraup_conv,           &
!$OMP&     qgraup_star,l_mcr_qcf2,qcf2_conv,qcf2_star)

!$OMP DO SCHEDULE(STATIC)
    DO k=1,tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          theta_conv(i,j,k) = theta_star(i,j,k)
        END DO ! i
      END DO ! j
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
    DO k=1,tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          q_conv(i,j,k)   = q_star(i,j,k)
          qcl_conv(i,j,k) = qcl_star(i,j,k)
          qcf_conv(i,j,k) = qcf_star(i,j,k)
          cf_liquid_conv(i,j,k) = cf_liquid_star(i,j,k)
          cf_frozen_conv(i,j,k) = cf_frozen_star(i,j,k)
          bulk_cf_conv(i,j,k)   = bulk_cf_star(i,j,k)
        END DO ! i
      END DO ! j
    END DO
!$OMP END DO NOWAIT

        ! Extra prognostics if present
    IF (l_mcr_qrain) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=1,tdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qrain_conv(i,j,k)   = qrain_star(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
!$OMP END DO NOWAIT
    END IF
    IF (l_mcr_qgraup) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=1,tdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qgraup_conv(i,j,k)   = qgraup_star(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
!$OMP END DO NOWAIT
    END IF
    IF (l_mcr_qcf2) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k=1,tdims%k_end
        DO j=1,rows
          DO i=1,row_length
            qcf2_conv(i,j,k)   = qcf2_star(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
!$OMP END DO NOWAIT
    END IF

!$OMP END PARALLEL

  END IF   ! test on l_rediagnosis

  ! ----------------------------------------------------------------------
  ! 1.6 Check initial input for convection ok if PC2
  ! ----------------------------------------------------------------------
  ! Prevent negative condensate problems by condensing vapour if needed.
  ! Input fields are updated without updating increments so change will
  ! only affect PC2 increments and not the top-level QX_STAR fields.
  ! ----------------------------------------------------------------------
  ! L_calc_dxek_if0:
  IF (l_calc_dxek) THEN

    CALL conv_pc2_init(rows, row_length, exner_theta_levels,            &
                       theta_conv, q_conv, qcl_conv, qcf_conv,          &
                       cf_liquid_conv, cf_frozen_conv, bulk_cf_conv)

  END IF  ! L_calc_dxek_if0

  ! ----------------------------------------------------------------------
  ! 1.7 Check for negative (less than a minimum) q being passed to convection
  ! ----------------------------------------------------------------------
  ! Note above PC2 code checking no negative qcl and qcf.
  ! Nothing currently preventing negative q going to convection
  ! Check for negative q if opt for safer convection or rediagnosis
  ! NOTE - NOT checking extra moist prognostics for negative values!
  !---------------------------------------------------------------------

  IF (l_safe_conv .OR. l_rediagnosis) THEN
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(k, i, j)                         &
!$OMP& SHARED(tdims,rows,row_length,q_conv,printstatus,                   &
!$OMP&     q_n,timestep_number,dq_add,mype)
    DO k=1,tdims%k_end
      DO j=1,rows
        DO i=1,row_length
          IF (q_conv(i,j,k) < qmin_conv) THEN
            IF (printstatus >= prstatus_normal) THEN
              ! Formatting for floating point double precision to get the same
              ! as unformatted print
              WRITE(umMessage,'(a,3(i0,1x),a,g26.18,a,g26.18,a,i0,1x,i0)') &
                    ' Negative q at start of convection: i,j,k ',i,j,k,    &
                    ' q_conv ', q_conv(i,j,k),' q_n ',q_n(i,j,k),' step ', &
                      timestep_number,mype
              CALL umPrint(umMessage,src='ni_conv_ctl')
            END IF
            ! store q added to ensure sensible profile
            dq_add(i,j,k) = qmin_conv - q_conv(i,j,k)
            ! Reset to qmin in non-conservative way
            q_conv(i,j,k) = qmin_conv
          ELSE
            dq_add(i,j,k) = 0.0
          END IF
        END DO ! i
      END DO ! j
    END DO ! k
!$OMP END PARALLEL DO 
  END IF  ! l_rediagnosis

  !-----------------------------------------------------------------------
  ! 1.9 Convective momemtum transport - l_mom=.true.
  !-----------------------------------------------------------------------
  ! U and V required by convection but on the p-grid.

  ! Convection uses wind at beginning of time step plus increments
  ! For new sub-stepping option require only part of increment to be
  ! added each sub-step

  ! Are U and V already on p grid in atmos_physics2? If in future this is
  ! true just pass down and only need to interpolate R_u etc
  !----------------------------------------------------------------------

  IF (l_mom) THEN

    SELECT CASE (model_type)
    CASE DEFAULT

      IF (l_rediagnosis) THEN

        DO k = 1,tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              u_conv(i,j,k) = u_p(i,j,k)+fraction_step*r_u_p(i,j,k)
              v_conv(i,j,k) = v_p(i,j,k)+fraction_step*r_v_p(i,j,k)
              ! Need to initialise u_inc_step and v_inc_step
              u_inc_step(i,j,k) = r_u_p(i,j,k)
              v_inc_step(i,j,k) = r_v_p(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k

      ELSE       ! original code

!$OMP  PARALLEL DO DEFAULT(NONE) SHARED(u_conv, v_conv, ustar_p, tdims, &
!$OMP& vstar_p, rows, row_length) PRIVATE(i, j, k)
        DO k = 1, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              u_conv(i,j,k) = ustar_p(i,j,k)
              v_conv(i,j,k) = vstar_p(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k
!$OMP END PARALLEL DO

      END IF   ! test on l_rediagnosis

    CASE (mt_single_column)

      IF (l_rediagnosis) THEN

        ! SCM u and v on p grid  (No problems with v having less rows)
        DO k = 1, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              u_conv(i,j,k) = u_p(i,j,k) + fraction_step*r_u_p(i,j,k)
              v_conv(i,j,k) = v_p(i,j,k) + fraction_step*r_v_p(i,j,k)
              ! Need to initialise u_inc_step and v_inc_step
              u_inc_step(i,j,k) = r_u_p(i,j,k)
              v_inc_step(i,j,k) = r_v_p(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k

      ELSE

        ! SCM u and v on p grid
        DO k = 1, tdims%k_end
          DO j = 1, rows
            DO i = 1, row_length
              u_conv(i,j,k) = ustar_p(i,j,k)
              v_conv(i,j,k) = vstar_p(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k

      END IF    ! test on l_rediagnosis
        
    END SELECT ! model_type

  END IF   !(l_mom)

  !-----------------------------------------------------------------------
  ! Section 1.10  Setup total tracer variables - Aerosol scheme or UKCA
  !-----------------------------------------------------------------------
  l_tracer = ( ( cycleno == numcycles ) .AND.                &
               ( l_soot       .OR.  l_co2_interactive  .OR.  &
                 l_dust       .OR.  l_biomass          .OR.  &
                 l_sulpc_so2  .OR.  l_use_cariolle     .OR.  &
                 l_ocff       .OR.  l_nitrate          .OR.  &
                 l_murk_conv  .OR.  tr_ukca > 0        .OR.  &
                 tr_vars > 0 ) )

  ! work with tracers only in final cycle

  IF ( cycleno == numcycles ) THEN

    ! Set up array dimensions for total tracer array (free + sulphur cycle
    ! tracers) so that convective transport of all tracers is done
    
    CALL tracer_total_var(tr_vars, tr_ukca, ntra_fld, ntra_lev )

    ! Allocate the space in tot_tracer
    ALLOCATE( tot_tracer(row_length, rows, ntra_lev, ntra_fld) )

    ! Copy the required tracers to tot_tracer array
    CALL tracer_copy(row_length, rows, ntra_fld, ntra_lev, tr_vars, tr_ukca,&
           aerosol,                                                         &
           dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,     &
           so2, so4_aitken, so4_accu, so4_diss,                             &
           dms, nh3, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd,    &
           bmass_cld, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss,    &
           co2,  free_tracers, tracer_ukca,                                 &
           ozone_tracer,                                                    &
           tot_tracer ) 

  ELSE
    ntra_fld = 1             ! can't have 0 sized arrays
    ntra_lev = 1
    ! Allocate dummy space in tot_tracer
    ALLOCATE( tot_tracer(row_length, rows, ntra_lev, ntra_fld) )
  END IF       ! test on cycleno == numcycles


  !--------------------------------------------------------------------------
  ! 1.11 Choices based on whether first call to ni_conv_ctl this model physics
  ! timestep
  !--------------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! Set up cloud decay lifetime
  !---------------------------------------------------------------------

  IF (rad_cloud_decay_opt /= rad_decay_off) THEN

    SELECT CASE(rad_cloud_decay_opt)

    CASE (rad_decay_conv_substep)
      decay_time   = timestep_conv
      ndecay_steps = 1

    CASE (rad_decay_full_timestep)
      decay_time   = timestep
      ndecay_steps = n_conv_calls

    END SELECT

    SELECT CASE(cld_life_opt)
    CASE (cld_life_func_hgt)
      !---------------------------------------------------------------
      ! Make lifetime a function of cloud size by making it a function
      ! of height: 30 minutes for shallow, 2 hours for high cloud
      !---------------------------------------------------------------
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(k, i, j)                    &
!$OMP& SHARED(n_cca_lev,rows,row_length,cld_life_3d,decay_time,      &
!$OMP&     fixed_cld_life,z_theta)
      DO k=1, n_cca_lev
        DO j=1, rows
          DO i=1, row_length
            cld_life_3d(i,j,k) = decay_time / (1800.0                 &
                               + (fixed_cld_life*0.5 - 900.0)         &
                               * (TANH((z_theta(i,j,k)/1000.0)-5.0)   &
                               + 1.0))
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

    CASE (cld_life_constant)
      !---------------------------------------------------------------
      ! Set Cloud_lifetime to a constant
      !---------------------------------------------------------------
!$OMP  PARALLEL DO DEFAULT(NONE) PRIVATE(k, i, j)                   &
!$OMP& SHARED(n_cca_lev,rows,row_length,cld_life_3d,decay_time,     &
!$OMP&     fixed_cld_life)
      DO k=1, n_cca_lev
        DO j=1, rows
          DO i=1, row_length
            cld_life_3d(i,j,k) = decay_time / fixed_cld_life
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO
    END SELECT    ! cld_life_opt

  END IF        !  Rad_cloud_decay_opt


  ! Set flag for Single-Column Model convection sub-step diagnostics
  l_scm_convss_dg = .FALSE.
  IF (model_type == mt_single_column) THEN
    IF (l_scmdiags(scmdiag_convss)) THEN
      ! We have the option of switching these diagnostics off on all but the
      ! first sub-step, but for now set to true for both sub-steps.
      l_scm_convss_dg = .TRUE.
    END IF
  END IF

  ! Allocate arrays for SCM convection sub-step diagnostics
  ! (done before substepping so that we only allocate / deallocate once)

  ! Allocate the array of structures for each SCM point (usually just
  ! (1,1) in the single-column model, but don't want to preclude
  ! multi-column functionality)
  ! Note: this is allocated in all runs since the call to glue_conv below
  ! requires the outer object to be indexed
  ALLOCATE( scm_convss_dg(row_length,rows) )

  IF (l_scm_convss_dg) THEN

    ! Allocate the fields within the structure for each grid-column
    DO j = 1, rows
      DO i = 1, row_length
        CALL scm_convss_dg_allocate( scm_convss_dg(i,j) )
      END DO
    END DO

  END IF


  !==========================================================================
  ! Section 2 Call Convection scheme.
  !==========================================================================
  !============================================================================
  !  Convection sub-stepping loop
  !============================================================================
  DO call_number = 1, n_conv_calls

    !--------------------------------------------------------------------------
    ! 2.1 initialisation of work arrays for convection substeps
    !--------------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(k, i, j)                        &
!$OMP& SHARED(it_ccw,it_ccw0,it_conv_rain_3d,it_conv_snow_3d,it_cca,    &
!$OMP&    it_cca0,it_w2p,it_lcca,it_lcbase,it_lctop,it_ccb,it_cct,      &
!$OMP&    it_cca_2d,it_cclwp,it_ccb0,it_cct0,it_cca0_2d,it_cclwp0,      &
!$OMP&    it_conv_rain,it_conv_snow,it_precip_dp,it_precip_sh,          &
!$OMP&    it_precip_md,it_cape_out,it_kterm_deep,it_kterm_shall,        &
!$OMP&    it_mid_level,it_dp_cfl_limited,it_md_cfl_limited,it_cca0_dp,  &
!$OMP&    it_cca0_md,it_cca0_sh, it_precip_cg, it_wstar_up, it_mb1,     &
!$OMP&    it_mb2,it_cg_term,rows,row_length,n_cca_lev,tdims,            &
!$OMP&    i_convection_vn,flg_w_eqn,l_cca_dp_prog,l_cca_md_prog,        &
!$OMP&    l_cca_sh_prog)

!$OMP DO SCHEDULE(STATIC)  
    DO k=1, tdims%k_end
      DO j=1, rows
        DO i=1, row_length
          it_ccw          (i,j,k) = 0.0
          it_ccw0         (i,j,k) = 0.0
          it_conv_rain_3d (i,j,k) = 0.0
          it_conv_snow_3d (i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)  
    DO k=1, n_cca_lev
      DO j=1, rows
        DO i=1, row_length
          it_cca  (i,j,k) = 0.0
          it_cca0 (i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    IF (flg_w_eqn) THEN
!$OMP DO SCHEDULE(STATIC)  
      DO k=1, tdims%k_end
        DO j=tdims%j_start, tdims%j_end
          DO i=tdims%i_start, tdims%i_end
            it_w2p(i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

!$OMP DO SCHEDULE(STATIC)  
    DO j = 1, rows
      DO i = 1, row_length
        it_lcca(i,j)   = 0.0
        it_lcbase(i,j) = 0
        it_lctop(i,j)  = 0

        it_ccb(i,j)    = 0
        it_cct(i,j)    = 0
        it_cca_2d(i,j) = 0.0
        it_cclwp(i,j)  = 0.0

        it_ccb0(i,j)   = 0
        it_cct0(i,j)   = 0
        it_cca0_2d(i,j)= 0.0
        it_cclwp0(i,j) = 0.0

        it_conv_rain(i,j) = 0.0
        it_conv_snow(i,j) = 0.0
        it_precip_dp(i,j) = 0.0
        it_precip_sh(i,j) = 0.0
        it_precip_md(i,j) = 0.0
        it_cape_out(i,j) = 0.0
        it_kterm_deep(i,j) = 0
        it_kterm_shall(i,j) = 0
        it_mid_level(i,j) = .FALSE.
        it_dp_cfl_limited(i,j) = 0.0
        it_md_cfl_limited(i,j) = 0.0

        it_precip_cg(i,j) = 0.0
        it_wstar_up(i,j) = 0.0
        it_mb1(i,j) = 0.0
        it_mb2(i,j) = 0.0
        it_cg_term(i,j) = 0
      END DO
    END DO
!$OMP END DO NOWAIT

    IF (l_cca_dp_prog) THEN
!$OMP DO SCHEDULE(STATIC)  
      DO k=1, n_cca_lev
        DO j=1, rows
          DO i=1, row_length
            it_cca0_dp (i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    IF (l_cca_md_prog) THEN
!$OMP DO SCHEDULE(STATIC)  
      DO k=1, n_cca_lev
        DO j=1, rows
          DO i=1, row_length
            it_cca0_md (i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    IF (l_cca_sh_prog) THEN
!$OMP DO SCHEDULE(STATIC)  
      DO k=1, n_cca_lev
        DO j=1, rows
          DO i=1, row_length
            it_cca0_sh (i,j,k) = 0.0
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

!$OMP END PARALLEL


    ! Intialise SCM convection sub-step diagnostics to be passed up from
    ! the convection scheme to zero (needs to be re-zeroed each sub-step)
    IF (l_scm_convss_dg) THEN
      CALL scm_convss_dg_initzero( scm_convss_dg )
    END IF


    !-----------------------------------------------------------------------------
    ! NB: increments to t and q are added on inside routine.
    ! Increments to qCL, qCF, CFl, CFf are calculated but only added at the
    ! control level (see below).
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! 2.2 segment the call to convection, to reduce the memory required in
    ! the subroutine convect.
    !-----------------------------------------------------------------------------


       !--------------------------------------------------------------------------
       ! 2.3  Rediagnosis ?
       !--------------------------------------------------------------------------

    IF (l_rediagnosis) THEN
      !-------------------------------------------------------------------------
      ! Do we want to rediagnose where shallow and deep convection occurs on
      ! subsequent convection sub-steps?
      !    call_number = 1  - conv_diag called from atmos_physics2 to
      !                       determine cumulus_1d, l_shallow etc
      !    call_number > 1  - conv_diag call here to rediagnose cumulus etc
      !-------------------------------------------------------------------------

      IF (call_number > 1) THEN

        ! Take copy as do not want to overwrite otherwise use
        ! value from previous sweep already held in zh_copy
        IF (call_number == 2) THEN
          DO j=1,rows
            DO i=1,row_length
              zh_copy(i,j) = zh(i,j)
            END DO
          END DO
        END IF
        !  Initialise conv_diag output arrays
        DO j = 1, rows
          DO i = 1, row_length
            ! ntml(i,j)=1     May want previous value
            ntpar(i,j)=1
            nlcl(i,j)=1
            cumulus(i,j)=.FALSE.
            l_shallow(i,j)=.FALSE.
            delthvu(i,j)=0.0
            ql_ad(i,j) = 0.0
            ztop(i,j)=0.0
            dzh(i,j) =-9.9e9     ! initialise to large and negative
            qcl_inv_top(i,j) = 0.0
            zlcl(i,j)=0.0
            zlcl_uv(i,j)=0.0
            wstar(i,j)=0.0
            wthvs(i,j)=0.0
          END DO
        END DO


        ! --------------------------------------------------------------------
        ! Prevent negative condensate problems by condensing vapour if needed.
        ! Input fields are updated without updating increments so change will
        ! only affect PC2 increments and not the top-level QX_STAR fields.
        ! --------------------------------------------------------------------
        ! L_calc_dxek_if0:
        IF (l_calc_dxek) THEN

          CALL conv_pc2_init(rows, row_length, exner_theta_levels,       &
                    theta_conv, q_conv, qcl_conv, qcf_conv,              &
                    cf_liquid_conv, cf_frozen_conv, bulk_cf_conv)

        END IF  ! L_calc_dxek_if0
        !---------------------------------------------------------------------
        ! Check for negative q being passed to convection
        ! Note above PC2 code checking no negative qcl and qcf
        ! Nothing currently preventing negative q going to convection
        !---------------------------------------------------------------------
        DO k=1,tdims%k_end
          DO j=1,rows
            DO i=1,row_length
              IF (q_conv(i,j,k) < qmin_conv) THEN
                IF (printstatus >= prstatus_normal) THEN
                  WRITE(umMessage,                                             &
                   '(a,3(i0,1x),a,g26.18,a,g26.18,a,g26.18,a,g26.18)')         &
                   ' Negative q at start of convection: i,j,k ',i,j,k,         &
                   ' q_conv ', q_conv(i,j,k),' q_n ',q_n(i,j,k),               &
                   ' dqbydt ', dqbydt(i,j,k),                                  &
                   ' inc ',fraction_step*q_inc_step(i,j,k)
                  CALL umPrint(umMessage,src='ni_conv_ctl')
                END IF
                ! store added q
                dq_add(i,j,k) = dq_add(i,j,k) +(qmin_conv - q_conv(i,j,k))
                ! Reset to qmin_conv in a non-conservative way
                q_conv(i,j,k) = qmin_conv
              END IF
            END DO ! i
          END DO ! j
        END DO ! k
        !---------------------------------------------------------------------

        IF (model_type == mt_single_column) THEN
          ! The SCM output diagnostic system needs to know what sub-step
          ! Note this method of fixing SCM substepping problems for SCM output
          ! with l_rediagnosis = .true. is probably not ideal but works.

          ! DEPENDS ON: scm_substep_start
          CALL scm_substep_start(call_number)
        END IF ! model_type

        SELECT CASE ( i_convection_vn )
        CASE ( i_convection_vn_5a )

          CALL conv_diag_5a(                                        &

          ! IN Parallel variables
           row_length, rows                                         &

          ! IN model dimensions.
          , bl_levels                                               &
          ,p, p_layer_centres(1,1,1)                                &
          ,exner_rho_levels                                         &
          ,rho_only, rho_theta, z_theta, z_rho                      &

          ! IN Model switches
          ,l_extra_call                                             &
          ,no_cumulus                                               &
          ! IN cloud data
          ,qcf_conv, qcl_conv, bulk_cf_conv                         &

          ! IN everything not covered so far :
          ,p_star, q_conv, theta_conv, exner_layer_centres(:,:,1:)  &
          ,u_conv, v_conv, u_0_p, v_0_p                             &
          ,tstar_land, tstar_sea, tstar_sice, z0msea                &
          ,flux_e, flux_h, ustar_in, l_spec_z0, z0m_scm, z0h_scm    &
          , t_surf, land_sea_mask, flandg, ice_fract, timestep      &
          , w, w_max, deep_flag, past_precip, past_conv_ht          &

          ! IN surface fluxes
          , fb_surf, ustar                                          &

          ! SCM Diagnostics (dummy values in full UM)
          , nscmdpkgs,l_scmdiags                                    &

          ! INOUT data required elsewhere in UM system :
          ,zh_copy,ztop,dzh,qcl_inv_top,zlcl,zlcl_uv,delthvu,ql_ad  &
          ,ntml,ntpar,nlcl,cumulus                                  &
          ,l_shallow,l_congestus,l_congestus2, conv_type            &
          ,cin_undilute,cape_undilute, wstar, wthvs                 &
          ,entrain_coef, qsat_lcl                                   &
          ,error_code                                               &
             )
        CASE ( i_convection_vn_6a )

          CALL conv_diag_6a(                                        &

          ! IN Parallel variables
           row_length, rows                                         &

          ! IN model dimensions.
          , bl_levels                                               &
          ,p, p_layer_centres(1,1,1)                                &
          ,exner_rho_levels                                         &
          ,rho_only, rho_theta, z_theta, z_rho                      &

          ! IN Model switches
          ,l_extra_call                                             &
          ,no_cumulus                                               &
          ! IN cloud data
          ,qcf_conv, qcl_conv, bulk_cf_conv                         &

          ! IN everything not covered so far :
          ,p_star, q_conv, theta_conv, exner_layer_centres(:,:,1:)  &
          ,u_conv, v_conv, u_0_p, v_0_p                             &
          ,tstar_land, tstar_sea, tstar_sice, z0msea                &
          ,flux_e, flux_h, ustar_in, l_spec_z0, z0m_scm, z0h_scm    &
          , t_surf, land_sea_mask, flandg, ice_fract                &
          , w, w_max, deep_flag, past_precip, past_conv_ht          &
          , conv_prog_precip                                        &
          , g_ccp, h_ccp                                            &

          ! IN surface fluxes
          , fb_surf, ustar                                          &

          ! SCM Diagnostics (dummy values in full UM)
          , nscmdpkgs,l_scmdiags                                    &

          ! INOUT data required elsewhere in UM system :
          ,zh_copy,ztop,dzh,qcl_inv_top,zlcl,zlcl_uv,delthvu,ql_ad  &
          ,ntml,ntpar,nlcl,cumulus                                  &
          ,l_shallow,l_congestus,l_congestus2, conv_type            &
          ,cin_undilute,cape_undilute, wstar, wthvs                 &
          ,entrain_coef, qsat_lcl                                   &
          ,error_code                                               &
             )

        CASE DEFAULT ! i_convection_vn

          errorstatus = 10
          WRITE (cmessage,'(A,A,I6)') 'Convection scheme version value',  &
                                  'i_convection_vn = ',i_convection_vn
          CALL Ereport ( routinename, errorstatus, cmessage)

        END SELECT ! i_convection_vn

        ! Reset l_mid(i,j) = .true. as mid can occur anywhere on
        ! rediagnosis
        DO j = 1, rows
          DO i = 1, row_length
            l_mid(i,j) = .TRUE.
          END DO
        END DO

      END IF  ! Call_number > 1
    END IF   ! test on l_rediagnosis

    !---------------------------------------------------------------------------
    ! 2.4 WORK OUT NUMBER OF CUMULUS POINTS IN EACH CONVECTION SEGMENT
    !---------------------------------------------------------------------------
    ! Also calculate number of shallow and deep points in each convection
    ! segment

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( rows, row_length, cumulus_1d, cumulus, l_shallow_1d,     &
!$OMP         l_shallow, l_congestus_1d, l_congestus, l_mid_1d, l_mid, &
!$OMP         l_no_cumulus_1d, no_cumulus )                            &
!$OMP PRIVATE( i, j, ii )
!DIR$ IVDEP
    DO j = 1, rows
!DIR$ IVDEP
      DO i = 1, row_length
        ii = row_length*(j-1) + i
        cumulus_1d(ii) = cumulus(i,j)
        l_shallow_1d(ii) = l_shallow(i,j)
        l_congestus_1d(ii) = l_congestus(i,j)
        l_mid_1d(ii) = l_mid(i,j)
        l_no_cumulus_1d(ii) = no_cumulus(i,j)
      END DO
    END DO
!$OMP END PARALLEL DO

    !===============================================================================
    ! Start of segmented region
    !===============================================================================

    !Number of threads is 1 at this point: we are not load-balancing between
    !threads, then segmenting here.  Rather, set up segmentation on one thread,
    !then distribute the segments among threads.
    num_parallel=1
    ipar=1

    !Set up automatic segment tuning.
    IF (l_autotune_segments .AND. .NOT. ALLOCATED(autotune_state)) THEN
      ALLOCATE(autotune_state)
      CALL autotune_init(                     &
        autotune_state,                       &
        name           = 'Convection',        &
        tag            = 'CONV',              &
        start_size     = a_convect_seg_size,  &
        calls_per_step = numcycles*n_conv_calls)
    END IF

    !If autotuning is active, retrieves a new trial segment size. If not, the
    !returned segment size will stay the same.
    segment_size = a_convect_seg_size
    IF (l_autotune_segments) THEN
      segment_size = autotune_get_trial_size(autotune_state)
    END IF

    !Set up segment meta-information.
    CALL segments_mod_seg_meta(meta_segments, ipar, num_parallel,      &
      rows*row_length, segment_size, a_convect_segments)

    ! Allocate space for segmentation arrays
    ALLOCATE(  segments( meta_segments%num_segments ) )

    ! Work out starting points and sizes of each segment individually
    CALL segments_mod_segments(segments, meta_segments, row_length, rows)

    IF (l_autotune_segments) THEN
      CALL autotune_start_region(autotune_state)
    END IF

!$OMP  PARALLEL  DEFAULT(SHARED)                             &
!$OMP& PRIVATE(ii, jj, i, j, errorstatus, n_cumulus, n_deep, &
!$OMP&         n_shallow, n_congestus, n_mid)

!$OMP DO SCHEDULE(DYNAMIC)
    DO i=1, meta_segments%num_segments

      n_cumulus   = 0
      n_deep      = 0
      n_shallow   = 0
      n_congestus = 0
      n_mid       = 0

      DO j=segments(i)%fp, segments(i)%fp+segments(i)%seg_points-1
        IF (cumulus_1d(j)) THEN
          n_cumulus=n_cumulus+1
          IF (iconv_deep >  0 .AND. .NOT. l_shallow_1d(j) .AND.      &
                      .NOT. l_congestus_1d(j)) THEN
            n_deep = n_deep+1
          END IF
          IF (iconv_shallow >  0 .AND. l_shallow_1d(j) .AND.         &
               .NOT. l_congestus_1d(j)) THEN
            n_shallow = n_shallow+1
          END IF
          IF (iconv_congestus >  0) THEN
            IF (l_congestus_1d(j)) THEN
              n_congestus = n_congestus+1
            END IF
          ELSE
            n_congestus = 1     ! may be required for dim
          END IF
        END IF
      END DO

      DO j=segments(i)%fp, segments(i)%fp+segments(i)%seg_points-1
        IF (l_mid_1d(j)) THEN
          n_mid=n_mid+1
        END IF
      END DO

    !---------------------------------------------------------------------------
    !  2.5 Loop over convection segments calling convection
    !---------------------------------------------------------------------------

      ii = segments(i)%first_x
      jj = segments(i)%first_y

      SELECT CASE ( i_convection_vn )

      CASE ( i_convection_vn_5a )

        CALL glue_conv_5a                                                 &
          ( rows*row_length, segments(i)%seg_points                       &
          , n_conv_levels, bl_levels                                      &
          , call_number, i                                                &
          , theta_conv(ii,jj,1), q_conv(ii,jj,1)                          &
          , qcl_conv(ii,jj,1), qcf_conv(ii,jj,1)                          &
          , qrain_conv(ii,jj,1), qgraup_conv(ii,jj,1), qcf2_conv(ii,jj,1) &
          , cf_liquid_conv(ii,jj,1), cf_frozen_conv(ii,jj,1)              &
          , bulk_cf_conv(ii,jj,1)                                         &
          , p_star(ii,jj), land_sea_mask(ii,jj)                           &
          , u_conv(ii,jj,1), v_conv(ii,jj,1), w(ii,jj,1)                  &
          , tot_tracer(ii,jj,1,1), dthbydt(ii,jj,1)                       &
          , dqbydt(ii,jj,1),   dqclbydt(ii,jj,1), dqcfbydt(ii,jj,1)       &
          , dcflbydt(ii,jj,1), dcffbydt(ii,jj,1), dbcfbydt(ii,jj,1)       &
          , dubydt_p(ii,jj,1), dvbydt_p(ii,jj,1)                          &
          , it_conv_rain(ii,jj),      it_conv_snow(ii,jj)                 &
          , it_conv_rain_3d(ii,jj,1), it_conv_snow_3d(ii,jj,1)            &
          , it_cca0_dp(ii,jj,1), it_cca0_md(ii,jj,1), it_cca0_sh(ii,jj,1) &
          , it_cca0(ii,jj,1),  it_ccb0(ii,jj),   it_cct0(ii,jj)           &
          , it_cclwp0(ii,jj),  it_ccw0(ii,jj,1), it_lcbase0(ii,jj)        &
          , it_cca0_2d(ii,jj), it_lctop(ii,jj),  it_lcca(ii,jj)           &
          , it_cca(ii,jj,1),   it_ccb(ii,jj),    it_cct(ii,jj)            &
          , it_cclwp(ii,jj),   it_ccw(ii,jj,1),  it_lcbase(ii,jj)         &
          , it_cca_2d(ii,jj), freeze_lev(ii,jj), it_dp_cfl_limited(ii,jj) &
          , it_md_cfl_limited(ii,jj)                                      &
          , it_mid_level(ii,jj), it_kterm_deep(ii,jj)                     &
          , it_kterm_shall(ii,jj)                                         &
          , it_precip_dp(ii,jj), it_precip_sh(ii,jj)                      &
          , it_precip_md(ii,jj), it_precip_cg(ii,jj)                      &
          , it_wstar_dn(ii,jj),  it_wstar_up(ii,jj)                       &
          , it_mb1(ii,jj), it_mb2(ii,jj), it_cg_term(ii,jj), n_cumulus    &
          , uw0_p(ii,jj), vw0_p(ii,jj), w_max(ii,jj)                      &
          , zlcl(ii,jj), zlcl_uv(ii,jj), ztop(ii,jj), entrain_coef(ii,jj) &
          , deep_flag(ii,jj), past_precip(ii,jj), past_conv_ht(ii,jj)     &
          , it_cape_out(ii,jj), n_deep, n_congestus, n_shallow            &
          , n_mid, r_rho_lev(ii,jj,1), r_theta_lev(ii,jj,0)               &
          , rho_only(ii,jj,1), rho_theta(ii,jj,1), rho_dry(ii,jj,1)       &
          , rho_dry_theta(ii,jj,1), delta_smag(ii,jj)                     &
          , exner_layer_boundaries(ii,jj,0), exner_layer_centres(ii,jj,0) &
          , p_layer_boundaries(ii,jj,0), p_layer_centres(ii,jj,0)         &
          , z_theta(ii,jj,1), z_rho(ii,jj,1), timestep_conv               &
          , t1_sd(ii,jj), q1_sd(ii,jj), ntml(ii,jj), ntpar(ii,jj)         &
          , conv_type(ii,jj), l_shallow(ii,jj), l_pc2_diag_sh_pts(ii,jj)  &
          , l_congestus(ii,jj), l_mid(ii,jj), cumulus(ii,jj)              &
          , wstar(ii,jj), wthvs(ii,jj), delthvu(ii,jj), ql_ad(ii,jj)      &
          , qsat_lcl(ii,jj), ftl(ii,jj), fqt(ii,jj), l_tracer, ntra_fld   &
          , ntra_lev, n_cca_lev, l_mcr_qrain                              &
          , l_mcr_qgraup, l_mcr_qcf2, l_calc_dxek                         &
          , l_q_interact, it_up_flux_half(ii,jj,1)                        &
          , it_up_flux(ii,jj,1),      it_dwn_flux(ii,jj,1)                &
          , it_entrain_up(ii,jj,1),   it_detrain_up(ii,jj,1)              &
          , it_entrain_dwn(ii,jj,1),  it_detrain_dwn(ii,jj,1)             &
          , it_uw_dp(ii,jj,1),        it_vw_dp(ii,jj,1)                   &
          , it_uw_shall(ii,jj,1),     it_vw_shall(ii,jj,1)                &
          , it_uw_mid(ii,jj,1),       it_vw_mid(ii,jj,1)                  &
          , it_wqt_flux(ii,jj,1),     it_wthetal_flux(ii,jj,1)            &
          , it_wthetav_flux(ii,jj,1), it_wql_flux(ii,jj,1)                &
          , it_mf_deep(ii,jj,1),      it_mf_congest(ii,jj,1)              &
          , it_mf_shall(ii,jj,1),     it_mf_midlev(ii,jj,1)               &
          , it_dt_deep(ii,jj,1),      it_dt_congest(ii,jj,1)              &
          , it_dt_shall(ii,jj,1),     it_dt_midlev(ii,jj,1)               &
          , it_dq_deep(ii,jj,1),      it_dq_congest(ii,jj,1)              &
          , it_dq_shall(ii,jj,1),     it_dq_midlev(ii,jj,1)               &
          , it_du_deep(ii,jj,1),      it_du_congest(ii,jj,1)              &
          , it_du_shall(ii,jj,1),     it_du_midlev(ii,jj,1)               &
          , it_dv_deep(ii,jj,1),      it_dv_congest(ii,jj,1)              &
          , it_dv_shall(ii,jj,1),     it_dv_midlev(ii,jj,1)               &
          , ind_cape_reduced(ii,jj),  cape_ts_used(ii,jj)                 &
          , it_ind_deep(ii,jj),       it_ind_shall(ii,jj)                 &
          , it_w2p(ii,jj,1), it_dt_dd(ii,jj,1), it_dq_dd(ii,jj,1)         &
          , it_du_dd(ii,jj,1), it_dv_dd(ii,jj,1), it_area_ud(ii,jj,1)     &
          , it_area_dd(ii,jj,1)                                           & 
                  )

      CASE ( i_convection_vn_6a )

        CALL glue_conv_6a                                                 &
          ( rows*row_length, segments(i)%seg_points                       &
          , n_conv_levels, bl_levels                                      &
          , call_number, i                                                &
          , theta_conv(ii,jj,1), q_conv(ii,jj,1)                          &
          , qcl_conv(ii,jj,1), qcf_conv(ii,jj,1)                          &
          , qrain_conv(ii,jj,1), qgraup_conv(ii,jj,1), qcf2_conv(ii,jj,1) &
          , cf_liquid_conv(ii,jj,1), cf_frozen_conv(ii,jj,1)              &
          , bulk_cf_conv(ii,jj,1)                                         &
          , p_star(ii,jj), land_sea_mask(ii,jj)                           &
          , u_conv(ii,jj,1), v_conv(ii,jj,1), w(ii,jj,1)                  &
          , tot_tracer(ii,jj,1,1), dthbydt(ii,jj,1)                       &
          , dqbydt(ii,jj,1),   dqclbydt(ii,jj,1), dqcfbydt(ii,jj,1)       &
          , dcflbydt(ii,jj,1), dcffbydt(ii,jj,1), dbcfbydt(ii,jj,1)       &
          , dubydt_p(ii,jj,1), dvbydt_p(ii,jj,1)                          &
          , it_conv_rain(ii,jj),      it_conv_snow(ii,jj)                 &
          , it_conv_rain_3d(ii,jj,1), it_conv_snow_3d(ii,jj,1)            &
          , it_cca0_dp(ii,jj,1), it_cca0_md(ii,jj,1), it_cca0_sh(ii,jj,1) &
          , it_cca0(ii,jj,1),  it_ccb0(ii,jj),   it_cct0(ii,jj)           &
          , it_cclwp0(ii,jj),  it_ccw0(ii,jj,1), it_lcbase0(ii,jj)        &
          , it_cca0_2d(ii,jj), it_lctop(ii,jj),  it_lcca(ii,jj)           &
          , it_cca(ii,jj,1),   it_ccb(ii,jj),    it_cct(ii,jj)            &
          , it_cclwp(ii,jj),   it_ccw(ii,jj,1),  it_lcbase(ii,jj)         &
          , it_cca_2d(ii,jj), freeze_lev(ii,jj), it_dp_cfl_limited(ii,jj) &
          , it_md_cfl_limited(ii,jj)                                      &
          , it_mid_level(ii,jj), it_kterm_deep(ii,jj)                     &
          , it_kterm_shall(ii,jj)                                         &
          , it_precip_dp(ii,jj), it_precip_sh(ii,jj)                      &
          , it_precip_md(ii,jj), it_precip_cg(ii,jj)                      &
          , it_wstar_dn(ii,jj),  it_wstar_up(ii,jj)                       &
          , it_mb1(ii,jj), it_mb2(ii,jj), it_cg_term(ii,jj), n_cumulus    &
          , uw0_p(ii,jj), vw0_p(ii,jj), w_max(ii,jj)                      &
          , zlcl(ii,jj), zlcl_uv(ii,jj), ztop(ii,jj), entrain_coef(ii,jj) &
          , conv_prog_precip_c(ii,jj,1)                                   &
          , deep_flag(ii,jj), past_precip(ii,jj), past_conv_ht(ii,jj)     &
          , it_cape_out(ii,jj), n_deep, n_congestus, n_shallow            &
          , n_mid, r_rho_lev(ii,jj,1), r_theta_lev(ii,jj,0)               &
          , rho_only(ii,jj,1), rho_theta(ii,jj,1), rho_dry(ii,jj,1)       &
          , rho_dry_theta(ii,jj,1), delta_smag(ii,jj)                     &
          , exner_rho_lev(ii,jj,1)                                        &
          , exner_layer_boundaries(ii,jj,0), exner_layer_centres(ii,jj,0) &
          , p_layer_boundaries(ii,jj,0), p_layer_centres(ii,jj,0)         &
          , z_theta(ii,jj,1), z_rho(ii,jj,1), timestep_conv               &
          , t1_sd(ii,jj), q1_sd(ii,jj), ntml(ii,jj), ntpar(ii,jj)         &
          , conv_type(ii,jj), l_shallow(ii,jj), l_pc2_diag_sh_pts(ii,jj)  &
          , l_congestus(ii,jj), l_mid(ii,jj), cumulus(ii,jj)              &
          , wstar(ii,jj), wthvs(ii,jj), delthvu(ii,jj), ql_ad(ii,jj)      &
          , qsat_lcl(ii,jj), ftl(ii,jj), fqt(ii,jj), l_tracer, ntra_fld   &
          , ntra_lev, n_cca_lev, l_mcr_qrain                              &
          , l_mcr_qgraup, l_mcr_qcf2, l_calc_dxek                         &
          , l_q_interact, it_up_flux_half(ii,jj,1)                        &
          , it_up_flux(ii,jj,1),      it_dwn_flux(ii,jj,1)                &
          , it_entrain_up(ii,jj,1),   it_detrain_up(ii,jj,1)              &
          , it_entrain_dwn(ii,jj,1),  it_detrain_dwn(ii,jj,1)             &
          , it_uw_dp(ii,jj,1),        it_vw_dp(ii,jj,1)                   &
          , it_uw_shall(ii,jj,1),     it_vw_shall(ii,jj,1)                &
          , it_uw_mid(ii,jj,1),       it_vw_mid(ii,jj,1)                  &
          , it_wqt_flux(ii,jj,1),     it_wthetal_flux(ii,jj,1)            &
          , it_wthetav_flux(ii,jj,1), it_wql_flux(ii,jj,1)                &
          , it_mf_deep(ii,jj,1),      it_mf_congest(ii,jj,1)              &
          , it_mf_shall(ii,jj,1),     it_mf_midlev(ii,jj,1)               &
          , it_dt_deep(ii,jj,1),      it_dt_congest(ii,jj,1)              &
          , it_dt_shall(ii,jj,1),     it_dt_midlev(ii,jj,1)               &
          , it_dq_deep(ii,jj,1),      it_dq_congest(ii,jj,1)              &
          , it_dq_shall(ii,jj,1),     it_dq_midlev(ii,jj,1)               &
          , it_du_deep(ii,jj,1),      it_du_congest(ii,jj,1)              &
          , it_du_shall(ii,jj,1),     it_du_midlev(ii,jj,1)               &
          , it_dv_deep(ii,jj,1),      it_dv_congest(ii,jj,1)              &
          , it_dv_shall(ii,jj,1),     it_dv_midlev(ii,jj,1)               &
          , ind_cape_reduced(ii,jj),  cape_ts_used(ii,jj)                 &
          , it_ind_deep(ii,jj),       it_ind_shall(ii,jj)                 &
          , it_w2p(ii,jj,1), it_dt_dd(ii,jj,1), it_dq_dd(ii,jj,1)         &
          , it_du_dd(ii,jj,1), it_dv_dd(ii,jj,1), it_area_ud(ii,jj,1)     &
          , it_area_dd(ii,jj,1)                                           &
          , scm_convss_dg(ii,jj), l_scm_convss_dg                         &
           )

      CASE DEFAULT ! i_convection_vn

        errorstatus = 10
        WRITE (cmessage,'(A,A,I6)') 'Convection scheme version value',     &
                  ' not recognised. i_convection_vn = ',i_convection_vn
        CALL Ereport ( routinename, errorstatus, cmessage)

      END SELECT ! i_convection_vn

    END DO  ! loop over number of segments
!$OMP END DO nowait
!$OMP END PARALLEL

    !If autotuning is active, decide what to do with the
    !segment size and report the current status.
    IF (l_autotune_segments) THEN
      CALL autotune_stop_region(autotune_state, row_length*rows)
      CALL autotune_advance(autotune_state)
      CALL autotune_report(autotune_state,      &
          quiet=(call_number < n_conv_calls .OR. cycleno < numcycles))
    END IF

    ! Deallocate segmentation variables
    DEALLOCATE( segments  )

    !===============================================================================
    ! End of segmented region
    !===============================================================================

    ! Update convective cloud diagnostics section 5
    ! and prognostics in section 0  after sub-step

    CALL update_conv_cloud (rows, row_length,                           &
       n_cca_lev, ndecay_steps, n_conv_calls,                           &
       call_number,                                                     &
       it_ccb, it_cct, it_lcbase, it_lctop,it_ccb0, it_cct0, it_lcbase0,&
       one_over_conv_calls, timestep,decay_time,                        &
       p_layer_boundaries, cld_life_3d,                                 &
       it_cca_2d, it_lcca, it_cca, it_ccw, it_cclwp,                    &
       it_cca0_2d, it_cca0,it_ccw0,                                     &
       it_cca0_dp, it_cca0_md, it_cca0_sh,                              &
       ! in/out
       ccb, cct, lcbase,lctop, ccb0, cct0, lcbase0,                     &
       ccb0_local, cct0_local,lcbase0_local,                            &
       lcca,  cca_2d,   cclwp,  cca,  ccw,                              &
       cca0_2d,  cclwp0, cca0, ccw0,                                    &
       cca0_dp, cca0_md, cca0_sh,                                       &
       cca0_local_dp, cca0_local_md, cca0_local_sh,                     &
       it_cclwp0, cclwp0_local, cca0_local, ccw0_local)

    ! Update other diagnostics after sub-step

    CALL update_conv_diags(rows, row_length,                              &
      n_conv_levels, call_number, n_conv_calls,                           &
      ntml,ntpar, freeze_lev, it_kterm_deep, it_kterm_shall, it_cg_term,  &
      cumulus, l_shallow, l_congestus, l_congestus2,it_dp_cfl_limited,    &
      it_md_cfl_limited, it_mid_level,                                    &
      one_over_conv_calls, timestep_conv , wstar,                         &
      exner_theta_levels, z_rho,                                          &
      ! Input (values output by latest convection call)
      it_cape_out, it_conv_rain, it_conv_snow,                            &
      it_precip_dp, it_precip_sh, it_precip_md, it_precip_cg,             &
      ind_cape_reduced, cape_ts_used, it_ind_deep, it_ind_shall,          &
      it_wstar_up, it_mb1, it_mb2,                                        &
      dubydt_p, dvbydt_p,                                                 &
      it_up_flux, it_up_flux_half, it_dwn_flux, it_entrain_up,            &
      it_entrain_dwn, it_detrain_up, it_detrain_dwn,                      &
      it_conv_rain_3d, it_conv_snow_3d,                                   &
      it_wqt_flux,it_wthetal_flux,it_wthetav_flux,it_wql_flux,            &
      it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall, it_uw_mid, it_vw_mid, &
      it_mf_deep, it_mf_congest, it_mf_shall, it_mf_midlev,               &
      it_dt_deep, it_dt_congest, it_dt_shall, it_dt_midlev,               &
      it_dq_deep, it_dq_congest, it_dq_shall, it_dq_midlev,               &
      it_du_deep, it_du_congest, it_du_shall, it_du_midlev,               &
      it_dv_deep, it_dv_congest, it_dv_shall, it_dv_midlev,               &
      it_w2p,     it_dt_dd, it_dq_dd, it_du_dd, it_dv_dd,                 &
      it_area_ud, it_area_dd,                                             &
      ! in/out diagnostics
      dubydt_pout, dvbydt_pout, conv_rain, conv_snow)


    ! Output SCM convection diagnostics that are done separately
    ! each substep, rather than output as a mean over the substeps
    IF (l_scm_convss_dg) THEN
      CALL scm_convss_dg_output( scm_convss_dg, call_number )
    END IF


    ! Add tendency increments to arrays in physics_tendencies_mod
    IF (l_retain_conv_tendencies) THEN
      ! Add convection increments for Temperature
!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(k, i, j)                        &
!$OMP& SHARED(n_conv_levels,tdims,dt_conv,dthbydt,timestep_conv,     &
!$OMP&     exner_theta_levels,dq_conv,dqbydt)

!$OMP DO SCHEDULE(STATIC)  
      DO k=1,n_conv_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            dt_conv(i,j,k )= dt_conv(i,j,k) +                             &
                               dthbydt(i,j,k) * timestep_conv             &
                               * exner_theta_levels(i,j,k)
          END DO
        END DO
      END DO      
!$OMP END DO NOWAIT

      ! Add convection increments for q
!$OMP DO SCHEDULE(STATIC)  
      DO k=1,n_conv_levels
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            dq_conv(i,j,k) = dq_conv(i,j,k) +                             &
                               dqbydt(i,j,k)*timestep_conv
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL
    END IF
    !-----------------------------------------------------------
    ! Mid-level conv is only possible on subsequent convective
    ! substeps if convection has occurred on the previous
    ! convective substep unless using rediagnosis option.
    !-----------------------------------------------------------
    IF (.NOT. l_rediagnosis) THEN
      DO j = 1, rows
        DO i= 1, row_length
          l_mid(i,j) = it_mid_level(i,j) .OR. cumulus(i,j)
        END DO
      END DO
    END IF


    ! ----------------------------------------------------------------------
    ! Section 2.8 Add on theta and q increments, qcl and qcf increments.
    ! ----------------------------------------------------------------------

        ! add on increments to theta and q for next convection call
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
    DO k = 1, n_conv_levels
      DO j = 1, rows
        DO i = 1, row_length
          theta_inc(i,j,k) = theta_inc(i,j,k)                     &
                           + dthbydt(i,j,k) * timestep_conv
          q_inc(i,j,k) = q_inc(i,j,k)                             &
                          + dqbydt(i,j,k) * timestep_conv
          qcl_inc(i,j,k) = qcl_inc(i,j,k) +                       &
                         (dqclbydt(i,j,k) * timestep_conv)
          qcf_inc(i,j,k) = qcf_inc(i,j,k) +                       &
                         (dqcfbydt(i,j,k) * timestep_conv)
          cf_liquid_inc(i,j,k) = cf_liquid_inc(i,j,k) +           &
                         (dcflbydt(i,j,k) * timestep_conv)
          cf_frozen_inc(i,j,k) = cf_frozen_inc(i,j,k) +           &
                         (dcffbydt(i,j,k) * timestep_conv)
          bulk_cf_inc(i,j,k)   = bulk_cf_inc(i,j,k) +             &
                         (dbcfbydt(i,j,k) * timestep_conv)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

        ! Final sweep  just add on convection increments from last sub-step

    IF (call_number == n_conv_calls) THEN

      IF (l_rediagnosis) THEN
        ! No longer need total increment arrays
        IF (l_mcr_qcf2) THEN
          DEALLOCATE( qcf2_inc_step )
        END IF
        IF (l_mcr_qgraup) THEN
          DEALLOCATE( qgraup_inc_step )
        END IF
        IF (l_mcr_qrain) THEN
          DEALLOCATE( qrain_inc_step )
        END IF
        DEALLOCATE( bulk_cf_inc_step )
        DEALLOCATE( cf_frozen_inc_step )
        DEALLOCATE( cf_liquid_inc_step )
        DEALLOCATE( qcf_inc_step )
        DEALLOCATE( qcl_inc_step )
        DEALLOCATE( q_inc_step )
        DEALLOCATE( theta_inc_step )
      END IF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k, i, j)
      IF (l_safe_conv) THEN
        ! remove any q added for a safe qmin input profile
!$OMP DO SCHEDULE(STATIC)
        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length
              q_conv(i,j,k)   = q_conv(i,j,k) - dq_add(i,j,k)              &
                                  + dqbydt(i,j,k) * timestep_conv
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO SCHEDULE(STATIC)
        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length
              q_conv(i,j,k)   = q_conv(i,j,k)                              &
                                     + dqbydt(i,j,k) * timestep_conv
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF

!$OMP DO SCHEDULE(STATIC)
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length
            theta_conv(i,j,k) = theta_conv(i,j,k)                        &
                          + dthbydt(i,j,k) * timestep_conv
            qcl_conv(i,j,k) = qcl_conv(i,j,k)                            &
                            + (dqclbydt(i,j,k) * timestep_conv)
            qcf_conv(i,j,k) = qcf_conv(i,j,k)                            &
                            + (dqcfbydt(i,j,k) * timestep_conv)
            cf_liquid_conv(i,j,k) = cf_liquid_conv(i,j,k)                &
                            + (dcflbydt(i,j,k) * timestep_conv)
            cf_frozen_conv(i,j,k) = cf_frozen_conv(i,j,k)                &
                            + (dcffbydt(i,j,k) * timestep_conv)
            bulk_cf_conv(i,j,k)  = bulk_cf_conv(i,j,k)                   &
                            + (dbcfbydt(i,j,k) * timestep_conv)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

    ELSE

      IF (l_rediagnosis) THEN

        ! not final sub-step so add on part of increment from
        ! slow physics and dynamics

        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length

              theta_conv(i,j,k) = theta_conv(i,j,k)                     &
                           + fraction_step*theta_inc_step(i,j,k)        &
                           + dthbydt(i,j,k) * timestep_conv
              q_conv(i,j,k) = q_conv(i,j,k)                             &
                            + fraction_step*q_inc_step(i,j,k)           &
                            + dqbydt(i,j,k) * timestep_conv
              qcl_conv(i,j,k) = qcl_conv(i,j,k)                         &
                           + fraction_step*qcl_inc_step(i,j,k)          &
                           +(dqclbydt(i,j,k) * timestep_conv)
              qcf_conv(i,j,k) = qcf_conv(i,j,k)                         &
                           + fraction_step*qcf_inc_step(i,j,k)          &
                           +(dqcfbydt(i,j,k) * timestep_conv)
              cf_liquid_conv(i,j,k) = cf_liquid_conv(i,j,k)             &
                           + fraction_step*cf_liquid_inc_step(i,j,k)    &
                           +(dcflbydt(i,j,k) * timestep_conv)
              cf_frozen_conv(i,j,k) = cf_frozen_conv(i,j,k)             &
                           + fraction_step*cf_frozen_inc_step(i,j,k)    &
                           + (dcffbydt(i,j,k) * timestep_conv)
              bulk_cf_conv(i,j,k)   = bulk_cf_conv(i,j,k)               &
                           + fraction_step*bulk_cf_inc_step(i,j,k)      &
                           +(dbcfbydt(i,j,k) * timestep_conv)

            END DO
          END DO
        END DO

        ! Extra prognostics if present
        IF (l_mcr_qrain) THEN
          DO k=1,tdims%k_end
            DO j=1,rows
              DO i=1,row_length
                qrain_conv(i,j,k)   = qrain_conv(i,j,k)                  &
                                     + fraction_step*qrain_inc_step(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
        END IF
        IF (l_mcr_qgraup) THEN
          DO k=1,tdims%k_end
            DO j=1,rows
              DO i=1,row_length
                qgraup_conv(i,j,k)   = qgraup_conv(i,j,k)                &
                                      + fraction_step*qgraup_inc_step(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
        END IF
        IF (l_mcr_qcf2) THEN
          DO k=1,tdims%k_end
            DO j=1,rows
              DO i=1,row_length
                qcf2_conv(i,j,k)   = qcf2_conv(i,j,k)                    &
                                    + fraction_step*qcf2_inc_step(i,j,k)
              END DO ! i
            END DO ! j
          END DO ! k
        END IF

        IF (l_mom) THEN
          DO k = 1, n_conv_levels
            DO j = 1, rows
              DO i = 1, row_length
                u_conv(i,j,k) = u_conv(i,j,k)                           &
                              + fraction_step*u_inc_step(i,j,k)         &
                              + dubydt_p(i,j,k) * timestep_conv
                v_conv(i,j,k) = v_conv(i,j,k)                           &
                              + fraction_step*v_inc_step(i,j,k)         &
                              + dvbydt_p(i,j,k) * timestep_conv

              END DO
            END DO
          END DO
        END IF   ! l_mom

      ELSE    ! original code
        ! Note no need to update extra prognostics if present

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

!$OMP DO SCHEDULE(STATIC)
        DO k = 1, n_conv_levels
          DO j = 1, rows
            DO i = 1, row_length
              theta_conv(i,j,k) = theta_conv(i,j,k)                     &
                           + dthbydt(i,j,k) * timestep_conv
              q_conv(i,j,k) = q_conv(i,j,k)                             &
                           + dqbydt(i,j,k) * timestep_conv
              qcl_conv(i,j,k) = qcl_conv(i,j,k)                         &
                           +(dqclbydt(i,j,k) * timestep_conv)
              qcf_conv(i,j,k) = qcf_conv(i,j,k)                         &
                           +(dqcfbydt(i,j,k) * timestep_conv)
              cf_liquid_conv(i,j,k) = cf_liquid_conv(i,j,k)             &
                           +(dcflbydt(i,j,k) * timestep_conv)
              cf_frozen_conv(i,j,k) = cf_frozen_conv(i,j,k)             &
                           + (dcffbydt(i,j,k) * timestep_conv)
              bulk_cf_conv(i,j,k)   = bulk_cf_conv(i,j,k)               &
                           +(dbcfbydt(i,j,k) * timestep_conv)
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT

        IF (l_mom) THEN

!$OMP DO SCHEDULE(STATIC)
          DO k = 1, n_conv_levels
            DO j = 1, rows
              DO i = 1, row_length
                u_conv(i,j,k) = u_conv(i,j,k)                           &
                             + dubydt_p(i,j,k) * timestep_conv
                v_conv(i,j,k) = v_conv(i,j,k)                           &
                             + dvbydt_p(i,j,k) * timestep_conv

              END DO
            END DO
          END DO
!$OMP END DO NOWAIT

        END IF    !l_mom

!$OMP END PARALLEL

      END IF    ! l_rediagnosis

    END IF    ! step number


    IF (i_pc2_conv_coupling /= pc2_conv_original) THEN
      ! Protect against generation of inconsistently low cloud
      ! fraction implying very high in-cloud condensate amounts.
      ! In-cloud condensate amounts above 2.0e-3 lead to
      ! cloud fraction being increased (up to a value of 1.0)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, tmpB4)
!$OMP DO SCHEDULE(STATIC)
      DO k = 1, n_conv_levels
        DO j = 1, rows
          DO i = 1, row_length

            ! Liquid cloud fraction
            IF (cf_liquid_conv(i,j,k) > 0.0) THEN
              IF ((qcl_conv(i,j,k)/cf_liquid_conv(i,j,k))>2.0e-3) THEN
                tmpB4 = cf_liquid_conv(i,j,k)
                cf_liquid_conv(i,j,k) = MIN(1.0,qcl_conv(i,j,k)/2.0e-3)
                cf_liquid_inc(i,j,k)  = cf_liquid_inc(i,j,k)          &
                  + cf_liquid_conv(i,j,k) - tmpB4
                IF (l_fixbug_pc2_mixph) THEN
                  bulk_cf_conv(i,j,k) = bulk_cf_conv(i,j,k)           &
                  + cf_liquid_conv(i,j,k) - tmpB4
                  bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)             &
                  + cf_liquid_conv(i,j,k) - tmpB4
                ELSE
                  dcflbydt(i,j,k)       = dcflbydt(i,j,k)             &
                  + (cf_liquid_conv(i,j,k) - tmpB4) / timestep_conv
                END IF
              END IF
            END IF

            ! Ice cloud fraction
            IF (cf_frozen_conv(i,j,k) > 0.0) THEN
              IF ((qcf_conv(i,j,k)/cf_frozen_conv(i,j,k))>2.0e-3) THEN
                tmpB4 = cf_frozen_conv(i,j,k)
                cf_frozen_conv(i,j,k) = MIN(1.0,qcf_conv(i,j,k)/2.0e-3)
                cf_frozen_inc(i,j,k)  = cf_frozen_inc(i,j,k)          &
                  + cf_frozen_conv(i,j,k) - tmpB4
                IF (l_fixbug_pc2_mixph) THEN
                  bulk_cf_conv(i,j,k) = bulk_cf_conv(i,j,k)           &
                  + cf_frozen_conv(i,j,k) - tmpB4
                  bulk_cf_inc(i,j,k) = bulk_cf_inc(i,j,k)             &
                  + cf_frozen_conv(i,j,k) - tmpB4
                ELSE
                  dcffbydt(i,j,k)       = dcffbydt(i,j,k)             &
                  + (cf_frozen_conv(i,j,k) - tmpB4) / timestep_conv
                END IF
              END IF
            END IF

          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

    END IF ! i_pc2_conv_coupling


    IF (isrfexcnvgust == ip_srfexwithcnv) THEN
      !           Store convective mass fluxes at cloud base
      !           if required for surface exchange.
      DO j = 1, rows
        DO i = 1, row_length
          IF (ccb(i,j) > 0) THEN
            ddmfx(i,j)=dwn_flux(i,j,ccb(i,j))
          ELSE
            ddmfx(i,j)=0.0
          END IF
        END DO
      END DO
    END IF

    ! diagnose number of convecting points
    IF (printstatus == prstatus_diag) THEN
      IF (n_conv_calls  >   1) THEN
        n_conv_points = 0
        DO j = 1, rows
          DO i = 1, row_length
            this_point = 0
            DO k = 1, n_conv_levels
              IF (ABS(dthbydt(i,j,k) * timestep)  >   0.0001) THEN
                this_point = 1
              END IF
            END DO
            n_conv_points = n_conv_points + this_point
          END DO
        END DO
        IF (n_proc  >   1) THEN
          CALL gc_isum(1, n_proc, info, n_conv_points )
        END IF
        IF (mype  ==  0) THEN
          WRITE(umMessage,'(a,i0,a,i0,a)')                            &
                     ' conv CALL ',call_number,' has ',               &
                     n_conv_points,' convecting points '
          CALL umPrint(umMessage,src='ni_conv_ctl')
        END IF

      END IF
    END IF

    !--------------------------------------------------------------------------
    ! Wind increments - interpolation back to uv grid now in atmos_physics2
    !                   so that swap-bound calls can be grouped with others
    !                   to save CPU.
    !--------------------------------------------------------------------------

  END DO ! loop over number of convection calls

  !============================================================================
  !  End of Convection sub-stepping loop
  !============================================================================

  ! 3D Convection prognostic of recent precip
  IF (l_conv_prog_precip) THEN

    decay_amount        = timestep/tau_conv_prog_precip
    theta_inc_threshold = dthetadt_conv_active_threshold * timestep_conv
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( rows, row_length, tot_conv_precip_2d, precip_shall,      &
!$OMP         precip_cong, precip_deep, precip_mid, theta_inc, tdims,  &
!$OMP         theta_inc_threshold, decay_amount, conv_prog_precip )    &
!$OMP PRIVATE( i, j, k, conv_active )
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        ! Total convective precip rate for the timestep / Kg m-2 s-1
        tot_conv_precip_2d(i,j) = MAX(precip_shall(i,j) +               &
                              precip_cong(i,j) +                        &
                              precip_deep(i,j)  + precip_mid(i,j),      &
                              conv_prog_precip_min_threshold         )
      END DO
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO k = 1, tdims%k_end
      DO j = 1, rows
        DO i = 1, row_length
          ! Define "active" convection to be where the magnitude of the 
          ! convective potential temperature tendency exceeds a threshold 
          ! (which is set in cv_param_mod)
          IF (ABS(theta_inc(i,j,k)) > theta_inc_threshold) THEN
            conv_active = 1.0
          ELSE
            conv_active = 0.0
          END IF
          ! Decay value held  - note will never reset to zero
          conv_prog_precip(i,j,k)                                         &
              = decay_amount * tot_conv_precip_2d(i,j) * conv_active      &
              + (1.0 - decay_amount) * conv_prog_precip(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

  END IF


  ! Deallocate arrays for SCM convection substep diagnostics
  ! (done after substepping so that we only allocate / deallocate once)
  IF (l_scm_convss_dg) THEN

    ! Deallocate the fields contained in the structures
    DO j = 1, rows
      DO i = 1, row_length
        CALL scm_convss_dg_deallocate( scm_convss_dg(i,j) )
      END DO
    END DO

  END IF

  ! Deallocate the array of structures
  DEALLOCATE( scm_convss_dg )


  !============================================================================
  ! Section 3. After convection substepping
  !============================================================================

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

  ! ----------------------------------------------------------------------
  ! 3.1 Check that CCA doesn't exceed 1.0
  ! ----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO k=1, n_cca_lev
    DO j=1, rows
      DO i=1, row_length
        cca(i,j,k) = MIN(cca(i,j,k), 1.0)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

  !==========================================================================
  ! Section 4 Treat tracers
  !==========================================================================

  !----------------------------------------------------------------------
  ! Section 4.1 Copy tracers back into variables from tot_tracers
  !----------------------------------------------------------------------
  IF ( l_tracer ) THEN

    CALL tracer_restore(row_length, rows, ntra_fld, ntra_lev,            &
          tr_vars, tr_ukca, tot_tracer,                                  &
          aerosol,                                                       &
          dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,   &
          so2, so4_aitken, so4_accu, so4_diss,                           &
          dms, nh3, soot_new, soot_agd, soot_cld, bmass_new, bmass_agd,  &
          bmass_cld, ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss,  &
          co2,  free_tracers, tracer_ukca,                               &
          ozone_tracer )

  END IF 

  ! Free up the memory again
  DEALLOCATE( tot_tracer )

  !----------------------------------------------------------------------
  ! Section 4.2 Scavenging for tracers
  !----------------------------------------------------------------------

  ! work with tracers only in final cycle
  IF ( cycleno == numcycles ) THEN

    CALL all_scav_calls(row_length, rows, ccb, cct,                    &
      conv_rain, conv_snow, rho, q_conv, qcl_conv, qcf_conv,           &
      aerosol,                                                         &
      dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,     &
      so2, so4_aitken, so4_accu, so4_diss,                             &
      dms, nh3, soot_agd, bmass_agd, ocff_agd, nitr_acc, nitr_diss  )

  END IF ! CycleNo == NumCycles

  ! ----------------------------------------------------------------------
  ! Copy convection increments into diagnostic arrays
  ! ----------------------------------------------------------------------
  IF ( l_apply_diag ) THEN
    CALL save_conv_diags( theta_inc, q_inc, qcl_inc, qcf_inc,          &
                          bulk_cf_inc, cf_liquid_inc, cf_frozen_inc,   &
                          exner_theta_levels, l_calc_dxek )
  END IF

  ! ----------------------------------------------------------------------
  ! Section 5 Calculate PC2 scheme increments to the increments.
  ! ----------------------------------------------------------------------

  ! L_calc_dxek_if2:
  IF (l_calc_dxek) THEN

    CALL pc2_from_conv_ctl(                                               &
                 row_length, rows, n_conv_levels,                         &
                 ! Model switches
                 l_calc_dxek, l_q_interact, l_mr_physics,               &
                 ! Coordinate info
                 exner_theta_levels, p_layer_centres,                     &  
                 ! in/out Convection increments
                 theta_inc, q_inc, qcl_inc, qcf_inc,                      &
                 cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,               &
                 ! in/out values with all increments added
                 theta_conv, q_conv, qcl_conv, qcf_conv,                  &
                 cf_liquid_conv, cf_frozen_conv, bulk_cf_conv )

  END IF  ! L_calc_dxek_if2


  ! ----------------------------------------------------------------------
  ! Section 6 Convective history - decay of convective precip
  ! ----------------------------------------------------------------------

  ! Convective history using 2D prognostic of recent precip
  IF (l_conv_hist) THEN

    decay_amount = timestep/decay_period
    DO  j = 1, rows
      DO i = 1, row_length
        ! total convective precip rate for the timestep
        tot_conv_precip = precip_shall(i,j) + precip_cong(i,j) +   &
                          precip_deep(i,j)  + precip_mid(i,j)
        ! If convective precip > 0.0 reset history
        IF (tot_conv_precip > 0.0) THEN
          past_precip(i,j) = tot_conv_precip
        ELSE
          ! Decay value held  - note will never reset to zero
          past_precip(i,j) = (1.0 - decay_amount)* past_precip(i,j)
        END IF
      END DO
    END DO

  END IF

  ! ----------------------------------------------------------------------
  ! Check moisture integrals if required, only call on final cycle 
  ! ----------------------------------------------------------------------
  IF (cycleno == numcycles .AND. l_check_moist_inc                         &
                                  .AND. model_type == mt_global) THEN

    ALLOCATE (dqt(row_length,rows,model_levels))

    DO j = 1, rows
      DO i=1,row_length
        ep4(i,j) = conv_rain(i,j)+conv_snow(i,j)
      END DO
    END DO
    DO k = 1, model_levels
      DO j = 1, rows
        DO i=1,row_length
          dqt(i,j,k) = (q_inc(i,j,k)+ qcl_inc(i,j,k)                       &
                            +qcf_inc(i,j,k))*recip_timestep
        END DO
      END DO
    END DO

    ! Working on specific humidities
    scheme_name =  'convection scheme'
    CALL check_dmoist_inc(row_length, rows, model_levels,                     &
                              delta_lambda,delta_phi, timestep ,              &
                              p_layer_boundaries,rho,                         &
                              dqt,                                            &
                              ep4, scheme_name,                               &
                              ep1, ep2, ep3)

    DEALLOCATE (dqt)
  END IF

  ! ----------------------------------------------------------------------
  ! Add the convection increments on to update the primary fields
  ! ----------------------------------------------------------------------

  ! If running the SCM with the convection in diagnostics-only mode
  ! (conv_mode=1), don't add on the convection increments.
  IF ( .NOT.( model_type==mt_single_column .AND. conv_mode == 1 ) ) THEN

!$OMP PARALLEL DEFAULT(NONE)                                        &
!$OMP SHARED(theta_star, theta_inc, q_star, q_inc, tdims,           &
!$OMP qcl_star, qcl_inc, qcf_star, qcf_inc,                         &
!$OMP cf_liquid_star, cf_frozen_star, bulk_cf_star,                 &
!$OMP cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,                    &
!$OMP L_calc_dxek, L_q_interact)                                    &
!$OMP PRIVATE(i,j,k)

!$OMP DO SCHEDULE(STATIC)
    DO k =                 1, tdims%k_end
      DO j =   tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          theta_star(i,j,k) = theta_star(i,j,k) + theta_inc(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
    DO k =                 1, tdims%k_end
      DO j =   tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_star(i,j,k) = q_star(i,j,k) + q_inc(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    ! ----------------------------------------------------------------------
    ! Protected loop. Update increments only when PC2 scheme is fully ON-ON.
    ! ----------------------------------------------------------------------
    !
    ! L_calc_dxek_if1:
    IF ( L_calc_dxek .AND. L_q_interact ) THEN

!$OMP DO SCHEDULE(STATIC)
      DO k =                 1, tdims%k_end
        DO j =   tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            qcl_star(i,j,k) = qcl_star(i,j,k) + qcl_inc(i,j,k)
            qcf_star(i,j,k) = qcf_star(i,j,k) + qcf_inc(i,j,k)
            cf_liquid_star(i,j,k) = cf_liquid_star(i,j,k)            &
                                  + cf_liquid_inc(i,j,k)
            cf_frozen_star(i,j,k) = cf_frozen_star(i,j,k)            &
                                  + cf_frozen_inc(i,j,k)
            bulk_cf_star(i,j,k)   = bulk_cf_star(i,j,k)              &
                                  + bulk_cf_inc(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF  ! L_calc_dxek_if1

  ! Other prognostics if present are not updated by convection so no
  ! copying back to qstar variables required.
  ! i.e. qrain, qgraup, qcf2
  ! Note that u, v are updated in atmos_physics2, after the convective
  ! momentum tendencies have been interpolated onto u / v points.

!$OMP END PARALLEL

  END IF  ! not SCM running convection in diagnostics-only mode


  !----------------------------------------------------------------------------
  ! end of routine NI_conv_ctl
  !----------------------------------------------------------------------------
END IF ! on error code equal to zero

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE ni_conv_ctl

END MODULE ni_conv_ctl_mod
