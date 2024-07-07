! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Interface to Atmospheric Physics radiation code

! Purpose:
!   This is the top-level radiation control routine. Major book-keeping
!   is carried out here. Separate calls are made to calculate SW and LW
!   radiative fluxes and heating rates on radiation time-steps.
!   Radiation increments are applied on all physics time-steps.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

!------------------------------------------------------------------------------
MODULE rad_ctl_mod


USE timer_mod, ONLY: timer

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'RAD_CTL_MOD'
CONTAINS

SUBROUTINE rad_ctl (                                                    &

! Parallel variables
  off_x, off_y, at_extremity, n_proc, global_cloud_top,                 &

! model dimensions.
  row_length, rows,                                                     &
  bl_levels,ozone_levels, cloud_levels, n_cca_levels,                   &
  ntiles, land_field, dust_dim1, dust_dim2,                             &
  biogenic_dim1,biogenic_dim2, sulp_dim1, sulp_dim2, soot_dim1,         &
  soot_dim2,bmass_dim1, bmass_dim2, salt_dim1, salt_dim2, salt_dim3,    &
  co2_dim_len, co2_dim_row, co2_dim2, arcl_dim1, arcl_dim2,             &
  n_arcl_species, n_arcl_compnts, i_arcl_compnts,                       &
  ocff_dim1, ocff_dim2, nitrate_dim1, nitrate_dim2,                     &
  ukca_dim1, ukca_dim2, n_ukca_mode, n_ukca_cpnt,                       &

! Model switches
  l_emcorr, l_snow_albedo,                                              &
  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                                &
  l_murk_rad, l_co2_interactive,                                        &
  l_use_arcl,                                                           &
  l_dust, l_sulpc_so2, l_soot, l_biomass, l_ocff, l_nitrate,            &
  l_cosp_in,                                                            &
  l_easyaerosol_sw, l_easyaerosol_lw, l_easyaerosol_cdnc,               &

! in coordinate information
  rho_r2, trindx,                                                       &

! in time stepping information.
  timestep_number, previous_time,                                       &

! diagnostic info
  stashwork1,                                                           &
  stashwork2,                                                           &

! SCM diagnostics switches (dummy in full UM)
  nscmdpkgs,l_scmdiags,                                                 &

! in data fields.
  p_star, p_layer_boundaries, p_layer_centres,                          &
  land_sea_mask_in, fland, land0p5_in, sd_orog_land,                    &
  t_surf, tstar_sea, tstar_sice_cat, area_cloud_fraction,               &
  dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,          &
  biogenic, so4_aitken, so4_accu, so4_diss,                             &
  soot_new, soot_agd, bmass_new, bmass_agd, bmass_cld,                  &
  ocff_new, ocff_agd, ocff_cld, nitr_acc, nitr_diss, aerosol, arcl,     &
  ukca_radaer, sea_salt_film, sea_salt_jet, ws_10m_sea,                 &
  co2_3d, frac_control, n_drop_pot,                                     &
  easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc,                     &

! chemical greenhouse gas fields
  ngrgas, grgas_field,                                                  &

! ancillary fields and fields needed to be kept from timestep to timestep
  snow_depth, snow_depth_sea_cat, ice_fract_in, ice_fract_cat,          &
  ice_thick_cat, chloro_sea, rgrain, soot, canht,                       &
  cca, ccb, cct, cclwp,ccw,lcbase,                                      &
  ozone, sw_incs, lw_incs,                                              &
  o3_trop_level, o3_trop_height,                                        &
  t_trop_level, t_trop_height, zh,                                      &
  land_index, albsoil, albobs_sw, albobs_vis, albobs_nir,               &
  lai, snow_tile, frac, tstar_tile, z0_tile,                            &
  dOLR_rts, lw_down, sw_tile_rts,                                       &
  land_alb,sice_alb,                                                    &
  es_space_interp, rad_mask,                                            &

! in/out
  t_n, q_n, qcl_n, qcf_n, cf_n, cfl_n, cff_n,                           &
  qcf2_n, qrain_n, qgraup_n,                                            &
  t_inc, q_inc, qcl_inc, cf_inc, cfl_inc,                               &
  sum_eng_fluxes, cos_zenith_angle,                                     &

! out.
  photosynth_act_rad, rad_hr, dOLR, sw_tile,                            &

! COSP arguments
  cosp_gbx, cosp_sgx, cosp_sgh,                                         &

! error information
  error_code  )

USE missing_data_mod, ONLY: rmdi
USE jules_surface_types_mod
USE planet_constants_mod, ONLY: l_planet_g, g, sc,                      &
  l_planet_grey_surface, l_planet_intrinsic_flux, planet_t_intrinsic,   &
  planet_emissivity, planet_radius, stellar_radius
USE conversions_mod, ONLY: pi
USE csigma,      ONLY: sbcon
USE trignometric_mod, ONLY:                                             &
  true_longitude, true_latitude,                                        &
  cos_theta_latitude
USE horiz_grid_mod, ONLY: cell_area_surface
USE rad_ccf, ONLY: astronomical_unit
USE level_heights_mod, ONLY: r_layer_centres, r_layer_boundaries
USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE nlstcall_mod, ONLY: ltimer
! Modules required for multiple calls to radiation
USE diagoffset
USE gen_phys_inputs_mod, ONLY: l_mr_physics

! Modules required for reading in spectral files and control options
! as structures
USE spec_sw_lw, ONLY: sw_spectrum, lw_spectrum
USE sw_control_struct, ONLY: sw_control
USE lw_control_struct, ONLY: lw_control
USE def_control, ONLY: deallocate_control

! Structure definitions for the core radiation code interface
USE def_dimen, ONLY: StrDim

! Modules required for radiative forcing
USE coradoca, ONLY: c2c_o2, c2c_o3, c2c_co2, c2c_n2o, c2c_ch4,          &
  c2c_cfc11, c2c_cfc12, c2c_c113, c2c_hcfc22, c2c_hfc125,               &
  c2c_hfc134, c2c_aerosol, c2c_sulpc_d, c2c_seas_d, c2c_soot_d,         &
  c2c_bmb_d, c2c_ocff_d, c2c_land_s, c2c_all,                           &
  c2c_wmg, c2c_nitr_d, c2c_dust_d, c2c_biog_d, c2c_ukca_d, c2c_easy_d,  &
  co2_mmr_scl, co2_mmr_add, n2o_mmr_scl, n2o_mmr_add,                   &
  ch4_mmr_scl, ch4_mmr_add, o2_mmr_scl, o2_mmr_add,                     &
  cfc11_mmr_scl, cfc11_mmr_add, cfc12_mmr_scl, cfc12_mmr_add,           &
  cfc113_mmr_scl, cfc113_mmr_add, hcfc22_mmr_scl, hcfc22_mmr_add,       &
  hfc125_mmr_scl, hfc125_mmr_add, hfc134a_mmr_scl, hfc134a_mmr_add
USE ukca_feedback_mod, ONLY: p_n2o, p_ch4, p_f11, p_f12, p_f113, p_f22

! Modules required for the orography scheme
USE solinc_data, ONLY:                                                  &
  orog_corr, f_orog, slope_aspect, slope_angle, horiz_ang, horiz_aspect,&
  n_horiz_layer, n_horiz_ang, l_orog, l_skyview

! Modules for diagnostics
USE sw_diag_mod, ONLY: sw_diag
USE lw_diag_mod, ONLY: lw_diag
USE def_diag, ONLY: init_diag_logic, allocate_diag, deallocate_diag

! Module for radiation switches
USE rad_input_mod, ONLY: l_rad_perturb, l_rad_szacor,                   &
  l_rad_snow_emis, l_t_land_nosnow, l_quad_t_coast,                     &
  l_t_rad_solid, l_forcing, l_timestep, l_radiance,                     &
  l_inhom_cloud, l_rad_deg, l_extra_top, l_t_bdy_surf,                  &
  l_use_dust, l_use_biogenic, l_use_sulpc_direct,                       &
  l_use_soot_direct, l_use_soot_indirect,                               &
  l_use_bmass_direct, l_use_bmass_indirect,                             &
  l_use_ocff_direct, l_use_ocff_indirect,                               &
  l_use_sulpc_indirect_sw, l_use_sulpc_indirect_lw,                     &
  l_use_nitrate_direct, l_use_nitrate_indirect,                         &
  l_use_seasalt_indirect, l_use_seasalt_direct,                         &
  l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,               &
  inhom_cloud_sw, inhom_cloud_lw, dp_corr_strat, dp_corr_conv,          &
  co2_mmr, o2mmr, n2ommr, ch4mmr, so2mmr, c11mmr, c12mmr,               &
  c113mmr, c114mmr, hcfc22mmr, hfc125mmr, hfc134ammr,                   &
  aero_bl_levels,                                                       &
  a_sw_radstep_diag, a_sw_radstep_prog,                                 &
  a_lw_radstep_diag, a_lw_radstep_prog,                                 &
  n_swcall, n_lwcall

USE tuning_segments_mod, ONLY:                                          &
  l_autotune_segments,                                                  &
  a_sw_segments, a_sw_seg_size,                                         &
  a_lw_segments, a_lw_seg_size

! Flags for radiation timesteps
USE set_rad_steps_mod, ONLY: l_rad_step_prog, l_rad_step_diag

! Radiation timestep length
USE timestep_mod, ONLY: radiation_tstep_diag, radiation_tstep_prog,     &
                        timestep, recip_timestep

! Solar spectral variation
USE solvar_mod, ONLY: solvar

! Sperical path for solar beam
USE solang_sph_mod, ONLY: solang_sph

! Modules for JULES
!Subroutines
USE tilepts_mod, ONLY: tilepts
USE surf_couple_radiation_mod, ONLY: surf_couple_radiation

!Variables
USE jules_surface_mod, ONLY: l_aggregate, l_flake_model
USE jules_sea_seaice_mod, ONLY: l_ctile
USE jules_radiation_mod, ONLY: l_dolr_land_black, l_embedded_snow

USE nvegparm, ONLY: emis_nvg
USE pftparm,  ONLY: emis_pft
USE jules_sea_seaice_mod,  ONLY: nice_use, emis_sea, emis_sice
USE jules_snow_mod, ONLY: rho_snow_const
USE ancil_info, ONLY:                                                   &
  ssi_pts, sice_pts, sice_pts_ncat,                                     &
  ssi_index, sice_index, sice_index_ncat, sea_index, sea_pts
USE lake_mod, ONLY: lake_albedo_gb

USE fluxes, ONLY: sw_sicat, sw_rts_sicat, swup_rts_sicat,               &
                  swdn_rts_sicat, alb_sicat, sw_sea, sw_rts_sea


USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

! Module for UKCA-MODE aerosol interaction with radiation
USE ukca_option_mod, ONLY: l_ukca_radaer
USE ukca_radaer_struct_mod, ONLY: ukca_radaer_struct
USE glomap_clim_option_mod, ONLY: l_glomap_clim_radaer

! Module for two-bin or six-bin dust switch
USE dust_parameters_mod, ONLY: l_twobin_dust

! Module for the aerosol climatologies
USE arcl_mod, ONLY: npd_arcl_compnts, npd_arcl_species,                 &
  ip_arcl_sulp_ac, ip_arcl_sulp_ak

USE ereport_mod, ONLY: ereport
USE um_parparams, ONLY: pnorth, psouth
USE Field_Types

! Modules for COSP
USE cosp_types_mod, ONLY: cosp_gridbox, cosp_subgrid, cosp_sghydro

USE segments_mod, ONLY:                                                 &
   segment_type, meta_segment_type,                                     &
   segments_mod_seg_meta, segments_mod_segments

USE autotune_mod, ONLY:       &
    autotune_type,            &
    autotune_init,            &
    autotune_get_trial_size,  &
    autotune_start_region,    &
    autotune_stop_region,     &
    autotune_advance,         &
    autotune_report

USE physics_tendencies_mod,  ONLY:                                      &
    l_retain_rad_tendencies,dt_sw, dq_sw, dt_lw, dq_lw

USE nlsizes_namelist_mod, ONLY: model_levels
USE stash_array_mod, ONLY: sf_calc, nitems, nsects
USE stparam_mod, ONLY: s_modl, s_sect, s_item
USE atm_fields_mod, ONLY: pond_frac_cat, pond_depth_cat
USE atm_fields_bounds_mod, ONLY: ScmRowLen, ScmRow
USE s_scmop_mod,           ONLY: default_streams,                       &
                                 t_avg, t_acc, t_mult, only_radsteps,   &
                                 d_sl, d_wet, d_all, d_cloud,           &
                                 scmdiag_rad, scmdiag_pc2
USE scmoutput_mod, ONLY: scmoutput

USE def_easyaerosol, ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

!$ USE omp_lib

USE errormessagelength_mod, ONLY: errormessagelength
USE model_domain_mod, ONLY: model_type, mt_global, mt_single_column
USE set_thermodynamic_mod, ONLY: set_thermodynamic

USE flux_diag_mod, ONLY: flux_diag
USE fill_missing_data_lw_mod, ONLY: fill_missing_data_lw
USE fill_missing_data_sw_mod, ONLY: fill_missing_data_sw
USE lw_rad_mod, ONLY: lw_rad
USE prelim_lwrad_mod, ONLY: prelim_lwrad
USE prelim_swrad_mod, ONLY: prelim_swrad
USE set_control_mod, ONLY: set_control
USE set_dimen_mod, ONLY: set_dimen
USE set_lwdiag_logic_mod, ONLY: set_lwdiag_logic
USE set_swdiag_logic_mod, ONLY: set_swdiag_logic
USE solang_mod, ONLY: solang
USE solinc_mod, ONLY: solinc
USE solpos_mod, ONLY: solpos
USE sw_rad_mod, ONLY: sw_rad
USE pc2_homog_plus_turb_mod, ONLY: pc2_homog_plus_turb
USE diagnostics_rad_mod, ONLY: diagnostics_rad

IMPLICIT NONE

! Segmentation variables
TYPE(segment_type),     ALLOCATABLE  :: segments(:)
TYPE(meta_segment_type)              :: meta_segments
INTEGER                              :: ipar
INTEGER                              :: num_parallel_sw
INTEGER                              :: num_parallel_lw

! Arguments with intent in. ie: input variables.

! Parallel setup variables
INTEGER ::                                                              &
  off_x,                                                                &
!   Size of small halo in i
  off_y,                                                                &
!   Size of small halo in j.
  n_proc,                                                               &
!   Total number of processors
  global_cloud_top
!   Global topmost cloudy level

LOGICAL ::                                                              &
  at_extremity(4)
!   Indicates if this processor is at north, south,
!   east or west of the processor grid

! Model dimensions
INTEGER ::                                                              &
  row_length,                                                           &
  rows,                                                                 &
  bl_levels,                                                            &
  ozone_levels,                                                         &
  cloud_levels,                                                         &
  n_cca_levels,                                                         &
  ntiles,                                                               &
  land_field,                                                           &
!   dimensions of mineral dust arrays
  dust_dim1,                                                            &
  dust_dim2,                                                            &
!   dimensions of biogenic aerosol arrays
  biogenic_dim1,                                                        &
  biogenic_dim2,                                                        &
!   dimensions of S Cyc arrays
  sulp_dim1,                                                            &
  sulp_dim2,                                                            &
!   dimensions of soot arrays
  soot_dim1,                                                            &
  soot_dim2,                                                            &
!   dimensions of biomass arrays
  bmass_dim1,                                                           &
  bmass_dim2,                                                           &
!   dimensions of OCFF arrays
  ocff_dim1,                                                            &
  ocff_dim2,                                                            &
!   dimensions of sea-salt arrays
  salt_dim1,                                                            &
  salt_dim2,                                                            &
  salt_dim3,                                                            &
!   dimensions of CO2 array
  co2_dim1,                                                             &
  co2_dim2,                                                             &
  co2_dim_len,                                                          &
  co2_dim_row,                                                          &
!   dimensions of aerosol clim for NWP
  arcl_dim1,                                                            &
  arcl_dim2,                                                            &
!   dimensions of nitrate arrays
  nitrate_dim1,                                                         &
  nitrate_dim2,                                                         &
!   dimensions of UKCA_RADAER arrays
  ukca_dim1,                                                            &
  ukca_dim2,                                                            &
  n_ukca_mode,                                                          &
  n_ukca_cpnt

LOGICAL ::                                                              &
  l_dust,                                                               &
!   mineral dust available for use (for direct effect or diagnostics)
  l_sulpc_so2,                                                          &
!   Sulphur C available for use (for direct/indirect or diagnostics)
  l_soot,                                                               &
!   Soot available for use (for direct/indirect or diagnostics)
  l_biomass,                                                            &
!   biomass smoke available for use (for direct/indirect or diagnostics)
  l_ocff,                                                               &
!   OCFF available for use (for direct/indirect or diagnostics)
  l_nitrate,                                                            &
!   nitrate available for use (for direct/indirect or diagnostics)
  l_emcorr,                                                             &
!   true if energy correction scheme is to be used.
  l_snow_albedo,                                                        &
!   True if spectral albedo scheme selected
  l_mcr_qcf2,                                                           &
!   Use second ice category
  l_mcr_qrain,                                                          &
!   Use prognostic rain
  l_mcr_qgraup,                                                         &
!   Use graupel
  l_murk_rad,                                                           &
!   True if using radiative effects of 'murk'.
  l_co2_interactive
!   Controls the use of 3D CO2 field

! Is COSP requested?
LOGICAL, INTENT(IN) :: l_cosp_in

! Diagnostics info
REAL ::                                                                 &
  stashwork1(*),                                                        &
  stashwork2(*)
!   STASH workspace


! Data arrays
REAL ::                                                                 &
  p_layer_boundaries(row_length, rows, 0:model_levels),                 &
!   pressure at layer boundaries. Same as p except at
!   bottom level = pstar, and at top = 0.
  p_layer_centres(row_length, rows, 0:model_levels),                    &
!   pressure at layer centres. Same as p_theta_levels
!   except zeroth level = pstar.
  p_star(row_length, rows)

REAL    :: ccw   (row_length, rows, model_levels)
INTEGER :: lcbase(row_length, rows)

REAL ::                                                                 &
  cca (row_length, rows, n_cca_levels),                                 &
  cclwp(row_length, rows),                                              &
!   condensed water path (KG/M**2)
  area_cloud_fraction(row_length, rows, model_levels)

INTEGER ::                                                              &
  ccb (row_length, rows),                                               &
  cct (row_length, rows)

! ancillary arrays and fields required to be saved from timestep to
! timestep.

REAL ::                                                                 &
  t_surf(row_length, rows),                                             &
  tstar_sea(row_length,rows),                                           &
!   IN Open sea sfc temperature (K).
  tstar_sice_cat(row_length,rows,nice_use),                             &
!   IN Sea-ice sfc temperature (K).
  chloro_sea(row_length, rows)
!   IN sea nr. surface chlorophyll

REAL ::                                                                 &
  ice_fract_cat(row_length, rows, nice_use),                            &
!   Area fraction of sea ice categories
  ice_thick_cat(row_length, rows, nice_use)
!   Effective thickness of each sea ice categories


LOGICAL ::                                                              &
  land_sea_mask_in(row_length, rows),                                   &
  land0p5_in(row_length, rows),                                         &
!   A mask set to .TRUE. if the fraction of land in the grid-box
!   exceeds 0.5.
  rad_mask(row_length, rows)
!   A mask which ensures a chequerboard pattern of radiation
!   calculations over the whole domain (not just one PE)

REAL ::                                                                 &
  fland(land_field)
!   Fractional amount of land at each land point

 REAL ::                                                                &
  sd_orog_land(land_field)
!     Standard deviation of the orography


REAL ::                                                                 &
  snow_depth (row_length, rows),                                        &
!   snow/qrclim.snow.(month)
  snow_depth_sea_cat (row_length, rows, nice_use)
!   snow depth on sea ice categories

REAL ::                                                                 &
  ice_fract_in (row_length, rows),                                      &
!   ice/qrclim.ice.(month)
  soot (row_length, rows)

! Input ancillary data:
REAL, INTENT(IN) ::                                                     &
  albsoil(land_field),                                                  &
  albobs_sw(land_field),                                                &
  albobs_vis(land_field),                                               &
  albobs_nir(land_field),                                               &
  lai(land_field, npft),                                                &
  canht(land_field, npft)

REAL ::                                                                 &
  rgrain(land_field, ntiles),                                           &
  snow_tile(land_field, ntiles),                                        &
  frac(land_field, ntype),                                              &
  frac_control(land_field, ntype),                                      &
  tstar_tile(land_field, ntiles),                                       &
  z0_tile(land_field, ntiles),                                          &
  dOLR_rts(row_length, rows),                                           &
!   TOA - surface upward LW
  lw_down(row_length, rows),                                            &
!   Surface downward LW
  sw_tile_rts(land_field, ntiles),                                      &
!   Surface net SW on land tiles
  land_alb(row_length, rows),                                           &
!   Mean land albedo
  sice_alb(row_length, rows),                                           &
!   Mean sea-ice albedo
  es_space_interp(4, row_length, rows)
!   Coeffs for spatial interpolation of radiation quantities

! Number of requested species within the climatology
INTEGER :: n_arcl_species

! Corresponding number of requested components
INTEGER :: n_arcl_compnts

! Model switch for each species
LOGICAL :: l_use_arcl(npd_arcl_species)

! Array indices of components
INTEGER :: i_arcl_compnts(npd_arcl_compnts)

! Mass-mixing ratios
REAL ::                                                                 &
  arcl(row_length, rows, model_levels, n_arcl_compnts)

! UKCA_RADAER structure: Interaction between UKCA-MODE aerosols and radiation
TYPE (ukca_radaer_struct), INTENT(IN) :: ukca_radaer

! EasyAerosol climatology
!
! Model switch for each spectrum
LOGICAL, INTENT(IN) :: l_easyaerosol_sw
LOGICAL, INTENT(IN) :: l_easyaerosol_lw

! Model switch for cloud droplet number concentrations
LOGICAL, INTENT(IN) :: l_easyaerosol_cdnc

! EasyAerosol distributions
TYPE (t_easyaerosol_rad), INTENT(IN) :: easyaerosol_sw
TYPE (t_easyaerosol_rad), INTENT(IN) :: easyaerosol_lw
TYPE (t_easyaerosol_cdnc), INTENT(IN) :: easyaerosol_cdnc

INTEGER ::                                                              &
  land_index(land_field)

REAL ::                                                                 &
  ozone(row_length, rows, ozone_levels),                                &
  o3_trop_level(row_length,rows),                                       &
  o3_trop_height(row_length,rows),                                      &
  t_trop_level(row_length,rows),                                        &
  t_trop_height(row_length,rows),                                       &
  sw_incs(row_length, rows, 0:model_levels+1),                          &
  lw_incs(row_length, rows, 0:model_levels),                            &
  zh(row_length, rows),                                                 &
!   boundary layer height
  co2_3d(row_length, rows, model_levels)

! Potential droplet number (includes values where cloud not present)
REAL, INTENT(IN) :: n_drop_pot(row_length, rows, model_levels)

! chemical greenhouse gas fields
INTEGER, INTENT(IN) :: ngrgas
REAL, INTENT(INOUT) :: grgas_field(row_length,rows,model_levels,ngrgas)

! Co-ordinate arrays
REAL ::                                                                 &
  rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)
!   Density*radius^2

! Allocatable arrays for diagnostic variables - required to save memory
! use during radiation routines
REAL, ALLOCATABLE :: t_incr_diagnostic(:,:,:)
!   temperature increment for STASH


! time information for current timestep
INTEGER ::                                                              &
  timestep_number,                                                      &
  previous_time(7)

! Diagnostic variables
! Level of tropopause
INTEGER, INTENT(IN) :: trindx(row_length, rows)

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER ::                                                              &
  nscmdpkgs
!   No of diagnostics packages

LOGICAL ::                                                              &
  l_scmdiags(nscmdpkgs)
!   Logicals for diagnostics packages

INTEGER ::                                                              &
  error_code

! Variables with intent (in/out)

REAL ::                                                                 &
  t_n(row_length, rows, model_levels),                                  &
  q_n(row_length, rows, model_levels),                                  &
  qcl_n(row_length, rows, model_levels),                                &
  qcf_n(row_length, rows, model_levels),                                &
  qcf2_n(row_length, rows, model_levels),                               &
!   2nd ice prog
  qrain_n(row_length, rows, model_levels),                              &
!   Rain prognostic
  qgraup_n(row_length, rows, model_levels),                             &
!   Graupel
  cf_n(row_length, rows, model_levels),                                 &
  cfl_n(row_length, rows, model_levels),                                &
  cff_n(row_length, rows, model_levels),                                &
  t_inc(row_length, rows, model_levels),                                &
  q_inc(row_length, rows, model_levels),                                &
  qcl_inc(row_length, rows, model_levels),                              &
  cf_inc(row_length, rows, model_levels),                               &
  cfl_inc(row_length, rows, model_levels),                              &
  sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                       &
  sea_salt_jet(salt_dim1, salt_dim2, salt_dim3),                        &
  ws_10m_sea(row_length, rows),                                         &
  sum_eng_fluxes(row_length, rows)

REAL, INTENT(INOUT) ::                                                  &
  so4_aitken(1-off_x:row_length+off_x, 1-off_y:rows+off_y,              &
                                               model_levels),           &
  so4_accu(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  so4_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  soot_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  soot_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  bmass_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                               model_levels),           &
  bmass_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                               model_levels),           &
  bmass_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                               model_levels),           &
  ocff_new(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  ocff_agd(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  ocff_cld(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  nitr_acc(1-off_x:row_length+off_x, 1-off_y:rows+off_y,                &
                                               model_levels),           &
  nitr_diss(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                               model_levels),           &
  dust_div1(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                              model_levels),            &
  dust_div2(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                              model_levels),            &
  dust_div3(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                              model_levels),            &
  dust_div4(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                              model_levels),            &
  dust_div5(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                              model_levels),            &
  dust_div6(1-off_x:row_length+off_x, 1-off_y:rows+off_y,               &
                                              model_levels),            &
  aerosol(row_length, rows, model_levels)

REAL, INTENT(IN) ::                                                     &
  biogenic(row_length, rows, model_levels)

! Variables with intent out
REAL ::                                                                 &
  photosynth_act_rad(row_length, rows),                                 &
!   Net downward shortwave radiation in band 1 (w/m2).
  rad_hr(row_length, rows, 2, bl_levels),                               &
!   BL radiative (LW,SW) heating rates
  dOLR(row_length, rows),                                               &
!   TOA - surface upward LW
  sw_tile(land_field, ntiles),                                          &
!   Surface net SW on land tiles
  cos_zenith_angle(row_length,rows)

! COSP variables
TYPE(cosp_gridbox),INTENT(INOUT) :: cosp_gbx
TYPE(cosp_subgrid),INTENT(INOUT) :: cosp_sgx
TYPE(cosp_sghydro),INTENT(INOUT) :: cosp_sgh


! local variables.
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RAD_CTL'
CHARACTER (LEN=errormessagelength) :: cmessage

LOGICAL :: TmpLogic ! Temporary logical
REAL :: dummy2d(1,1), dummy3d(1,1,1)

! loop counters
INTEGER ::                                                              &
  i, j, k, l, n, point

INTEGER :: ij

INTEGER ::                                                              &
  ii, jj

LOGICAL ::                                                              &
  land_sea_mask(row_length, rows),                                      &
  land0p5(row_length, rows)
!   A mask set to .TRUE. if the fraction of land in the grid-box
!   exceeds 0.5.

REAL ::                                                                 &
  ice_fract(row_length, rows)
!   Sea-ice fraction

INTEGER ::                                                              &
  daylight_points,                                                      &
  lit_points,                                                           &
  start_point,                                                          &
  first_point,                                                          &
  ptr_local,                                                            &
!   Pointer for LW arrays used to reduce length
!   of lines, repeating current element of segments(i)%fp
  list_lw_points(row_length*rows),                                      &
  lw_points,                                                            &
!   Variables to enable scatter/gather of LW (like SW): creates
!   a LIST of points to do a calculation, and also
!   facilitates segmenting.
  first_data_interp,                                                    &
!   The first data point, in data co-ords, which needs to be
!   interpolated on a PE (for use in interpolation routine)
  first_row, last_row,                                                  &
  tot_daylight_points,                                                  &
!   Total number of daylight points in whole domain
  max_n_swbands
!   Maximum number of sw_bands, in all sw calls:

! Local arrays holding information to be passed between physics
! routines.

! The following arrays are required for calculating increments
! and heating rates for SCM diagnostic output
REAL :: TmpScm2d  (ScmRowLen,ScmRow)
REAL :: TmpScm3d_1(ScmRowLen,ScmRow,model_levels)
REAL :: TmpScm3d_2(ScmRowLen,ScmRow,model_levels)


! The following arrays are required by the PC2 cloud scheme
REAL ::                                                                 &
  t_latest(row_length,rows,model_levels),                               &
  q_latest(row_length,rows,model_levels),                               &
  qcl_latest(row_length,rows,model_levels),                             &
  cf_latest(row_length,rows,model_levels),                              &
  cfl_latest(row_length,rows,model_levels),                             &
  cff_latest(row_length,rows,model_levels),                             &
  delta_t(row_length,rows,model_levels),                                &
  zeros(row_length,rows,model_levels)

! For a total ice, used to include both qcf2 species
REAL :: qcf_total(row_length, rows, model_levels)

! For calculation of flux diagnostics weighted by sea ice fraction:
REAL ::                                                                 &
  swup_sicat(ssi_pts,nice_use),                                         &
  swdn_sicat(ssi_pts,nice_use)

! Fields calculated by set_thermodynamic:
REAL :: t_layer_boundaries(row_length, rows, 0:model_levels)
!   Temperature at layer boundaries
REAL :: p_extra_layer(row_length, rows)
!   Pressure at centre of extra top layer
REAL :: t_extra_layer(row_length, rows)
!   Temperature at centre of extra top layer
REAL :: d_mass(row_length, rows, model_levels+1)
!   Mass of layer (kg m-2)
REAL :: density(row_length, rows, model_levels+1)
!   Density of layer (kg m-3)
REAL :: layer_heat_capacity(row_length, rows, model_levels)
!   Heat capacity of layer

! For albedo scaling diagnostic:
REAL ::                                                                 &
  albobs_sc(row_length, rows, ntiles, 2)

! Radiation fields 1. SW & common with LW.

REAL, SAVE :: solcon_rts
REAL, ALLOCATABLE, SAVE :: cos_zen_rts(:,:)
REAL, ALLOCATABLE, SAVE :: sol_azm_rts(:,:)
REAL, ALLOCATABLE, SAVE :: day_frac_rts(:,:)
REAL, ALLOCATABLE, SAVE :: surf_down_sw_rts(:,:,:)
REAL, ALLOCATABLE, SAVE :: alb_tile(:,:,:)
REAL :: surf_down_sw(row_length,rows,4)
REAL :: surf_down_sw_sum(row_length,rows)
REAL :: surfdir_rts(row_length,rows)
REAL :: eps
REAL :: sza_cor(row_length,rows)

! Fields for treatment of spherical geometry for the solar beam
REAL, ALLOCATABLE, SAVE :: cos_zen_sph_rts(:,:,:)
REAL, ALLOCATABLE, SAVE :: day_frac_sph_rts(:,:,:)

! Fields for zenith angle correction diagnostics:
LOGICAL :: l_surfsw_cor, l_toasw_cor
LOGICAL :: l_surfdir_cor, l_surfdif_cor
REAL :: surfsw_cor(row_length,rows)
REAL :: toasw_cor(row_length,rows)
REAL :: surfdir_cor(row_length,rows)
REAL :: surfdif_cor(row_length,rows)

LOGICAL :: l_swdn_sice_wt_cat, l_swdn_sice_wt, l_swup_sice_wt_cat,    &
           l_swup_sice_wt, l_alb_sice_wt_cat, l_alb_sice_wt,          &
           l_lwdn_sice_wt_cat, l_lwdn_sice_wt

REAL ::                                                                        &
  day_fraction(row_length,rows),                                               &
  mean_cos_zenith_angle(row_length,rows),                                      &
  solar_constant,                                                              &
  flandg(row_length, rows),                                                    &
  land_albedo(row_length, rows, 4, n_swcall),                                  &
  sea_ice_albedo(row_length, rows, 4, n_swcall),                               &
  netsw(row_length, rows),                                                     &
!   Net short-wave absorbed by planet
  swsea(row_length, rows),                                                     &
!   Net short-wave absorbed by sea
  dust_1(dust_dim1,dust_dim2),                                          &
  dust_2(dust_dim1,dust_dim2),                                          &
  dust_3(dust_dim1,dust_dim2),                                          &
  dust_4(dust_dim1,dust_dim2),                                          &
  dust_5(dust_dim1,dust_dim2),                                          &
  dust_6(dust_dim1,dust_dim2),                                          &
  local_biogenic(biogenic_dim1,biogenic_dim2),                          &
  accum_sulphate(sulp_dim1,sulp_dim2),                                  &
  aitken_sulphate(sulp_dim1,sulp_dim2),                                 &
  diss_sulphate(sulp_dim1,sulp_dim2),                                   &
  fresh_soot(soot_dim1,soot_dim2),aged_soot(soot_dim1,soot_dim2),       &
  fresh_bmass(bmass_dim1,bmass_dim2),                                   &
  aged_bmass(bmass_dim1,bmass_dim2),                                    &
  cloud_bmass(bmass_dim1,bmass_dim2),                                   &
  fresh_ocff(ocff_dim1, ocff_dim2),                                     &
  aged_ocff(ocff_dim1, ocff_dim2),                                      &
  cloud_ocff(ocff_dim1, ocff_dim2),                                     &
  accum_nitrate(nitrate_dim1, nitrate_dim2),                            &
  diss_nitrate(nitrate_dim1, nitrate_dim2),                             &
  sw_net_land(row_length, rows),                                        &
!   SW net local flux over land
  sw_net_sice(row_length, rows)
!   SW net local flux over sea-ice

! SW and LW diagnostics weighted by sea ice concentration:
REAL, ALLOCATABLE ::                                                    &
   sw_up_sice_weighted_cat(:,:,:),                                      &
   sw_down_sice_weighted_cat(:,:,:),                                    &
   sw_up_sice_weighted(:,:),                                            &
   sw_down_sice_weighted(:,:),                                          &
   albedo_sice_weighted_cat(:,:,:),                                     &
   albedo_sice_weighted(:,:),                                           &
   lw_down_sice_weighted_cat(:,:,:),                                    &
   lw_down_sice_weighted(:,:)

! Aerosol climatology for NWP
REAL ::                                                                 &
  local_arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)
REAL ::                                                                 &
  arcl_multf
!   Mass-mixing ratio conversion factor

! UKCA_RADAER: Interaction UKCA-MODE aerosols and radiation
REAL :: local_ukca_mmr(ukca_dim1, ukca_dim2, MAX(1,n_ukca_cpnt))
REAL :: local_ukca_cvl(ukca_dim1, ukca_dim2, MAX(1,n_ukca_cpnt))
REAL :: local_ukca_dry(ukca_dim1, ukca_dim2, MAX(1,n_ukca_mode))
REAL :: local_ukca_wet(ukca_dim1, ukca_dim2, MAX(1,n_ukca_mode))
REAL :: local_ukca_rho(ukca_dim1, ukca_dim2, MAX(1,n_ukca_mode))
REAL :: local_ukca_vol(ukca_dim1, ukca_dim2, MAX(1,n_ukca_mode))
REAL :: local_ukca_wtv(ukca_dim1, ukca_dim2, MAX(1,n_ukca_mode))
REAL :: local_ukca_nbr(ukca_dim1, ukca_dim2, MAX(1,n_ukca_mode))

INTEGER ::                                                              &
  list_daylight_points(row_length*rows,n_swcall),                       &
  list_daylight_points_start(row_length*rows),                          &
  diag_row_list(row_length*rows,MAX(n_swcall,n_lwcall)),                &
!   List of row indices of points where diagnostics are calculated.
  diag_col_list(row_length*rows,MAX(n_swcall,n_lwcall))
!   List of column indices of points where diagnostics are calculated.

LOGICAL ::                                                              &
  switch(row_length,rows,MAX(n_swcall,n_lwcall))

REAL ::                                                                 &
  sol_azimuth(row_length,rows),                                         &
!   Solar azimuth angles (radians clockwise from grid north)
  cosz_beg(row_length,rows), cosz_end(row_length,rows),                 &
!   cos(zenith_angle) at beginning and end of timestep
  seconds_since_midnight,                                               &
  sindec_obs,                                                           &
!   sin(observer declination)
  eqt_obs,                                                              &
!   The equation of time towards a distant observer
  obs_fraction(row_length,rows),                                        &
!   Fraction of the timestep that the observer is above the horizon
  obs_azimuth(row_length,rows),                                         &
!   Azimuth angle of the observer (radians clockwise from grid north)
  obs_solid_angle(row_length,rows), scale_area, scale_area_to_1au,      &
!   Solid angle subtended by grid-box for an observer at 1 AU
  trans_solid_angle(row_length,rows),                                   &
!   Solid angle subtended by grid-box for a transit observer at 1 AU
  dir_flux_to_trans(row_length,rows),                                   &
!   Conversion factor from direct flux to flux at transit observer
  cosz_obs(row_length,rows)
!   cos(zenith angle of observer)

REAL, SAVE :: sindec
!   sin(solar declination)
REAL, SAVE :: scs
!   solar constant scaling factor
REAL, SAVE :: eq_time
!   The Equation of time

LOGICAL :: l_planet_obs
!   Calculate solid angles and phase angles for distant observer

INTEGER ::                                                              &
  first_point_dust_a,                                                   &
  first_point_dust_b,                                                   &
  first_point_biogenic,                                                 &
  first_point_sulpc,                                                    &
  first_point_soot,                                                     &
  first_point_biomass,                                                  &
  first_point_ocff,                                                     &
  first_point_arcl,                                                     &
  first_point_nitrate,                                                  &
  first_point_ukca,                                                     &
  n_rad_layers

LOGICAL :: l_complete_north
!   Flag to complete field on the northern polar row
LOGICAL :: l_complete_south
!   Flag to complete field on the southern polar row
LOGICAL :: l_complete_deg
!   Flag to complete field because of degradation

LOGICAL :: l_co2_3d
!   Controls use of 3D co2 field

LOGICAL :: l_t_incr_sw
!   Flag for SW temperature increments

LOGICAL :: l_cosp
!   Local flag for COSP

! SW diagnostics not on stash flag
REAL ::                                                                 &
  itoasw(row_length,rows),                                              &
  surfsw(row_length,rows)

REAL ::                                                                 &
  flux_below_690nm_surf(row_length,rows)

! Direct Photosynthetically Active Radiation diagnostic for all tsteps
LOGICAL :: l_direct_par
REAL :: flxdirparsurf(row_length, rows)


! Radiation fields 2. LW

REAL ::                                                                 &
  lwsea(row_length, rows),                                              &
  olr(row_length, rows),                                                &
  surflw(row_length,rows),                                              &
  top_absorption(row_length, rows)
!   needed by energy correction
REAL ::                                                                 &
  net_atm_flux (row_length, rows)

LOGICAL :: l_t_incr_lw
!   Flag for LW temperature increments

! needed for 8A boundary layer
REAL ::                                                                 &
  frac_tile_alb(land_field, ntype),                                     &
  t_rad_surf(row_length, rows),                                         &
!   Effective radiative temperature over whole grid-box
  t_rad_land(row_length, rows),                                         &
!   Effective radiative temperature over land portion of grid-box
  t_rad_sice(row_length, rows),                                         &
!   Effective radiative temperature over sea ice
  t_rad_solid(row_length, rows),                                        &
!   Radiative temperature over solid surface (now deprecated)
  tstar_tile_tmp(land_field),                                           &
!   Temporary copy of tiled surface temperature to preserve existing
!   scientific behaviour within reorganized code
  surf_emission(row_length, rows),                                      &
!   LW emitted from the surface to calculate dOLR divided by sbcon
  tile_frac(land_field,ntiles),                                         &
  emis_tiles(ntype),                                                    &
                                   ! Emissivity for each surface type
  emis_here,                                                            &
!   Emissivity of the current tile in the current grid-box
  emis_nosnow,                                                          &
!   Emissivity of the current tile in the current grid-box
!   ignoring snow
  emis_land(row_length, rows)
                                   ! JULES mean land emissivity (zero
                                   ! for ocean only grid-boxes).
INTEGER ::                                                              &
  land_index_i(land_field),                                             &
  land_index_j(land_field),                                             &
  ssi_index_i(ssi_pts),                                                 &
  ssi_index_j(ssi_pts),                                                 &
  type_pts(ntype),                                                      &
!   Number of points contining each functional category
  type_index(land_field,ntype),                                         &
!   Index over functional types for aggregating albedo and emissivity
  tile_pts(ntype),                                                      &
!   Number of points contining tiles of the current category
  tile_index(land_field,ntiles)
!   Index over tiles

! Local variables required for multiple time-stepping
LOGICAL :: l_call_swrad
LOGICAL :: l_call_lwrad

! Local variables required for advancing the model. These are required
! for the improved timestepping algorithm.
REAL, ALLOCATABLE, SAVE :: sw_incs_local(:,:,:,:)
REAL, ALLOCATABLE, SAVE :: lw_incs_local(:,:,:,:)
REAL, ALLOCATABLE, SAVE :: netsw_local(:,:,:)
REAL, ALLOCATABLE, SAVE :: swsea_local(:,:,:)
REAL, ALLOCATABLE, SAVE :: lwsea_local(:,:,:)
REAL, ALLOCATABLE, SAVE :: top_abs_sw(:,:,:)
REAL, ALLOCATABLE, SAVE :: top_abs_lw(:,:,:)
REAL, ALLOCATABLE, SAVE :: olr_local(:,:,:)
REAL, ALLOCATABLE, SAVE :: lw_down_local(:,:,:)
REAL, ALLOCATABLE, SAVE :: flux_b690nm_local(:,:,:)
REAL, ALLOCATABLE, SAVE :: surf_down_sw_local(:,:,:,:)
REAL, ALLOCATABLE :: open_sea_albedo(:,:,:,:,:)

! Looping variables for extra calls to radiation
INTEGER :: j_sw
INTEGER :: j_lw

! Offset in STASH item number for different sets of diagnostics
INTEGER :: i_off


! Local variables relating to the radiance code

! Dynamically defined dimensions used by the radiance code
INTEGER :: nd_field_flux_diag
INTEGER :: nd_field_rad_diag


! Local Variables required to calculate radiative forcings.

! Temporary fields holding the mass mixing ratios from the grgas_field array
REAL, ALLOCATABLE :: n2o_mmr_tmp(:,:,:)
REAL, ALLOCATABLE :: ch4_mmr_tmp(:,:,:)
REAL, ALLOCATABLE :: cfc11_mmr_tmp(:,:,:)
REAL, ALLOCATABLE :: cfc12_mmr_tmp(:,:,:)
REAL, ALLOCATABLE :: cfc113_mmr_tmp(:,:,:)
REAL, ALLOCATABLE :: hcfc22_mmr_tmp(:,:,:)

! Variables holding the mass mixing ratios of either the
! reference call or the time advancing call.
REAL ::                                                                 &
  j_co2_mmr,                                                            &
  j_o2_mmr,                                                             &
  j_n2o_mmr,                                                            &
  j_ch4_mmr,                                                            &
  j_so2_mmr,                                                            &
  j_cfc11_mmr,                                                          &
  j_cfc12_mmr,                                                          &
  j_c113_mmr,                                                           &
  j_c114_mmr,                                                           &
  j_hcfc22_mmr,                                                         &
  j_hfc125_mmr,                                                         &
  j_hfc134_mmr

! The variable j_ozone contains the reference ozone field for the
! diagnostic call and the 'normal' ozone field for the time-advancing
! call.
REAL :: j_ozone(row_length, rows, ozone_levels)

! Same as above for logicals defining the use of aerosols
LOGICAL ::                                                              &
  j_l_sulpc_so2,                                                        &
  j_l_use_sulpc_direct,                                                 &
  j_l_use_seasalt_direct,                                               &
  j_l_soot,                                                             &
  j_l_use_soot_direct,                                                  &
  j_l_biomass,                                                          &
  j_l_use_bmass_direct,                                                 &
  j_l_ocff,                                                             &
  j_l_use_ocff_direct,                                                  &
  j_l_nitrate,                                                          &
  j_l_use_nitrate_direct,                                               &
  j_l_murk_rad,                                                         &
  j_l_dust,                                                             &
  j_l_use_dust,                                                         &
  j_l_use_biogenic,                                                     &
  j_l_climat_aerosol,                                                   &
  j_l_use_arcl(npd_arcl_species),                                       &
  j_l_use_ukca_radaer,                                                  &
  j_l_use_glomap_clim_radaer,                                           &
  j_l_use_easyaerosol

INTEGER :: j_n_arcl_species

!Automatic segment size tuning
TYPE(autotune_type), ALLOCATABLE, SAVE :: sw_autotune_state
TYPE(autotune_type), ALLOCATABLE, SAVE :: lw_autotune_state
INTEGER :: segment_size

! Structures for the core radiation code interface
TYPE(StrDim) :: dimen

! STASH section numbers
INTEGER, PARAMETER :: sw_sect=1
INTEGER, PARAMETER :: lw_sect=2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
eps  = EPSILON(1.0)

! ---------------------------------------------------------------------
! Section: CHECKS. Check that flags and timesteps are in agreement
! ---------------------------------------------------------------------
IF (l_timestep) THEN
  IF ( (a_sw_radstep_prog < a_sw_radstep_diag) .OR.                     &
       (a_lw_radstep_prog < a_lw_radstep_diag) ) THEN
    cmessage =                                                          &
      "If improved timestepping is required the diagnostic " //         &
      "timestep has to be smaller than the prognostic "   //            &
      "timestep."
    error_code=20
    GO TO 9999
  END IF
  IF (MOD(a_sw_radstep_prog, a_sw_radstep_diag) /= 0) THEN
    cmessage =                                                          &
      "The prognostic(slow) timestep needs to be " //                   &
      "a multiple of the diagnostic (fast) timestep."
    error_code=20
    GO TO 9999
  END IF
END IF

! ----------------------------------------------------------------------
! Section INI. Initialisation of stash output variables.
! ----------------------------------------------------------------------

l_co2_3d = l_co2_interactive

! Map the stashflags to the internal controlling logicals for
! diagnostics that are available on all timesteps:
SELECT CASE (model_type)

CASE (mt_single_column)
  l_t_incr_sw   = .TRUE.
  l_t_incr_lw   = .TRUE.
  l_direct_par  = .FALSE.
  l_surfsw_cor  = l_rad_szacor
  l_toasw_cor   = l_rad_szacor
  l_surfdir_cor = l_rad_szacor
  l_surfdif_cor = l_rad_szacor
  l_swdn_sice_wt_cat = .FALSE.
  l_swdn_sice_wt     = .FALSE.
  l_swup_sice_wt_cat = .FALSE.
  l_swup_sice_wt     = .FALSE.
  l_alb_sice_wt_cat  = .FALSE.
  l_alb_sice_wt      = .FALSE.
  l_lwdn_sice_wt_cat = .FALSE.
  l_lwdn_sice_wt     = .FALSE.
  l_planet_obs       = .FALSE.

CASE DEFAULT
  l_t_incr_sw   = (sf_calc(181,1) .OR. sf_calc(161,1) .OR. sf_calc(232,1))
  l_t_incr_lw   = (sf_calc(181,2) .OR. sf_calc(161,2) .OR. sf_calc(232,2))
  l_direct_par  = sf_calc(291,1)
  l_surfsw_cor  = sf_calc(202,1)
  l_toasw_cor   = sf_calc(205,1)
  l_surfdir_cor = sf_calc(215,1)
  l_surfdif_cor = sf_calc(216,1)
  l_swdn_sice_wt_cat = sf_calc(500,1)
  l_swdn_sice_wt     = sf_calc(501,1)
  l_swup_sice_wt_cat = sf_calc(502,1)
  l_swup_sice_wt     = sf_calc(503,1)
  l_alb_sice_wt_cat  = sf_calc(504,1)
  l_alb_sice_wt      = sf_calc(505,1)
  l_lwdn_sice_wt_cat = sf_calc(500,2)
  l_lwdn_sice_wt     = sf_calc(501,2)
  l_planet_obs = sf_calc(513,1) .OR. sf_calc(514,1) .OR. &
                 sf_calc(713,1) .OR. sf_calc(714,1) .OR. &
                 sf_calc(513,2) .OR. sf_calc(514,2) .OR. &
                 sf_calc(713,2) .OR. sf_calc(714,2)

END SELECT ! model_type


! Set up the logicals for diagnostics that are only available on
! radiation timesteps, contained in the structures SW_diag and LW_diag.
CALL init_diag_logic(sw_diag(1), sw_spectrum(1))
CALL init_diag_logic(lw_diag(1), lw_spectrum(1))

IF (l_rad_step_prog) THEN

  i_off=0

  ! The LW logicals need to be set before the SW ones since
  ! the latter ones depends on the LW ones.

  CALL set_lwdiag_logic(sf_calc, nitems, nsects, 1, i_off)
  CALL set_swdiag_logic(sf_calc, nitems, nsects, 1, i_off)

END IF

IF (l_radiance) THEN

  DO j_lw = 2, n_lwcall

    i_off=90+(diagnostic_offset/20)*(j_lw-1)

    CALL init_diag_logic(lw_diag(j_lw), lw_spectrum(j_lw))
    CALL set_lwdiag_logic(sf_calc, nitems, nsects, j_lw, i_off)

  END DO
  DO j_sw = 2, n_swcall

    i_off=90+(diagnostic_offset/20)*(j_sw-1)

    CALL init_diag_logic(sw_diag(j_sw), sw_spectrum(j_sw))
    CALL set_swdiag_logic(sf_calc, nitems, nsects, j_sw, i_off)

  END DO
END IF

IF (l_forcing .OR. l_timestep) THEN

  IF (l_timestep) i_off=0
  IF (l_forcing)  i_off=diagnostic_offset

  CALL init_diag_logic(sw_diag(2), sw_spectrum(2))
  CALL init_diag_logic(lw_diag(2), lw_spectrum(2))

  IF (l_rad_step_diag) THEN

    CALL set_lwdiag_logic(sf_calc, nitems, nsects, 2, i_off)
    CALL set_swdiag_logic(sf_calc, nitems, nsects, 2, i_off)

  END IF
END IF


IF (l_planet_grey_surface) THEN

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, flandg, land_sea_mask, land0p5, ice_fract )   &
!$OMP& PRIVATE( i, j )
  DO j = 1, rows
    DO i = 1, row_length
      flandg(i,j)=0.0
      land_sea_mask(i,j)=.FALSE.
      land0p5(i,j)=.FALSE.
      ice_fract(i,j)=0.0
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, land_sea_mask, land_sea_mask_in, land0p5,     &
!$OMP&         land0p5_in, ice_fract, ice_fract_in, flandg )                   &
!$OMP& PRIVATE( i, j )
  DO j = 1, rows
    DO i = 1, row_length
      land_sea_mask(i,j)=land_sea_mask_in(i,j)
      land0p5(i,j)=land0p5_in(i,j)
      ice_fract(i,j)=ice_fract_in(i,j)
      flandg(i,j)=0.0
    END DO
  END DO
!$OMP END PARALLEL DO

  ! --------------------------------------------------------------------
  ! Set tile_pts, type_index and tile_index.
  ! --------------------------------------------------------------------
  ! Type index over functional types.
  CALL tilepts(land_field,frac,type_pts,type_index)
  IF (l_aggregate) THEN
    tile_pts(1) = land_field
    DO l = 1, land_field
      tile_index(l,1) = l
      tile_frac(l,1) = 1.0
    END DO
  ELSE
    DO n = 1, ntype
      tile_pts(n)=type_pts(n)
      DO j = 1, type_pts(n)
        ! Here tiles match types.
        tile_index(j,n)=type_index(j,n)
        l = tile_index(j,n)
        tile_frac(l,n) = frac(l,n)
      END DO
    END DO
  END IF

  ! Set land_index.
  DO l = 1, land_field
    j = (land_index(l)-1)/row_length + 1
    land_index_i(l) = land_index(l) - (j-1)*row_length
    land_index_j(l) = j
  END DO

  ! Set ssi_index.
  DO l = 1, ssi_pts
    j = (ssi_index(l)-1)/row_length + 1
    ssi_index_i(l) = ssi_index(l) - (j-1)*row_length
    ssi_index_j(l) = j
  END DO

  ! --------------------------------------------------------------------
  ! Set global land fraction
  ! --------------------------------------------------------------------
  DO l = 1, land_field
    i = land_index_i(l)
    j = land_index_j(l)
    flandg(i,j)=fland(l)
  END DO
END IF

! set CO2 array dimensions for passing to SW & LW schemes
! If used:
! CO2_DIM1 = row_length*rows (calc'd from CO2_DIM_LEN & CO2_DIM_ROW)
! CO2_DIM2 = model_levels (already passed in from CO2_DIM_LEV)
co2_dim1 = co2_dim_len * co2_dim_row

! Set the number of layers seen in the radiation code.
! This may optionally be 1 greater than the number used in the rest of
! the model to avoid spurious effects resulting from the upper boundary.
IF (l_extra_top) THEN
  n_rad_layers=model_levels+1
ELSE
  n_rad_layers=model_levels
END IF

! ----------------------------------------------------------------------
! Section RAD.0 Set up required fields and constants.
! -----------------------------------------------------------------------

IF (error_code  ==  0) THEN
  IF (ltimer) CALL timer ('AP1R SW Rad  ',5)

  ! Calculate number of seconds since midnight to the beginning of
  ! the timetsep.
  seconds_since_midnight = REAL( previous_time(4) * 3600                &
           + previous_time(5) * 60  + previous_time(6))

  ! Take a temporary copy of the tiled surface temperature if using
  ! aggregated surface properties. Subsequently the actual tiled
  ! temperature will be over-written with by t_surf to enable the
  ! existing scientific behaviour to be replicated by the reorganized
  ! code that works with tstar_tile. (Currently the temperatures may
  ! differ after reconfiguration, though logically they should be
  ! equal, and will be equal once they have been through the
  ! surface scheme.) tstar_tile is restored from the temporary copy
  ! at the end of the routine. Once an acceptable method for making
  ! these temperatures consistent in reconfiguration has been
  ! implemented this temporary copy can be removed.
  IF (l_aggregate) tstar_tile_tmp=tstar_tile(:,1)

  ! Set radiation fields to zero
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, flux_below_690nm_surf, sw_net_land,           &
!$OMP&         sw_net_sice )                                                   &
!$OMP& PRIVATE( i, j )
  DO j = 1, rows
    DO i = 1, row_length
      flux_below_690nm_surf(i,j) = 0.0
      sw_net_land(i,j) = 0.0
      sw_net_sice(i,j) = 0.0
    END DO
  END DO
!$OMP END PARALLEL DO

  IF (l_orog) THEN
    ALLOCATE(f_orog(row_length,rows))
    f_orog=0.0
  ELSE
    ALLOCATE(f_orog(1,1))
  END IF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j, k, l, n )
!$OMP DO SCHEDULE(STATIC)
  DO n = 1, ntiles
    DO l = 1, land_field
      sw_tile(l,n) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO l = 1, ssi_pts
    sw_sea(l) = 0.0
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO n = 1, nice_use
    DO l = 1, ssi_pts
      sw_sicat(l,n) = 0.0
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        t_latest(i,j,k)= t_n(i,j,k)
        q_latest(i,j,k)= q_n(i,j,k)
        qcl_latest(i,j,k)=qcl_n(i,j,k)
        cf_latest(i,j,k)=cf_n(i,j,k)
        cfl_latest(i,j,k)=cfl_n(i,j,k)
        cff_latest(i,j,k)=cff_n(i,j,k)
        zeros(i,j,k)=0.0
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF (l_rad_step_diag .OR. l_rad_step_prog) THEN

    ! Code required to initialize and update snow soot content
    IF ( l_snow_albedo .OR. l_embedded_snow ) THEN
      DO j = 1, rows
        DO i = 1, row_length
          soot(i,j) = 0.0
        END DO
      END DO
    END IF

    ! Code for the mineral dust aerosol scheme. We copy
    ! dust aerosol into local arrays if the direct effect is
    ! switched on.
    IF (l_dust .AND. l_rad_step_prog) THEN
      IF (dust_dim1  ==  rows*row_length .AND.                          &
                         dust_dim2  ==  model_levels) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE ( i, j, k, ij )
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              dust_1(ij,k) = dust_div1(i,j,k)
              dust_2(ij,k) = dust_div2(i,j,k)
              IF (.NOT. l_twobin_dust) THEN
                dust_3(ij,k) = dust_div3(i,j,k)
                dust_4(ij,k) = dust_div4(i,j,k)
                dust_5(ij,k) = dust_div5(i,j,k)
                dust_6(ij,k) = dust_div6(i,j,k)
              END IF
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
      ELSE
        cmessage = 'DUST_DIM INCONSISTENT WITH L_DUST'
        error_code=1
        GO TO 9999
      END IF
    END IF

    ! Code for the biogenic aerosol. Copy into local arrays if
    ! this aerosol was requested.
    IF (l_use_biogenic .AND. l_rad_step_prog) THEN
      IF (biogenic_dim1  ==  rows*row_length .AND.                             &
          biogenic_dim2  ==  model_levels) THEN
!$OMP  PARALLEL DO DEFAULT(NONE)                                               &
!$OMP& SHARED( model_levels, rows, row_length, local_biogenic, biogenic )      &
!$OMP& PRIVATE( i, j, k, ij )
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              local_biogenic(ij,k) = biogenic(i,j,k)
            END DO
          END DO
        END DO
!$OMP END PARALLEL DO
      ELSE
        cmessage = 'BIOGENIC_DIM INCONSISTENT WITH L_USE_BIOGENIC'
        error_code=1
        GO TO 9999
      END IF
    END IF

    ! Code for the Sulphur Cycle. We multiply by 4.125 to convert from
    ! mass mixing ratio of sulphur atoms to mass mixing ratio of
    ! ammonium sulphate.
    IF ((l_sulpc_so2 .AND. l_rad_step_prog) .OR.                        &
      l_use_sulpc_indirect_sw .OR. l_use_sulpc_indirect_lw) THEN
      IF (sulp_dim1  ==  rows*row_length .AND.                          &
                           sulp_dim2  ==  model_levels) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              accum_sulphate(ij,k)=so4_accu(i,j,k)*4.125
              aitken_sulphate(ij,k)=so4_aitken(i,j,k)*4.125
              diss_sulphate(ij,k)=so4_diss(i,j,k)*4.125
            END DO
          END DO
        END DO
      ELSE
        cmessage = 'SULP_DIM INCONSISTENT WITH L_SULPC'
        error_code = 1
        GO TO 9999
      END IF
    END IF

    ! Code for the Soot Scheme. As for the Sulphur Cycle (above), we
    ! copy soot into local arrays if the direct effect is switched on,
    ! but no multiplication is required.
    IF (l_soot .AND. l_rad_step_prog) THEN
      IF (soot_dim1  ==  rows*row_length .AND.                          &
                         soot_dim2  ==  model_levels) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              fresh_soot(ij,k) = soot_new(i,j,k)
              aged_soot(ij,k) = soot_agd(i,j,k)
            END DO
          END DO
        END DO
      ELSE
        cmessage = 'SOOT_DIM INCONSISTENT WITH L_SOOT'
        error_code=1
        GO TO 9999
      END IF
    END IF

    ! Code for the biomass aerosol scheme. As for soot, we copy biomass
    ! smoke aerosol into local arrays if the direct/indirect effect is
    ! switched on.
    IF ((l_biomass .AND. l_rad_step_prog)                               &
                             .OR. l_use_bmass_indirect) THEN
      IF (bmass_dim1  ==  rows*row_length .AND.                         &
                         bmass_dim2  ==  model_levels) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              fresh_bmass(ij,k) = bmass_new(i,j,k)
              aged_bmass(ij,k)  = bmass_agd(i,j,k)
              cloud_bmass(ij,k) = bmass_cld(i,j,k)
            END DO
          END DO
        END DO
      ELSE
        cmessage = 'BMASS_DIM INCONSISTENT WITH L_BIOMASS'
        error_code=1
        GO TO 9999
      END IF
    END IF

    ! Code for the fossil-fuel organic carbon aerosol scheme. As for
    ! biomass, we copy fossil-fuel organic carbon aerosol into local
    ! arrays if the direct/indirect effect is switched on.
    IF ((l_ocff .AND. l_rad_step_prog) .OR. l_use_ocff_indirect) THEN
      IF (ocff_dim1  ==  rows*row_length .AND.                          &
          ocff_dim2  ==  model_levels) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              fresh_ocff(ij,k) = ocff_new(i,j,k)
              aged_ocff(ij,k)  = ocff_agd(i,j,k)
              cloud_ocff(ij,k) = ocff_cld(i,j,k)
            END DO
          END DO
        END DO
      ELSE
        cmessage = 'OCFF_DIM INCONSISTENT WITH L_USE_OCFF'
        error_code=1
        GO TO 9999
      END IF
    END IF

    ! Code for the NWP aerosol climatology. As above, we copy the input
    ! into local arrays if the climatology was requested.
    IF (n_arcl_species > 0 .AND. l_rad_step_prog) THEN

      DO l = 1, npd_arcl_compnts

        IF (i_arcl_compnts(l) /= -1) THEN
          IF (arcl_dim1 == rows * row_length .AND.                      &
              arcl_dim2 == model_levels) THEN

            ! Sulphur mass-mixing ratios have to be converted
            ! into ammonium sulphate.
            IF ((l == ip_arcl_sulp_ac) .OR. (l == ip_arcl_sulp_ak)) THEN
              arcl_multf = 4.125
            ELSE
              arcl_multf = 1.0
            END IF

!$OMP  PARALLEL DO DEFAULT(NONE)                                               &
!$OMP& SHARED( model_levels, rows, row_length, local_arcl, arcl_multf,         &
!$OMP&         i_arcl_compnts, l, arcl )                                       &
!$OMP& PRIVATE( i, j, k, ij )
            DO k = 1, model_levels
              DO j = 1, rows
                DO i = 1, row_length
                  ij=i+(j-1)*row_length
                  local_arcl(ij,k,i_arcl_compnts(l)) =                  &
                        arcl_multf *                                    &
                        arcl(i,j,k,i_arcl_compnts(l))
                END DO ! i
              END DO ! j
            END DO ! k
!$OMP END PARALLEL DO
          ELSE
            cmessage = 'ARCL_DIM IS INCONSISTENT WITH L_USE_ARCL'
            error_code=1
            GO TO 9999
          END IF
        END IF

      END DO
    END IF ! n_arcl_species

    ! Code for the nitrate aerosol scheme. We copy
    ! nitrate aerosol into local arrays if the direct/indirect effect is
    ! switched on and we multiply by 5.714 to convert from
    ! mass mixing ratio of nitrogen atoms to mass mixing ratio of
    ! ammonium nitrate.
    IF ((l_nitrate .AND. l_rad_step_prog)                               &
                               .OR. l_use_nitrate_indirect) THEN
      IF (nitrate_dim1  ==  rows*row_length .AND.                       &
          nitrate_dim2  ==  model_levels) THEN
        DO k = 1, model_levels
          DO j = 1, rows
            DO i = 1, row_length
              ij=i+(j-1)*row_length
              accum_nitrate(ij,k) = nitr_acc(i,j,k)*5.714
              diss_nitrate(ij,k)  = nitr_diss(i,j,k)*5.714
            END DO
          END DO
        END DO
      ELSE
        cmessage = 'NITRATE_DIM INCONSISTENT WITH L_NITRATE'
        error_code=1
        GO TO 9999
      END IF
    END IF

    ! Code for UKCA_RADAER. As for all other aerosols, we copy
    ! UKCA aerosols into local arrays if they are switched on.
    IF (l_ukca_radaer .OR. l_glomap_clim_radaer) THEN

      IF (ukca_dim1 == rows*row_length .AND.                            &
          ukca_dim2 == model_levels) THEN

        ! Component-related variables.
        DO l = 1, n_ukca_cpnt
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij = i + (j-1)*row_length
                local_ukca_mmr(ij, k, l) =                              &
                      ukca_radaer%mix_ratio(i,j,k,l)
                local_ukca_cvl(ij, k, l) =                              &
                      ukca_radaer%comp_vol(i,j,k,l)
              END DO ! i
            END DO ! j
          END DO ! k
        END DO ! l

        ! Mode-related variables.
        DO l = 1, n_ukca_mode
          DO k = 1, model_levels
            DO j = 1, rows
              DO i = 1, row_length
                ij = i + (j-1)*row_length
                local_ukca_dry(ij, k, l) =                              &
                      ukca_radaer%dry_diam(i,j,k,l)
                local_ukca_wet(ij, k, l) =                              &
                      ukca_radaer%wet_diam(i,j,k,l)
                local_ukca_rho(ij, k, l) =                              &
                      ukca_radaer%modal_rho(i,j,k,l)
                local_ukca_vol(ij, k, l) =                              &
                      ukca_radaer%modal_vol(i,j,k,l)
                local_ukca_wtv(ij, k, l) =                              &
                      ukca_radaer%modal_wtv(i,j,k,l)
                local_ukca_nbr(ij, k, l) =                              &
                      ukca_radaer%modal_nbr(i,j,k,l)
              END DO ! i
            END DO ! j
          END DO ! k
        END DO ! l

      ELSE
        cmessage = 'UKCA_DIM INCONSISTENT WITH L_UKCA_RADAER'
        error_code=1
        GO TO 9999
      END IF
    END IF

    ! Configure cloud ice total amounts. For the case where there is a single
    ! ice prognostic, l_mcr_qcf2 is .FALSE. and qcf_total should be qcf_n.
    ! With two ice prognostics, qcf_total is the sum of both species.

    IF ( l_mcr_qcf2 ) THEN
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, model_levels, qcf_n, qcf2_n, qcf_total )      &
!$OMP& PRIVATE( i, j, k )
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            qcf_total(i,j,k) = qcf_n(i,j,k) + qcf2_n(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

    ELSE ! not l_mcr_qcf2
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, model_levels, qcf_n, qcf_total )              &
!$OMP& PRIVATE( i, j, k )
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            qcf_total(i,j,k) = qcf_n(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

    END IF ! l_mcr_qcf2

    ! Set layer masses, heat capacities and boundary temperatures
    CALL set_thermodynamic(row_length, rows, off_x, off_y, t_n, rho_r2,        &
      p_layer_centres, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n,           &
      l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                                   &
      p_layer_boundaries, t_layer_boundaries, p_extra_layer, t_extra_layer,    &
      d_mass, density, layer_heat_capacity)

  END IF ! l_rad_step_diag .OR. l_rad_step_prog


  ! -----------------------------------------------------------------------
  ! Section Rad.0.1 Allocate SW variables and diagnostics
  ! -----------------------------------------------------------------------

  ! work out the maximum number of sw bands in all calls:
  max_n_swbands = 0
  DO j_sw = 1, n_swcall
    IF (sw_spectrum(j_sw)%basic%n_band > max_n_swbands) THEN
      max_n_swbands = sw_spectrum(j_sw)%basic%n_band
    END IF
  END DO
  ALLOCATE(open_sea_albedo(row_length, rows, 2, max_n_swbands, n_swcall))

  IF (l_rad_step_prog) THEN

    CALL allocate_diag(sw_diag(1), sw_spectrum(1),                      &
      row_length, rows, model_levels, cloud_levels, ntiles)

    IF (l_orog) ALLOCATE(orog_corr(row_length,rows))
    ALLOCATE(cos_zen_rts(row_length, rows))
    ALLOCATE(sol_azm_rts(row_length, rows))
    ALLOCATE(day_frac_rts(row_length, rows))

    IF (ANY(sw_control(1:n_swcall)%l_spherical_solar)) THEN
      ALLOCATE(cos_zen_sph_rts(row_length, rows, 0:n_rad_layers+1))
      ALLOCATE(day_frac_sph_rts(row_length, rows, 0:n_rad_layers+1))
    ELSE
      ALLOCATE(cos_zen_sph_rts(row_length, rows, 0:1))
      ALLOCATE(day_frac_sph_rts(row_length, rows, 0:1))
    END IF

    ALLOCATE(surf_down_sw_rts(row_length, rows, 4))
    ALLOCATE(alb_tile(land_field, ntiles, 4))

    IF (l_timestep) THEN
      ALLOCATE(sw_incs_local(row_length, rows, 0:model_levels+1,n_swcall))
      ALLOCATE(surf_down_sw_local(row_length,rows,4,n_swcall))
      ALLOCATE(netsw_local(row_length, rows, n_swcall))
      ALLOCATE(swsea_local(row_length, rows ,n_swcall))
      ALLOCATE(top_abs_sw(row_length, rows, n_swcall))
      ALLOCATE(flux_b690nm_local(row_length,rows,n_swcall))
    END IF

    IF (l_radiance) THEN
      DO j_sw = 2, n_swcall
        CALL allocate_diag(sw_diag(j_sw), sw_spectrum(j_sw),            &
          row_length, rows, model_levels, cloud_levels, ntiles)
      END DO
    END IF

  END IF

  IF (l_rad_step_diag) THEN
    IF (l_timestep .OR. l_forcing) THEN
      CALL allocate_diag(sw_diag(2), sw_spectrum(2),                    &
        row_length, rows, model_levels, cloud_levels, ntiles)
    END IF
  END IF

  ! -----------------------------------------------------------------------
  ! Section Rad.0.2 Astronomy
  ! -----------------------------------------------------------------------
  IF (l_rad_step_prog) THEN
    ! Calculates sine of the solar declination and the scaling
    ! factor for solar intensity from the day number and year.
    CALL solpos (previous_time(7), previous_time(1), seconds_since_midnight,   &
         radiation_tstep_prog, l_planet_obs,                                   &
         eq_time, sindec, scs, sindec_obs, eqt_obs)

    ! Set solar zenith angles
    CALL solang(                                                               &
      sindec, seconds_since_midnight, radiation_tstep_prog, eq_time,           &
      true_latitude, true_longitude, row_length*rows,                          &
      day_frac_rts, cos_zen_rts, sol_azm_rts, cosz_beg, cosz_end)

    ! Set rounding-error size values to zero - the criterion depends
    ! on the frequency of full SW calculations because on the physics
    ! timesteps which are not SW timesteps a test has to be done to
    ! avoid using the unset data for such points.
    DO j = 1, rows
      DO i = 1, row_length
        IF ( cos_zen_rts(i,j)*day_frac_rts(i,j) <                              &
             (1.0e-10*timestep/radiation_tstep_prog) ) THEN
          cos_zen_rts(i,j) = 0.0
          day_frac_rts(i,j) = 0.0
        END IF
      END DO
    END DO

    ! Solar variability
    IF (sw_control(1)%l_solvar) THEN
      ! Apply solar spectral variation
      CALL solvar(previous_time, sw_spectrum(1), solcon_rts)
    ELSE
      ! Otherwise set solar constant from input namelist value
      solcon_rts=sc
    END IF

    ! Calculate orography correction:
    IF (l_orog .AND. .NOT. sw_control(1)%l_spherical_solar) THEN
       CALL solinc(row_length, rows, cos_zen_rts, sol_azm_rts,                 &
         cosz_beg, cosz_end)
      IF (l_rad_szacor) THEN
        WHERE (ABS(orog_corr-0.5) < eps)
          orog_corr=0.5+eps
        END WHERE
      END IF
    END IF

    IF (ANY(sw_control(1:n_swcall)%l_spherical_solar) .OR. l_planet_obs) THEN
      scale_area = r_layer_boundaries(1,1,n_rad_layers)**2 / planet_radius**2
      scale_area_to_1au = scale_area / astronomical_unit**2
    END IF

    IF (ANY(sw_control(1:n_swcall)%l_spherical_solar)) THEN
      CALL solang_sph(                                                         &
        sindec, seconds_since_midnight, radiation_tstep_prog, eq_time,         &
        true_latitude, true_longitude, row_length*rows, n_rad_layers,          &
        r_layer_centres, r_layer_boundaries, slope_aspect, slope_angle,        &
        horiz_ang, horiz_aspect, n_horiz_layer, n_horiz_ang, l_orog, l_skyview,&
        day_frac_sph_rts, cos_zen_sph_rts, sol_azm_rts)
      DO j = 1, rows
        DO i = 1, row_length
          k = n_rad_layers+1
          ! Solid angle subtended by the gridbox as seen from an observer
          ! viewing a transit at a distance of 1 AU. (Note cos_zen_sph_rts is
          ! negative for these points which are lit from below.)
          trans_solid_angle(i,j) = cell_area_surface(i,j) * scale_area_to_1au  &
            * (-cos_zen_sph_rts(i,j,k)) * day_frac_sph_rts(i,j,k)
          ! Conversion factor from direct flux to the flux seen from an
          ! observer viewing a transit at a distance of 1 AU.
          dir_flux_to_trans(i,j) = cell_area_surface(i,j) * scale_area         &
            * day_frac_sph_rts(i,j,k) / (scs*pi*stellar_radius**2)
        END DO
      END DO
    END IF

    IF (l_planet_obs) THEN
      ! Find the zenith angle of the observer for diagnostics
      CALL solang(                                                             &
        sindec_obs, seconds_since_midnight, radiation_tstep_prog, eqt_obs,     &
        true_latitude, true_longitude, row_length*rows,                        &
        obs_fraction, cosz_obs, obs_azimuth, cosz_beg, cosz_end )

      ! Find the solid angle subtended by the gridbox
      ! as seen from an observer at a distance of 1 AU.
      DO j = 1, rows
        DO i = 1, row_length
          obs_solid_angle(i,j) = cell_area_surface(i,j) * scale_area_to_1au    &
            * cosz_obs(i,j) * obs_fraction(i,j)
        END DO
      END DO
    END IF
  END IF


  ! ------------------------------------------------------------
  ! Section Rad.1 Short Wave Radiation Code
  ! ------------------------------------------------------------
  IF (l_rad_step_prog .OR. l_rad_step_diag) THEN
    DO j_sw = n_swcall, 1, -1
      ! Decide whether it is necessary to call radiation
      l_call_swrad=.FALSE.
      IF (j_sw==1 .AND. l_rad_step_prog) THEN
        l_call_swrad=.TRUE.
      ELSE IF (l_timestep .AND. j_sw == 2) THEN
        l_call_swrad=.TRUE.
      ELSE IF (l_forcing .AND. j_sw > 1) THEN
        IF (sw_diag(j_sw)%l_diag_call) THEN
          l_call_swrad=.TRUE.
        END IF
      END IF

      ! If call to SW Radiation required
      IF (l_call_swrad) THEN

        ! Set COSP flag only on prognostic radiation steps
        IF (l_cosp_in .AND. (j_sw == 1) .AND. l_rad_step_prog) THEN
          l_cosp=.TRUE.
        ELSE
          l_cosp=.FALSE.
        END IF

        ! For each SW call set the correct mass mixing ratios and switch the
        ! aerosols on/off.
        ! Note that in the case of the aerosols the reference state is
        ! 'no aerosols' and hence they are switched off in the diagnostic call.

        ! Since the "default" for each variable is to have it set in the
        ! diagnostic call to the same as in the prognostic, we set them all
        ! unconditionally to the prognostic & then re-set as needed.
        j_co2_mmr = co2_mmr
        j_o2_mmr  = o2mmr
        j_n2o_mmr = n2ommr
        j_ch4_mmr = ch4mmr
        j_so2_mmr = so2mmr
        j_ozone = ozone  ! Note this is an array assignment
        j_l_sulpc_so2          = l_sulpc_so2
        j_l_use_sulpc_direct   = l_use_sulpc_direct
        j_l_use_seasalt_direct = l_use_seasalt_direct
        j_l_soot               = l_soot
        j_l_use_soot_direct    = l_use_soot_direct
        j_l_biomass            = l_biomass
        j_l_use_bmass_direct   = l_use_bmass_direct
        j_l_ocff               = l_ocff
        j_l_use_ocff_direct    = l_use_ocff_direct
        j_l_nitrate            = l_nitrate
        j_l_use_nitrate_direct = l_use_nitrate_direct
        j_l_dust               = l_dust
        j_l_use_dust           = l_use_dust
        j_l_use_biogenic       = l_use_biogenic
        j_l_climat_aerosol     = l_climat_aerosol
        j_l_murk_rad           = l_murk_rad
        j_l_use_arcl           = l_use_arcl
        j_n_arcl_species       = n_arcl_species
        j_l_use_ukca_radaer    = l_ukca_radaer
        j_l_use_glomap_clim_radaer = l_glomap_clim_radaer
        frac_tile_alb          = frac
        solar_constant         = solcon_rts
        j_l_use_easyaerosol    = l_easyaerosol_sw

        ! If incremental timestepping required
        IF (l_timestep .AND. j_sw == 2) THEN

          ! Turn aerosols off for the "cloud only" call in order to minimise
          ! computational cost. Clear-sky and other diagnostics have already
          ! been turned off in set_swdiag_logic.
          IF (l_rad_perturb) THEN
            j_l_sulpc_so2          = .FALSE.
            j_l_use_sulpc_direct   = .FALSE.
            j_l_use_seasalt_direct = .FALSE.
            j_l_soot               = .FALSE.
            j_l_use_soot_direct    = .FALSE.
            j_l_biomass            = .FALSE.
            j_l_use_bmass_direct   = .FALSE.
            j_l_ocff               = .FALSE.
            j_l_use_ocff_direct    = .FALSE.
            j_l_nitrate            = .FALSE.
            j_l_use_nitrate_direct = .FALSE.
            j_l_dust               = .FALSE.
            j_l_use_dust           = .FALSE.
            j_l_use_biogenic       = .FALSE.
            j_l_climat_aerosol     = .FALSE.
            j_l_murk_rad           = .FALSE.
            j_l_use_arcl           = .FALSE.
            j_n_arcl_species       = 0
            j_l_use_ukca_radaer    = .FALSE.
            j_l_use_glomap_clim_radaer = .FALSE.
            j_l_use_easyaerosol    = .FALSE.
          END IF

        ! If diagnostic call required
        ELSE IF (l_forcing .AND. j_sw > 1) THEN
          ! Solar variability
          IF (sw_control(j_sw)%l_solvar) THEN
            ! Apply solar spectral variation
            CALL solvar(previous_time, sw_spectrum(j_sw), solar_constant)
          ELSE
            ! Otherwise set solar constant from input namelist value
            solar_constant=sc
          END IF

          ! Gases
          IF ( c2c_co2 ) j_co2_mmr = co2_mmr*co2_mmr_scl + co2_mmr_add
          IF ( c2c_o2 )  j_o2_mmr  = o2mmr*o2_mmr_scl + o2_mmr_add
          IF ( c2c_n2o ) THEN
            j_n2o_mmr = n2ommr*n2o_mmr_scl + n2o_mmr_add
            IF (ngrgas >= p_n2o) THEN
              ALLOCATE(n2o_mmr_tmp(row_length, rows, model_levels))
              n2o_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_n2o)
              grgas_field(:,:,:,p_n2o) = &
                grgas_field(:,:,:,p_n2o)*n2o_mmr_scl + n2o_mmr_add
            END IF
          END IF
          IF ( c2c_ch4 ) THEN
            j_ch4_mmr = ch4mmr*ch4_mmr_scl + ch4_mmr_add
            IF (ngrgas >= p_ch4) THEN
              ALLOCATE(ch4_mmr_tmp(row_length, rows, model_levels))
              ch4_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_ch4)
              grgas_field(:,:,:,p_ch4) = &
                grgas_field(:,:,:,p_ch4)*ch4_mmr_scl + ch4_mmr_add
            END IF
          END IF

          ! Ozone: Note that aerosol contains the MURK array. Care needs to be
          ! taken so that MURK has the correct dimensions before it arrives in
          ! RAD_CTL2.
          IF ( c2c_o3 ) j_ozone = aerosol

          ! Aerosols
          IF ( c2c_sulpc_d ) j_l_sulpc_so2          = .FALSE.
          IF ( c2c_sulpc_d ) j_l_use_sulpc_direct   = .FALSE.
          IF ( c2c_seas_d )  j_l_use_seasalt_direct = .FALSE.
          IF ( c2c_soot_d )  j_l_soot               = .FALSE.
          IF ( c2c_soot_d )  j_l_use_soot_direct    = .FALSE.
          IF ( c2c_bmb_d )   j_l_biomass            = .FALSE.
          IF ( c2c_bmb_d )   j_l_use_bmass_direct   = .FALSE.
          IF ( c2c_ocff_d )  j_l_ocff               = .FALSE.
          IF ( c2c_ocff_d )  j_l_use_ocff_direct    = .FALSE.
          IF ( c2c_nitr_d )  j_l_nitrate            = .FALSE.
          IF ( c2c_nitr_d )  j_l_use_nitrate_direct = .FALSE.
          IF ( c2c_dust_d )  j_l_dust               = .FALSE.
          IF ( c2c_dust_d )  j_l_use_dust           = .FALSE.
          IF ( c2c_biog_d )  j_l_use_biogenic       = .FALSE.
          IF ( c2c_ukca_d )  j_l_use_ukca_radaer    = .FALSE.
          IF ( c2c_easy_d )  j_l_use_easyaerosol    = .FALSE.
          IF ( c2c_ukca_d )  j_l_use_glomap_clim_radaer = .FALSE.

          ! Land surface
          IF (c2c_land_s) frac_tile_alb = frac_control
        END IF

        IF (.NOT. l_planet_grey_surface) THEN
          CALL surf_couple_radiation(                                   &
          ! Fluxes INTENT(IN)
            t_surf,                                                     &
          ! Misc INTENT(IN)
            ws_10m_sea, chloro_sea,                                     &
            sw_spectrum(j_sw)%basic%n_band,                             &
            max_n_swbands,                                              &
            sw_spectrum(j_sw)%basic%wavelength_short,                   &
            sw_spectrum(j_sw)%basic%wavelength_long,                    &
          ! Misc INTENT(OUT)
            sea_ice_albedo(1,1,1,j_sw),                                 &
          ! Fluxes INTENT(OUT)
            alb_tile, land_albedo(1,1,1,j_sw),                          &
          ! UM-only args: INTENT(IN)
            pond_frac_cat, pond_depth_cat,                              &
          ! (ancil_info mod)
            ntiles, land_field, land_index, type_pts, type_index,       &
            row_length, rows, ice_fract, ice_fract_cat, frac_tile_alb,  &
          ! (p_s_parms mod)
            cos_zen_rts, albobs_sw, albobs_vis, albobs_nir,&
            z0_tile, albsoil,                                           &
          ! (coastal mod)
            flandg, tstar_sice_cat,                                     &
          ! (prognostics mod)
            snow_depth, snow_depth_sea_cat, ice_thick_cat,              &
            lai, canht, rgrain,                                         &
            snow_tile, soot, tstar_tile, sd_orog_land,                  &
          ! UM-only args: INTENT(OUT)
            albobs_sc, open_sea_albedo(1,1,1,1,j_sw) )

          ! Write the albedo scalings to sw_diag, if requested:
          IF (SW_diag(j_sw)%l_vis_albedo_sc) THEN
            DO k =1, ntiles
              DO j = 1, rows
                DO i = 1, row_length
                  sw_diag(j_sw)%vis_albedo_sc(i,j,k)=albobs_sc(i,j,k,1)
                END DO
              END DO
            END DO
          END IF
          IF (SW_diag(j_sw)%l_nir_albedo_sc) THEN
            DO k =1, ntiles
              DO j = 1, rows
                DO i = 1, row_length
                  sw_diag(j_sw)%nir_albedo_sc(i,j,k)=albobs_sc(i,j,k,2)
                END DO
              END DO
            END DO
          END IF
        END IF ! l_planet_grey_surface

        CALL prelim_swrad(error_code,                                   &
          at_extremity, n_proc,                                         &
        ! Model Dimensions
          row_length, rows, n_rad_layers,                               &
        ! Model Switches
          l_rad_deg, sw_control(j_sw)%l_subsample,                      &
          sw_control(j_sw)%l_geostationary,                             &
          sw_control(j_sw)%l_spherical_solar,                           &
        ! Time stepping Information
          timestep_number,a_sw_radstep_prog,                            &
        ! Geometry
          sw_control(j_sw)%min_view_lon,                                &
          sw_control(j_sw)%max_view_lon,                                &
          sw_control(j_sw)%min_view_lat,                                &
          sw_control(j_sw)%max_view_lat,                                &
        ! Number of Call
          j_sw,                                                         &
        ! Other variables
          true_latitude,true_longitude,                                 &
          seconds_since_midnight,                                       &
          tot_daylight_points,daylight_points,                          &
          day_frac_rts, day_frac_sph_rts,                               &
          list_daylight_points(1,j_sw), rad_mask,                       &
          switch(1,1,j_sw), first_data_interp,                          &
          diag_row_list(1,j_sw), diag_col_list(1,j_sw) )

        ! Zero the calculated outputs: this will later simplify the
        ! treatment of the case when there are no lit points.
!$OMP  PARALLEL DEFAULT(NONE)                                                  &
!$OMP& SHARED( model_levels, rows, row_length, sw_incs, surf_down_sw_rts,      &
!$OMP&         netsw, swsea, l_extra_top, top_absorption, l_cosp,              &
!$OMP&         day_frac_rts, j_sw, cosp_gbx )                                  &
!$OMP& PRIVATE ( i, j, k )
!$OMP  DO SCHEDULE(STATIC)
        DO k = 0, model_levels+1
          DO j = 1, rows
            DO i = 1, row_length
              sw_incs(i,j,k)=0.0
            END DO
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            netsw (i,j) = 0.0
            swsea (i,j) = 0.0
            surf_down_sw_rts (i,j,:) = 0.0
          END DO
        END DO
!$OMP END DO NOWAIT
        IF (l_extra_top) THEN
!$OMP DO SCHEDULE(STATIC)
          DO j = 1, rows
            DO i = 1, row_length
              top_absorption(i,j) = 0.0
            END DO
          END DO
!$OMP END DO NOWAIT
        END IF

        ! Sunlit points for COSP
        IF (l_cosp) THEN
!$OMP DO SCHEDULE(STATIC)
          DO j = 1, rows
            DO i = 1, row_length
               IF (day_frac_rts(i,j) >  0.0)                                   &
                cosp_gbx%sunlit((j-1)*row_length + i) = 1.0
            END DO
          END DO
!$OMP END DO
        END IF
!$OMP END PARALLEL

        !Set up automatic segment tuning for SW radiation
        IF (l_autotune_segments .AND. .NOT. ALLOCATED(sw_autotune_state)) THEN
          ALLOCATE(sw_autotune_state)
          CALL autotune_init(                   &
            sw_autotune_state,                  &
            name           = 'SW_Radiation',    &
            tag            = 'SW-RAD',          &
            start_size     = a_sw_seg_size,     &
            calls_per_step = n_swcall)
        END IF

        !If autotuning is active, get a new trial segment size.
        segment_size = a_sw_seg_size
        IF (l_autotune_segments) THEN
          segment_size = autotune_get_trial_size(sw_autotune_state)
          CALL autotune_start_region(sw_autotune_state)
        END IF

        IF ( daylight_points  >   0 ) THEN

          CALL set_control(sw_control(j_sw), sw_spectrum(j_sw),                &
          ! Control flags
            j_l_climat_aerosol, j_l_use_sulpc_direct, j_l_use_soot_direct,     &
            j_l_use_biogenic, j_l_use_dust, j_l_use_bmass_direct,              &
            j_l_use_ocff_direct, j_l_use_nitrate_direct,                       &
            j_l_use_seasalt_direct, j_l_murk_rad, l_use_sulpc_indirect_sw,     &
            j_n_arcl_species, j_l_use_ukca_radaer, j_l_use_glomap_clim_radaer, &
            j_l_use_easyaerosol,                                               &
          ! Diagnostic options
            sw_diag(j_sw), l_cosp)


          ! Parallelise over segments.
!$OMP  PARALLEL DEFAULT(SHARED)                                              &
!$OMP& PRIVATE(i, j, lit_points,start_point,first_point,first_point_dust_a,  &
!$OMP& first_point_dust_b,                                                   &
!$OMP& first_point_sulpc,first_point_soot,first_point_biomass,               &
!$OMP& first_point_biogenic, first_point_ocff, first_point_arcl,             &
!$OMP& first_point_nitrate, dimen,                                           &
!$OMP& nd_field_flux_diag, nd_field_rad_diag, first_point_ukca,              &
!$OMP& meta_segments, segments, ipar,                                        &
!$OMP& ii, jj)

!$OMP SINGLE
          ! Determine the number of threads in this parallel region
          num_parallel_sw=1
!$ num_parallel_sw = omp_get_num_threads()
!$OMP END SINGLE
          ! Implicit barrier

          ! Determine the team member number
          ! NB: this is numbered from 1 to follow the Fortran convention.
          ipar=1
!$ ipar = omp_get_thread_num()+1

          ! Set up meta-information for segments
          CALL segments_mod_seg_meta(meta_segments, ipar, num_parallel_sw,&
                  daylight_points, segment_size, a_sw_segments)

          ! Allocate storage
          ALLOCATE( segments ( meta_segments%num_segments ) )

          ! Find specific points in each segment
          IF (model_type == mt_single_column) THEN
            CALL segments_mod_segments(segments,meta_segments,                 &
                                       row_length, rows)

          ELSE
            CALL segments_mod_segments(segments,meta_segments,                 &
                                       row_length, rows,                       &
                                       list_points=list_daylight_points(:,j_sw))
          END IF

          DO i = 1, meta_segments%num_segments
            lit_points  = segments(i)%use_points
            start_point = segments(i)%start_index
            first_point = segments(i)%fp
            ii=segments(i)%first_x
            jj=segments(i)%first_y

            ! Additional work on list_points. Note that although work is done
            ! on this array, each thread and each segment should only affect
            ! one element.
            DO j = start_point, segments(i)%end_index
              list_daylight_points_start(j) = list_daylight_points(j,j_sw)&
                                             -segments(i)%start+1
            END DO

            ! Set the first point of the dust arrays to be used.
            IF (l_dust) THEN
              IF (l_twobin_dust) THEN
                first_point_dust_a = first_point
                first_point_dust_b = 1
              ELSE
                first_point_dust_a = first_point
                first_point_dust_b = first_point
              END IF
            ELSE
              first_point_dust_a = 1
              first_point_dust_b = 1
            END IF

            ! Set the first point of the biogenic array.
            IF (l_use_biogenic) THEN
              first_point_biogenic=first_point
            ELSE
              first_point_biogenic=1
            END IF

            ! Set the first point of the array of sulphate to be used.
            ! A separate assignment is necessary since this array will
            ! not be of the full size unless the sulphur cycle is on.
            IF (l_sulpc_so2 .OR. l_use_sulpc_indirect_sw) THEN
              first_point_sulpc=first_point
            ELSE
              first_point_sulpc=1
            END IF

            IF (l_soot) THEN
              first_point_soot=first_point
            ELSE
              first_point_soot=1
            END IF

            IF (l_biomass .OR. l_use_bmass_indirect) THEN
              first_point_biomass=first_point
            ELSE
              first_point_biomass=1
            END IF

            IF (l_ocff .OR. l_use_ocff_indirect) THEN
              first_point_ocff=first_point
            ELSE
              first_point_ocff=1
            END IF

            IF (n_arcl_species > 0) THEN
              first_point_arcl = first_point
            ELSE
              first_point_arcl = 1
            END IF

            IF (l_nitrate .OR. l_use_nitrate_indirect) THEN
              first_point_nitrate=first_point
            ELSE
              first_point_nitrate=1
            END IF

            IF (l_ukca_radaer .OR. l_glomap_clim_radaer) THEN
              first_point_ukca=first_point
            ELSE
              first_point_ukca=1
            END IF

            ! Set dynamically determined array sizes required for
            ! radiation calculations.
            CALL set_dimen(                                             &
              sw_control(j_sw), dimen, sw_spectrum(j_sw),               &
              row_length*rows, lit_points,                              &
              n_rad_layers, global_cloud_top, n_ukca_mode,              &
              j_l_use_easyaerosol,                                      &
              nd_field_flux_diag, nd_field_rad_diag)


            CALL sw_rad(                                                       &
            ! Mixing Ratios
              q_n(ii,jj,1), j_co2_mmr,j_ozone(ii,jj,1),                        &
              j_o2_mmr, co2_dim1, co2_dim2, co2_3d(ii,jj,1),                   &
              l_co2_3d, j_n2o_mmr, j_ch4_mmr, j_so2_mmr,                       &

            ! Chemical greenhouse gas fields
              ngrgas, grgas_field(ii,jj,1,1),                                  &

            ! Pressures and Temperatures
              p_star(ii,jj), p_layer_boundaries(ii,jj,0),                      &
              p_layer_centres(ii,jj,0),t_n(ii,jj,1),                           &
              t_layer_boundaries(ii,jj,0),                                     &
              p_extra_layer(ii,jj), t_extra_layer(ii,jj),                      &
              d_mass(ii,jj,1), density(ii,jj,1),                               &
              layer_heat_capacity(ii,jj,1),                                    &
              r_layer_centres(ii,jj,1), r_layer_boundaries(ii,jj,0),           &

            ! Options for COSP
              l_cosp,                                                          &

            ! Options for treating clouds
              l_inhom_cloud, inhom_cloud_sw, dp_corr_strat, dp_corr_conv,      &

            ! Stratiform Cloud Fields
              area_cloud_fraction(ii,jj,1),                                    &
              cf_n(ii,jj,1),                                                   &
              qcl_n(ii,jj,1), qcf_total(ii,jj,1),                              &
              n_drop_pot(ii,jj,1),                                             &

            ! Convective Cloud Fields
              cca(ii,jj,1), cclwp(ii,jj),                                      &
              ccw(ii,jj,1), lcbase(ii,jj),                                     &
              ccb(ii,jj), cct(ii,jj),                                          &

            ! Surface Fields
              land_albedo(ii,jj,1,j_sw),                                       &
              flandg(ii,jj), sea_ice_albedo(ii,jj,1,j_sw),                     &
              open_sea_albedo(ii,jj,1,1,j_sw),                                 &
              ice_fract(ii,jj), land_sea_mask(ii,jj),                          &
              land0p5(ii,jj), snow_depth(ii,jj),                               &

            ! Solar Fields
              cos_zen_rts(ii,jj), day_frac_rts(ii,jj), solar_constant,         &
              list_daylight_points_start(start_point), scs, sindec,            &
              cos_zen_sph_rts(ii,jj,0), day_frac_sph_rts(ii,jj,0),             &

            ! Aerosol Fields
              j_l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,        &
              zh(ii,jj), aero_bl_levels,                                       &
              j_l_dust, j_l_use_dust, dust_dim1, dust_dim2,                    &
              dust_1(first_point_dust_a,1), dust_2(first_point_dust_a,1),      &
              dust_3(first_point_dust_b,1), dust_4(first_point_dust_b,1),      &
              dust_5(first_point_dust_b,1), dust_6(first_point_dust_b,1),      &
              j_l_use_biogenic, biogenic_dim1, biogenic_dim2,                  &
              local_biogenic(first_point_biogenic, 1),                         &
              j_l_sulpc_so2, j_l_use_sulpc_direct, l_use_sulpc_indirect_sw,    &
              sulp_dim1, sulp_dim2,                                            &
              accum_sulphate(first_point_sulpc, 1),                            &
              aitken_sulphate(first_point_sulpc, 1),                           &
              diss_sulphate(first_point_sulpc, 1),                             &
              sea_salt_film, sea_salt_jet,                                     &
              l_use_seasalt_indirect, j_l_use_seasalt_direct,                  &
              salt_dim1, salt_dim2, salt_dim3, j_l_soot, j_l_use_soot_direct,  &
              soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),           &
              aged_soot(first_point_soot, 1), j_l_biomass,                     &
              j_l_use_bmass_direct, bmass_dim1, bmass_dim2,                    &
              fresh_bmass(first_point_biomass, 1),                             &
              aged_bmass(first_point_biomass, 1),                              &
              cloud_bmass(first_point_biomass, 1), l_use_bmass_indirect,       &
              j_l_ocff, j_l_use_ocff_direct, ocff_dim1, ocff_dim2,             &
              fresh_ocff(first_point_ocff, 1),aged_ocff(first_point_ocff, 1),  &
              cloud_ocff(first_point_ocff, 1), l_use_ocff_indirect,            &
              j_l_nitrate, j_l_use_nitrate_direct, nitrate_dim1, nitrate_dim2, &
              accum_nitrate(first_point_nitrate, 1),                           &
              diss_nitrate(first_point_nitrate, 1), l_use_nitrate_indirect,    &
              j_l_use_arcl, arcl_dim1, arcl_dim2, j_n_arcl_species,            &
              n_arcl_compnts,i_arcl_compnts,local_arcl(first_point_arcl,1,1),  &
              aerosol(ii,jj,1), j_l_murk_rad,                                  &
              j_l_use_ukca_radaer, j_l_use_glomap_clim_radaer,                 &
              ukca_radaer, ukca_dim1, ukca_dim2,                               &
              local_ukca_mmr(first_point_ukca, 1, 1),                          &
              local_ukca_cvl(first_point_ukca, 1, 1),                          &
              local_ukca_dry(first_point_ukca, 1, 1),                          &
              local_ukca_wet(first_point_ukca, 1, 1),                          &
              local_ukca_rho(first_point_ukca, 1, 1),                          &
              local_ukca_vol(first_point_ukca, 1, 1),                          &
              local_ukca_wtv(first_point_ukca, 1, 1),                          &
              local_ukca_nbr(first_point_ukca, 1, 1),                          &
              j_l_use_easyaerosol, easyaerosol_sw,                             &
              l_easyaerosol_cdnc, easyaerosol_cdnc,                            &

            ! time
              previous_time, seconds_since_midnight,                           &

            ! grid-dependent arrays
              true_latitude(ii,jj),                                            &
              true_longitude(ii,jj),                                           &
              obs_solid_angle(ii,jj),                                          &
              trans_solid_angle(ii,jj),                                        &
              dir_flux_to_trans(ii,jj),                                        &

            ! Level of tropopause
              trindx(ii,jj),                                                   &

            ! Spectrum
              sw_spectrum(j_sw),                                               &

            ! Algorithmic options
              sw_control(j_sw), timestep,                                      &

            ! All diagnostics and associated arrays
              sw_diag(j_sw),                                                   &
              diag_row_list(start_point,j_sw),diag_col_list(start_point,j_sw), &

            ! Dimensions
              dimen,                                                           &
              lit_points,segments(i)%seg_points,n_rad_layers,                  &
              global_cloud_top,ozone_levels,row_length,                        &
              rows,rows*row_length,max_n_swbands,nd_field_flux_diag,           &
              nd_field_rad_diag, n_cca_levels, n_ukca_mode, n_ukca_cpnt,       &

            ! Output data
              surf_down_sw_rts(ii,jj,1),                                       &
              flux_below_690nm_surf(ii,jj),                                    &
              netsw(ii,jj),                                                    &
              top_absorption(ii,jj),                                           &
              swsea(ii,jj),                                                    &
              sw_incs(ii,jj,0),                                                &

            ! COSP input arguments
              cosp_gbx, cosp_sgx)

          END DO              ! Loop over segments

          ! Deallocate the segmentation arrays
          DEALLOCATE(segments)

!$OMP END PARALLEL

          ! Deallocate band-by-band control options
          CALL deallocate_control(sw_control(j_sw))

        END IF                 ! If daylight_points>0

        !If autotuning is active, decide what to do with the
        !trial segment size and report the current status.
        IF (l_autotune_segments) THEN
          CALL autotune_stop_region(sw_autotune_state, daylight_points)
          CALL autotune_advance(sw_autotune_state)
          CALL autotune_report(sw_autotune_state, quiet=(j_sw > 1))
        END IF

        IF (model_type /= mt_single_column) THEN

          ! Radiative fluxes may not have been calculated at all
          ! points: we now fill in as required. In the case of
          ! spatial degradation calls to SWAPBOUNDS are made and
          ! these must be made on all PEs, which requires that
          ! this code shall executed even when the PE contains
          ! no daylit points.


          l_complete_north=.FALSE.
          l_complete_south=.FALSE.

          ! When spatial degradation is performed fields must be
          ! filled in at alternate points.
          l_complete_deg = ( l_rad_deg ) .AND.                          &
                   ( tot_daylight_points > 0 )

          ! Set addressing limits for spatial degradation.
          IF ( l_complete_deg ) THEN
            first_row=1
            last_row=rows
          END IF

          ! Call appropriate subroutines to fill in missing data
          ! as required.
          IF ( l_complete_north .OR. l_complete_south .OR.              &
               l_complete_deg ) THEN

            CALL fill_missing_data_sw(                                  &
              off_x, off_y, row_length, rows, model_levels, ntiles,     &
              salt_dim1, salt_dim2, salt_dim3,                          &
              cloud_levels,first_row,last_row,                          &
              first_data_interp, es_space_interp,                       &
              l_complete_north, l_complete_south,                       &
              l_complete_deg, sw_control(j_sw)%l_blue_flux_surf,        &
              sw_control(j_sw)%n_channel, j_sw, l_extra_top,            &
              sw_incs, netsw, swsea, flux_below_690nm_surf,             &
              top_absorption, surf_down_sw_rts,                         &
              sea_salt_film, sea_salt_jet)
          END IF

        END IF ! model_type

        IF (l_rad_perturb .AND. (j_sw==1)) THEN
          sw_incs_local(:,:,:,1) = sw_incs                                     &
                                 - sw_incs_local(:,:,:,2)
          netsw_local(:,:,1) = netsw                                           &
                             - netsw_local(:,:,2)
          swsea_local(:,:,1) = swsea                                           &
                             - swsea_local(:,:,2)
          surf_down_sw_local(:,:,:,1) = surf_down_sw_rts                       &
                                      - surf_down_sw_local(:,:,:,2)
          flux_b690nm_local(:,:,1) = flux_below_690nm_surf                     &
                                   - flux_b690nm_local(:,:,2)
          IF (l_extra_top) THEN
            top_abs_sw(:,:,1) = top_absorption                                 &
                              - top_abs_sw(:,:,2)
          END IF
        ELSE IF (l_timestep) THEN
          sw_incs_local(:,:,:,j_sw)=sw_incs
          netsw_local(:,:,j_sw)=netsw
          swsea_local(:,:,j_sw)=swsea
          surf_down_sw_local(:,:,:,j_sw)=surf_down_sw_rts
          flux_b690nm_local(:,:,j_sw)=flux_below_690nm_surf
          IF (l_extra_top) THEN
            top_abs_sw(:,:,j_sw)=top_absorption
          END IF
        END IF

        IF (l_forcing .AND. j_sw > 1) THEN
          ! Reset grgas_field array to original values
          IF (c2c_n2o .AND. ngrgas >= p_n2o) THEN
            grgas_field(:,:,:,p_n2o) = n2o_mmr_tmp(:,:,:)
            DEALLOCATE(n2o_mmr_tmp)
          END IF
          IF (c2c_ch4 .AND. ngrgas >= p_ch4) THEN
            grgas_field(:,:,:,p_ch4) = ch4_mmr_tmp(:,:,:)
            DEALLOCATE(ch4_mmr_tmp)
          END IF
        END IF

      END IF ! l_call_swrad
    END DO ! radiation calls

    IF (l_rad_perturb .AND. l_rad_step_prog) THEN
      ! Do nothing: for this case the full fluxes are already set.
    ELSE IF (l_timestep .AND. l_rad_step_diag) THEN
      sw_incs(:,:,:) = sw_incs_local(:,:,:,1)                                  &
                     + sw_incs_local(:,:,:,2)
      netsw(:,:)     = netsw_local(:,:,1)                                      &
                     + netsw_local(:,:,2)
      swsea(:,:)     = swsea_local(:,:,1)                                      &
                     + swsea_local(:,:,2)
      surf_down_sw_rts(:,:,:) = surf_down_sw_local(:,:,:,1)                    &
                              + surf_down_sw_local(:,:,:,2)
      flux_below_690nm_surf(:,:) = flux_b690nm_local(:,:,1)                    &
                                 + flux_b690nm_local(:,:,2)
      IF (l_extra_top) THEN
        top_absorption(:,:) = top_abs_sw(:,:,1)                                &
                            + top_abs_sw(:,:,2)
      END IF
    END IF

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, mean_cos_zenith_angle,                        &
!$OMP&         cos_zen_rts, day_frac_rts )                                     &
!$OMP& PRIVATE( i, j )
    DO j = 1, rows
      DO i = 1, row_length
        mean_cos_zenith_angle(i,j) = cos_zen_rts(i,j)*day_frac_rts(i,j)
      END DO
    END DO
!$OMP END PARALLEL DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j )
    IF (l_rad_perturb .AND. (.NOT. l_rad_step_prog)) THEN
      IF (sw_control(1)%l_spherical_solar) THEN
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            IF ( sw_incs(i,j,0) < 0.0 ) THEN
              netsw(i,j) = netsw(i,j) - sw_incs(i,j,0)
              sw_incs(i,j,0) = 0.0
              swsea(i,j) = 0.0
            END IF
          END DO
        END DO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            IF ( sw_incs(i,j,0) < 0.0 ) THEN
              netsw(i,j) = netsw(i,j)-sw_incs(i,j,0)*mean_cos_zenith_angle(i,j)
              sw_incs(i,j,0) = 0.0
            END IF
          END DO
        END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            IF ( swsea(i,j) < 0.0 ) THEN
              netsw(i,j) = netsw(i,j) - (1.0-flandg(i,j))*swsea(i,j)
              swsea(i,j) = 0.0
            END IF
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          surf_down_sw_rts(i,j,1) = MAX(surf_down_sw_rts(i,j,1), 0.0)
          surf_down_sw_rts(i,j,2) = MAX(surf_down_sw_rts(i,j,2), 0.0)
          surf_down_sw_rts(i,j,3) = MAX(surf_down_sw_rts(i,j,3), 0.0)
          surf_down_sw_rts(i,j,4) = MAX(surf_down_sw_rts(i,j,4), 0.0)
          flux_below_690nm_surf(i,j) = MAX(flux_below_690nm_surf(i,j), 0.0)
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    ! Calculate net surface SW for diagnostic
    IF (sw_control(1)%l_spherical_solar) THEN
      IF (l_ctile) THEN
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            surfsw(i,j) = sw_incs(i,j,0)
            IF (flandg(i,j) == 1.0) swsea(i,j)=rmdi
          END DO
        END DO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            surfsw(i,j) = sw_incs(i,j,0)
            IF (land_sea_mask(i,j)) swsea(i,j)=rmdi
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF
    ELSE
      IF (l_ctile) THEN
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            surfsw(i,j) = sw_incs(i,j,0) * mean_cos_zenith_angle(i,j)          &
              + (1.0-flandg(i,j))*swsea(i,j)
            IF (flandg(i,j) == 1.0) swsea(i,j)=rmdi
          END DO
        END DO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            surfsw(i,j) = sw_incs(i,j,0) * mean_cos_zenith_angle(i,j)          &
              + swsea(i,j)
            IF (land_sea_mask(i,j)) swsea(i,j)=rmdi
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF
!$OMP BARRIER
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          IF (mean_cos_zenith_angle(i,j) > eps) THEN
            sw_incs(i,j,0)=surfsw(i,j)/mean_cos_zenith_angle(i,j)
          ELSE
            sw_incs(i,j,0)=0.0
          END IF
        END DO
      END DO
!$OMP END DO
    END IF
!$OMP END PARALLEL


    IF (.NOT. l_planet_grey_surface) THEN
      ! Calculate net surface SW on land tiles
!$OMP  PARALLEL DO IF(ntiles > 1) DEFAULT(NONE) SCHEDULE(STATIC)               &
!$OMP& SHARED( ntiles, land_field, sw_tile_rts, tile_pts, tile_index,          &
!$OMP&         land_index_i, land_index_j, surf_down_sw_rts, alb_tile )        &
!$OMP& PRIVATE( i, j, k, l, n, point )
      DO n = 1, ntiles
        DO l = 1, land_field
          sw_tile_rts(l,n) = 0.0
        END DO
        DO point = 1, tile_pts(n)
          l = tile_index(point,n)
          i = land_index_i(l)
          j = land_index_j(l)
          DO k = 1, 4
            sw_tile_rts(l,n) = sw_tile_rts(l,n)                         &
               + (1.0 - alb_tile(l,n,k))*surf_down_sw_rts(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Calculate net surface SW on open sea
      DO l = 1, ssi_pts
        sw_rts_sea(l) = 0.0
      END DO
      DO point = 1, sea_pts
        l = sea_index(point)
        i = ssi_index_i(l)
        j = ssi_index_j(l)
        sw_rts_sea(l) = swsea(i,j)
      END DO

      ! Calculate net, upward and downward surface SW on sea ice categories
!$OMP  PARALLEL DO IF(nice_use > 1) DEFAULT(NONE) SCHEDULE(STATIC)             &
!$OMP& SHARED( nice_use, ssi_pts, sw_rts_sicat, sice_pts_ncat, sice_index_ncat,&
!$OMP&         ssi_index_i, ssi_index_j, alb_sicat, surf_down_sw_rts,          &
!$OMP&         swup_rts_sicat,swdn_rts_sicat  )                                &
!$OMP& PRIVATE( i, j, k, l, n, point )
      DO n = 1, nice_use
        DO l = 1, ssi_pts
          sw_rts_sicat(l,n) = 0.0
          swup_rts_sicat(l,n) = 0.0
          swdn_rts_sicat(l,n) = 0.0
        END DO
        DO point = 1, sice_pts_ncat(n)
          l = sice_index_ncat(point,n)
          i = ssi_index_i(l)
          j = ssi_index_j(l)
          DO k = 1, 4
            sw_rts_sicat(l,n) = sw_rts_sicat(l,n)                       &
              + (1.0 - alb_sicat(l,n,k))*surf_down_sw_rts(i,j,k)
            swup_rts_sicat(l,n) = swup_rts_sicat(l,n)                   &
                             + alb_sicat(l,n,k)*surf_down_sw_rts(i,j,k)
            swdn_rts_sicat(l,n) = swdn_rts_sicat(l,n)                   &
                                              + surf_down_sw_rts(i,j,k)
          END DO
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Set indexing for albedo arrays appropriate to timestep
      IF (l_rad_step_prog) THEN
        j_sw = 1
      ELSE
        j_sw = 2
      END IF

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, land_alb, sice_alb, surf_down_sw_sum,         &
!$OMP&         surf_down_sw_rts )                                              &
!$OMP& PRIVATE( i, j )
      DO j = 1, rows
        DO i = 1, row_length
          land_alb(i,j)=0.0
          sice_alb(i,j) = 0.0
          surf_down_sw_sum(i,j)=SUM(surf_down_sw_rts(i,j,:))
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Calculate mean land albedo
      DO l = 1, land_field
        i = land_index_i(l)
        j = land_index_j(l)
        IF (surf_down_sw_sum(i,j) > eps) THEN
          DO k = 1, 4
            land_alb(i,j) = land_alb(i,j)                               &
              + land_albedo(i,j,k,j_sw)*surf_down_sw_rts(i,j,k)
          END DO
          land_alb(i,j) = land_alb(i,j)/surf_down_sw_sum(i,j)
        END IF
      END DO

      ! Calculate mean albedo of lake tile
      IF (l_flake_model .AND. (.NOT. l_aggregate)) THEN
        ! First zero the array
        lake_albedo_gb(:) = 0.0
        n = lake
        DO point = 1, tile_pts(n)
          l = tile_index(point,n)
          i = land_index_i(l)
          j = land_index_j(l)
          IF (surf_down_sw_sum(i,j) > eps) THEN
            lake_albedo_gb(l)=1.0-(sw_tile_rts(l,n)/surf_down_sw_sum(i,j))
          END IF
        END DO
      END IF

      ! Calculate mean sea ice albedo
      DO point = 1, sice_pts
        l = sice_index(point)
        i = ssi_index_i(l)
        j = ssi_index_j(l)
        IF (surf_down_sw_sum(i,j) > eps) THEN
          DO k = 1, 4
            sice_alb(i,j) = sice_alb(i,j)                               &
              + sea_ice_albedo(i,j,k,j_sw)*surf_down_sw_rts(i,j,k)
          END DO
          sice_alb(i,j) = sice_alb(i,j)/surf_down_sw_sum(i,j)
        END IF
      END DO
    END IF ! l_planet_grey_surface

  END IF ! on a radiation timestep


  ! Calculate day fraction and mean cos(solar zenith angle while
  ! the sun is up) for each grid point for this physics timestep:
  ! (if in fact full SW calculations are being done every timestep, this
  ! is of course unnecessary, as are various calculations later on)
  CALL solang(                                                          &
  ! input constants
        sindec, seconds_since_midnight,                                 &
        timestep, eq_time,                                              &
  ! row and column dependent constants
        true_latitude,                                                  &
        true_longitude,                                                 &
  ! size variables
        row_length*rows,                                                &
  ! output fields
        day_fraction, cos_zenith_angle,                                 &
        sol_azimuth, cosz_beg, cosz_end)

  DO j = 1, rows
    DO i = 1, row_length
      surf_down_sw_sum(i,j) = SUM(surf_down_sw_rts(i,j,:))
      surfdir_rts(i,j) = surf_down_sw_rts(i,j,1) + surf_down_sw_rts(i,j,3)
    END DO
  END DO

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j )
  IF (l_rad_szacor .AND. .NOT. sw_control(1)%l_spherical_solar) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        IF (cos_zen_rts(i,j)*day_frac_rts(i,j) < eps) THEN
          cos_zenith_angle(i,j) = 0.0
          day_fraction(i,j) = 0.0
        END IF
      END DO
    END DO
!$OMP END DO NOWAIT
    IF (l_orog) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          IF ( (surfdir_rts(i,j) > eps) .AND.                                  &
               (surfdir_rts(i,j) < orog_corr(i,j)*scs*solcon_rts) .AND.        &
               (cos_zen_rts(i,j) > eps) .AND.                                  &
               (cos_zenith_angle(i,j) > eps) .AND.                             &
               (orog_corr(i,j) > SQRT(eps)) ) THEN
            sza_cor(i,j) = (orog_corr(i,j) - 0.5) * (scs*solcon_rts            &
                   * EXP(LOG(surfdir_rts(i,j)/(orog_corr(i,j)*scs*solcon_rts)) &
                   *(cos_zen_rts(i,j)/cos_zenith_angle(i,j)))                  &
                   - surfdir_rts(i,j)/orog_corr(i,j))/surf_down_sw_sum(i,j)
          ELSE
            sza_cor(i,j) = 0.0
          END IF
        END DO
      END DO
!$OMP END DO NOWAIT
    ELSE
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          IF ( (surfdir_rts(i,j) > eps) .AND.                                  &
               (surfdir_rts(i,j) < scs*solcon_rts) .AND.                       &
               (cos_zen_rts(i,j) > eps) .AND.                                  &
               (cos_zenith_angle(i,j) > eps) ) THEN
            sza_cor(i,j) = 0.5*(scs*solcon_rts                                 &
                   * EXP(LOG(surfdir_rts(i,j)/(scs*solcon_rts))                &
                   *(cos_zen_rts(i,j)/cos_zenith_angle(i,j)))                  &
                   - surfdir_rts(i,j))/surf_down_sw_sum(i,j)
          ELSE
            sza_cor(i,j) = 0.0
          END IF
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF
  ELSE
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        sza_cor(i,j) = 0.0
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      ! Combine the two terms to give the mean cos zenith angle over the
      ! whole of the physics timestep.
      cos_zenith_angle(i,j) = cos_zenith_angle(i,j) * day_fraction(i,j)
      ! Calculate incoming SW at top of atmosphere
      itoasw(i,j) = cos_zenith_angle(i,j) * scs * solcon_rts
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

  IF (sw_control(1)%l_spherical_solar) THEN

    ! Currently no timestep correction is done for the increments
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          delta_t(i,j,k) = sw_incs(i,j,k)
        END DO
      END DO
    END DO

    ! Calculate surface downward SW flux corrected to model timestep
    DO j = 1, rows
      DO i = 1, row_length
        sza_cor(i,j) = cos_zenith_angle(i,j) &
          / MAX(cos_zen_rts(i,j)*day_frac_rts(i,j), eps)
        surf_down_sw(i,j,1)=surf_down_sw_rts(i,j,1) * sza_cor(i,j)
        surf_down_sw(i,j,2)=surf_down_sw_rts(i,j,2)
        surf_down_sw(i,j,3)=surf_down_sw_rts(i,j,3) * sza_cor(i,j)
        surf_down_sw(i,j,4)=surf_down_sw_rts(i,j,4)
      END DO
    END DO

    ! Set up photosynthetically active surface radiation
    DO j = 1, rows
      DO i = 1, row_length
        photosynth_act_rad(i,j) = surf_down_sw(i,j,1) + surf_down_sw(i,j,2)
      END DO
    END DO

    ! Calculate corrected net surface SW for diagnostic:
    IF (l_surfsw_cor .OR. l_planet_intrinsic_flux) THEN
      DO j = 1, rows
        DO i = 1, row_length
          surfsw_cor(i,j) = sw_incs(i,j,0) * SUM(surf_down_sw(i,j,:)) &
            / MAX(surf_down_sw_sum(i,j), eps)
        END DO
      END DO
    END IF

    ! Calculate corrected TOA outgoing SW for diagnostic:
    IF (l_toasw_cor) THEN
      DO j = 1, rows
        DO i = 1, row_length
          toasw_cor(i,j) = -sw_incs(i,j,model_levels+1)
        END DO
      END DO
    END IF

    ! Calculate corrected direct surface SW for diagnostic:
    IF (l_surfdir_cor) THEN
      DO j = 1, rows
        DO i = 1, row_length
          surfdir_cor(i,j) = surf_down_sw(i,j,1) + surf_down_sw(i,j,3)
        END DO
      END DO
    END IF

    ! Calculate corrected diffuse surface SW for diagnostic:
    IF (l_surfdif_cor) THEN
      DO j = 1, rows
        DO i = 1, row_length
          surfdif_cor(i,j) = surf_down_sw(i,j,2) + surf_down_sw(i,j,4)
        END DO
      END DO
    END IF

    IF (.NOT. l_planet_grey_surface) THEN
      ! Set up net surface SW on land tiles
      DO n = 1, ntiles
        DO point = 1, tile_pts(n)
          l = tile_index(point,n)
          i = land_index_i(l)
          j = land_index_j(l)
          sw_tile(l,n) &
            = (1.0 - alb_tile(l,n,1))*surf_down_sw(i,j,1) &
            + (1.0 - alb_tile(l,n,2))*surf_down_sw(i,j,2) &
            + (1.0 - alb_tile(l,n,3))*surf_down_sw(i,j,3) &
            + (1.0 - alb_tile(l,n,4))*surf_down_sw(i,j,4)
          sw_net_land(i,j) = sw_net_land(i,j) &
            + sw_tile(l,n) * tile_frac(l,n)
        END DO
      END DO

      ! Set up net surface SW on sea ice categories
      DO n = 1, nice_use
        DO point = 1, sice_pts_ncat(n)
          l = sice_index_ncat(point,n)
          i = ssi_index_i(l)
          j = ssi_index_j(l)
          swdn_sicat(l,n) = SUM(surf_down_sw(i,j,:))
          swup_sicat(l,n) &
            = alb_sicat(l,n,1)*surf_down_sw(i,j,1) &
            + alb_sicat(l,n,2)*surf_down_sw(i,j,2) &
            + alb_sicat(l,n,3)*surf_down_sw(i,j,3) &
            + alb_sicat(l,n,4)*surf_down_sw(i,j,4)
          sw_sicat(l,n) = swdn_sicat(l,n) - swup_sicat(l,n)
          sw_net_sice(i,j) = sw_net_sice(i,j)                             &
            + sw_sicat(l,n) * ice_fract_cat(i,j,n)
        END DO
      END DO
    END IF

  ELSE

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j, k )
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          delta_t(i,j,k) = sw_incs(i,j,k) * cos_zenith_angle(i,j)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

    ! Calculate surface downward SW flux corrected to model timestep
    DO k = 1, 4
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          surf_down_sw(i,j,k)=surf_down_sw_rts(i,j,k) * cos_zenith_angle(i,j)
        END DO
      END DO
!$OMP END DO NOWAIT
    END DO

    ! Set up photosynthetically active surface radiation
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        photosynth_act_rad(i,j) = ( surf_down_sw_rts(i,j,1) + &
          surf_down_sw_rts(i,j,2) ) * cos_zenith_angle(i,j)
      END DO
    END DO
!$OMP END DO NOWAIT

    ! Calculate corrected net surface SW for diagnostic:
    IF (l_surfsw_cor .OR. l_planet_intrinsic_flux) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          surfsw_cor(i,j) = &
            (1.0+sza_cor(i,j))*cos_zenith_angle(i,j)*sw_incs(i,j,0)
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    ! Calculate corrected TOA outgoing SW for diagnostic:
    IF (l_toasw_cor) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          toasw_cor(i,j) = MAX( cos_zenith_angle(i,j)                          &
            * ( scs * solcon_rts - sw_incs(i,j,model_levels+1)                 &
              - sza_cor(i,j) * sw_incs(i,j,0) ), 0.0)
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    ! Calculate corrected direct surface SW for diagnostic:
    IF (l_surfdir_cor .OR. l_surfdif_cor) THEN
      IF (l_orog .AND. l_rad_szacor) THEN
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            surfdir_cor(i,j) = MAX(cos_zenith_angle(i,j)                       &
              * (surfdir_rts(i,j) + orog_corr(i,j) * sza_cor(i,j)              &
              * surf_down_sw_sum(i,j)/(orog_corr(i,j)-0.5) ), 0.0)
          END DO
        END DO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            surfdir_cor(i,j) = MAX(cos_zenith_angle(i,j)*(surfdir_rts(i,j) +   &
              2.0 * sza_cor(i,j) * surf_down_sw_sum(i,j)), 0.0)
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF
    END IF

    ! Calculate corrected diffuse surface SW for diagnostic:
    IF (l_surfdif_cor) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          surfdif_cor(i,j) = MAX( (1.0+sza_cor(i,j)) * cos_zenith_angle(i,j)   &
            * surf_down_sw_sum(i,j) - surfdir_cor(i,j), 0.0)
        END DO
      END DO
!$OMP END DO
    END IF
!$OMP END PARALLEL

    IF (.NOT. l_planet_grey_surface) THEN
      ! Set up net surface SW on land tiles
      DO n = 1, ntiles
        DO point = 1, tile_pts(n)
          l = tile_index(point,n)
          i = land_index_i(l)
          j = land_index_j(l)
          sw_tile(l,n) = sw_tile_rts(l,n)                                 &
             * cos_zenith_angle(i,j) * (1.0+sza_cor(i,j))
          sw_net_land(i,j) = sw_net_land(i,j)                             &
             + sw_tile(l,n) * tile_frac(l,n)
        END DO
      END DO

      ! Set up net surface SW on sea ice categories
      DO n = 1, nice_use
        DO point = 1, sice_pts_ncat(n)
          l = sice_index_ncat(point,n)
          i = ssi_index_i(l)
          j = ssi_index_j(l)
          sw_sicat(l,n) = sw_rts_sicat(l,n)                               &
             * cos_zenith_angle(i,j) * (1.0+sza_cor(i,j))
          swdn_sicat(l,n) = swdn_rts_sicat(l,n)                           &
             * cos_zenith_angle(i,j) * (1.0+sza_cor(i,j))
          swup_sicat(l,n) = swup_rts_sicat(l,n)                           &
             * cos_zenith_angle(i,j) * (1.0+sza_cor(i,j))
          sw_net_sice(i,j) = sw_net_sice(i,j)                             &
             + sw_sicat(l,n) * ice_fract_cat(i,j,n)
        END DO
      END DO
    END IF ! l_planet_grey_surface

  END IF


  ! Set up radiative heating rates for 6A boundary layer code
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE( i, j, k )
  DO k = 1, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        rad_hr(i,j,2,k) = delta_t(i,j,k) * recip_timestep
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO


  ! Is the PC2 cloud scheme being used?
  IF (i_cld_vn == i_cld_pc2) THEN

    ! Reset _latest values to _n values
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE ( i, j, k )
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_latest(i,j,k)   = t_n(i,j,k)
          q_latest(i,j,k)   = q_n(i,j,k)
          qcl_latest(i,j,k) = qcl_n(i,j,k)
          cf_latest(i,j,k)  = cf_n(i,j,k)
          cfl_latest(i,j,k) = cfl_n(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------
    ! Homogeneous forcing. Note the temperature increment from the shortwave
    ! is added in this routine
    ! ----------------------------------------------------------------------

    CALL pc2_homog_plus_turb(p_layer_centres(1,1,1),                    &
      model_levels,                                                     &
      timestep, t_latest, cf_latest, cfl_latest,                        &
      cff_latest, q_latest, qcl_latest, delta_t(1,1,1),                 &
      zeros, zeros, zeros, 0.0, 0.0,                                    &
      l_mr_physics)

    ! Add increments from the homogeneous forcing to the increment variables
!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE ( i, j, k )
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_inc(i,j,k) = t_inc(i,j,k) + t_latest(i,j,k)-t_n(i,j,k)

          q_inc(i,j,k) = q_inc(i,j,k) + q_latest(i,j,k)-q_n(i,j,k)
          qcl_inc(i,j,k) = qcl_inc(i,j,k)                               &
                           + qcl_latest(i,j,k)-qcl_n(i,j,k)
          cf_inc(i,j,k) = cf_inc(i,j,k)                                 &
                           + cf_latest(i,j,k)-cf_n(i,j,k)
          cfl_inc(i,j,k) = cfl_inc(i,j,k)                               &
                           + cfl_latest(i,j,k)-cfl_n(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
    IF (l_retain_rad_tendencies) THEN
      DO k = 1, model_levels
        DO j = 1, rows
          DO i = 1, row_length
            dq_sw(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
          END DO
        END DO
      END DO
    END IF !end if over l_retain_rad_tendencies

  ELSE  ! i_cld_pc2

    ! add SW radiative heating to temperatures
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_inc(i,j,k) = t_inc(i,j,k) + delta_t(i,j,k)
          t_latest(i,j,k) = t_n(i,j,k) + delta_t(i,j,k)
        END DO
      END DO
    END DO

  END IF  ! i_cld_pc2

  ! Copy SW radiation tendency to physics_tendencies_mod
  IF (l_retain_rad_tendencies) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          dt_sw(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
        END DO
      END DO
    END DO
  END IF

  ! Get T_incr for output as STASH diagnostic

  IF ( l_t_incr_sw ) THEN
    ! Increments will be calculated if to be diagnosed directly
    ! or to be used to determine heating rates.
    ALLOCATE ( t_incr_diagnostic(row_length,rows,model_levels) )
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_incr_diagnostic(i,j,k) =  delta_t(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k


    IF ( l_scmdiags(scmdiag_pc2) .AND.                                  &
         model_type == mt_single_column ) THEN

      ! Calculate equivalent to stash code 1,181 - the total
      ! temperature increment including the pc2 scheme
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            TmpScm3d_2(i,j,k) = t_latest(i,j,k) - t_n(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

      ! stash item = 1,181  ! temperature increment minus pc2
      CALL scmoutput(TmpScm3d_2,'dt_swpc2',                             &
           'SW heating rate incl PC2','K/timestep',                     &
           t_acc,d_all,default_streams,'',routinename)
      CALL scmoutput(TmpScm3d_2,'sw2pc2',                               &
           'SW heating rate incl PC2','K/day',                          &
           t_mult,d_all,default_streams,'ntspday',routinename)
    END IF ! scmdiag_pc2 / model_type


  ELSE
    ALLOCATE ( t_incr_diagnostic(1,1,1) )
  END IF                    ! on STASHflag / TmpLogic


  ! Set up photosynthetically active surface radiation
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j )
  IF (l_direct_par) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        flxdirparsurf(i,j) = surf_down_sw(i,j,1)
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF
  IF (l_rad_szacor .AND. .NOT. sw_control(1)%l_spherical_solar) THEN
    IF (l_direct_par) THEN
      IF (l_orog) THEN
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            flxdirparsurf(i,j) = flxdirparsurf(i,j) + orog_corr(i,j)           &
                                   * sza_cor(i,j) * photosynth_act_rad(i,j)    &
                                   /(orog_corr(i,j)-0.5)
          END DO
        END DO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO SCHEDULE(STATIC)
        DO j = 1, rows
          DO i = 1, row_length
            flxdirparsurf(i,j) = flxdirparsurf(i,j) +                          &
                                 2.0 * sza_cor(i,j) * photosynth_act_rad(i,j)
          END DO
        END DO
!$OMP END DO NOWAIT
      END IF
    END IF
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        photosynth_act_rad(i,j) = photosynth_act_rad(i,j) * (1.0+sza_cor(i,j))
      END DO
    END DO
!$OMP END DO
  END IF
!$OMP END PARALLEL

  IF (.NOT. l_planet_grey_surface) THEN
    ! Set up net surface SW on open sea
    DO point = 1, sea_pts
      l = sea_index(point)
      i = ssi_index_i(l)
      j = ssi_index_j(l)
      sw_sea(l) = sw_rts_sea(l)
    END DO

    ! SW diagnostics weighted by sea ice fraction
    IF (l_swup_sice_wt) THEN
      ALLOCATE ( sw_up_sice_weighted(row_length,rows) )
      sw_up_sice_weighted(:,:) = 0.0
    END IF
    IF (l_swup_sice_wt_cat) THEN
      ALLOCATE ( sw_up_sice_weighted_cat(row_length,rows,nice_use) )
      sw_up_sice_weighted_cat(:,:,:) = 0.0
    END IF
    IF (l_swdn_sice_wt) THEN
      ALLOCATE ( sw_down_sice_weighted(row_length,rows) )
      sw_down_sice_weighted(:,:) = 0.0
    END IF
    IF (l_swdn_sice_wt_cat) THEN
      ALLOCATE ( sw_down_sice_weighted_cat(row_length,rows,nice_use) )
      sw_down_sice_weighted_cat(:,:,:) = 0.0
    END IF
    IF (l_alb_sice_wt) THEN
      ALLOCATE ( albedo_sice_weighted(row_length,rows) )
      albedo_sice_weighted(:,:) = 0.0
    END IF
    IF (l_alb_sice_wt_cat) THEN
      ALLOCATE ( albedo_sice_weighted_cat(row_length,rows,nice_use) )
      albedo_sice_weighted_cat(:,:,:) = 0.0
    END IF
    DO n = 1, nice_use
      DO point = 1, sice_pts_ncat(n)
        l = sice_index_ncat(point,n)
        i = ssi_index_i(l)
        j = ssi_index_j(l)
        IF (l_swdn_sice_wt_cat) sw_down_sice_weighted_cat(i,j,n) =      &
                                 swdn_sicat(l,n) * ice_fract_cat(i,j,n)
        IF (l_swdn_sice_wt) sw_down_sice_weighted(i,j) =                &
                                          sw_down_sice_weighted(i,j) +  &
                                 swdn_sicat(l,n) * ice_fract_cat(i,j,n)
        IF (l_swup_sice_wt_cat) sw_up_sice_weighted_cat(i,j,n) =        &
                                 swup_sicat(l,n) * ice_fract_cat(i,j,n)
        IF (l_swup_sice_wt) sw_up_sice_weighted(i,j) =                  &
                                           sw_up_sice_weighted(i,j) +   &
                                 swup_sicat(l,n) * ice_fract_cat(i,j,n)
        IF (l_alb_sice_wt_cat) THEN
          IF (swdn_sicat(l,n) >= 0.8) THEN
            albedo_sice_weighted_cat(i,j,n) = swup_sicat(l,n)           &
                           * ice_fract_cat(i,j,n) / swdn_sicat(l,n)
          ELSE
            albedo_sice_weighted_cat(i,j,n) = 0.0
          END IF
        END IF
        IF (l_alb_sice_wt) THEN
          IF (swdn_sicat(l,n) >= 0.8) THEN
            albedo_sice_weighted(i,j) = albedo_sice_weighted(i,j) +     &
                swup_sicat(l,n) * ice_fract_cat(i,j,n) / swdn_sicat(l,n)
          END IF
        END IF
      END DO
    END DO
  END IF ! l_planet_grey_surface

  IF (ltimer) CALL timer ('AP1R SW Rad  ',6)

  ! ----------------------------------------------------------------------
  ! Section RAD.1.1 Short Wave Radiation Energy correction code
  ! -----------------------------------------------------------------------

  IF ((l_rad_step_prog .OR. l_rad_step_diag) .AND. l_emcorr) THEN

    ! Sum short wave fluxes into the atmosphere and
    ! add into the net diabatic fluxes into the
    ! atmosphere for use in the energy correction
    ! procedure

    IF (l_extra_top) THEN
      ! The energy absorbed above the top of the model in
      ! the radiation scheme does not contribute to the
      ! energy absorbed, but the diagnostics are calculated
      ! at the top of the atmosphere, so the net atmospheric
      ! flux must be adjusted.
      DO j = 1, rows
        DO i = 1, row_length
          net_atm_flux(i,j) = netsw(i,j) - surfsw(i,j)                  &
                            - top_absorption(i,j)
        END DO
      END DO
    ELSE
      DO j = 1, rows
        DO i = 1, row_length
          net_atm_flux(i,j) = netsw(i,j) - surfsw(i,j)
        END DO
      END DO
    END IF

    IF (l_orog .AND. .NOT. sw_control(1)%l_spherical_solar) THEN
      net_atm_flux = net_atm_flux + f_orog
    END IF

    CALL flux_diag(net_atm_flux, cos_theta_latitude,                    &
      row_length, rows ,off_x,off_y, 1.0,                               &
      sum_eng_fluxes,radiation_tstep_diag)

  END IF

  ! ----------------------------------------------------------------------
  ! Section RAD.1.2 Short Wave Radiation diagnostics
  ! -----------------------------------------------------------------------

  ! Minimal allocation for unused diagnostics
  IF (.NOT. ALLOCATED( sw_up_sice_weighted) ) &
            ALLOCATE ( sw_up_sice_weighted(1,1) )
  IF (.NOT. ALLOCATED( sw_up_sice_weighted_cat) ) &
            ALLOCATE ( sw_up_sice_weighted_cat(1,1,1) )
  IF (.NOT. ALLOCATED( sw_down_sice_weighted) ) &
            ALLOCATE ( sw_down_sice_weighted(1,1) )
  IF (.NOT. ALLOCATED( sw_down_sice_weighted_cat) ) &
            ALLOCATE ( sw_down_sice_weighted_cat(1,1,1) )
  IF (.NOT. ALLOCATED( albedo_sice_weighted) ) &
            ALLOCATE ( albedo_sice_weighted(1,1) )
  IF (.NOT. ALLOCATED( albedo_sice_weighted_cat) ) &
            ALLOCATE ( albedo_sice_weighted_cat(1,1,1) )

  DO j_sw = 1, n_swcall

    ! Check that sw diagnostics requested this timestep


    SELECT CASE (model_type)
    CASE (mt_single_column)
      TmpLogic = ( error_code == 0 )
    CASE DEFAULT
      TmpLogic = ( ( error_code == 0 ) .AND. sf_calc(0,1) )
    END SELECT


    IF (TmpLogic) THEN

      IF (l_timestep) THEN

        i_off=0

        IF (j_sw == n_swcall) THEN

          IF (l_rad_perturb .AND. l_rad_step_prog) THEN

            IF (sw_diag(2)%l_flux_up)                                   &
                sw_diag(1)%flux_up =                                    &
                sw_diag(1)%flux_up -                                    &
                sw_diag(2)%flux_up

            IF (sw_diag(2)%l_flux_down)                                 &
                sw_diag(1)%flux_down =                                  &
                sw_diag(1)%flux_down -                                  &
                sw_diag(2)%flux_down

            IF (sw_diag(2)%l_solar_out_toa)                             &
                sw_diag(1)%solar_out_toa =                              &
                sw_diag(1)%solar_out_toa -                              &
                sw_diag(2)%solar_out_toa

            IF (sw_diag(2)%l_surface_down_flux)                         &
                sw_diag(1)%surface_down_flux =                          &
                sw_diag(1)%surface_down_flux -                          &
                sw_diag(2)%surface_down_flux

            IF (sw_diag(2)%l_net_flux_trop)                             &
                sw_diag(1)%net_flux_trop =                              &
                sw_diag(1)%net_flux_trop -                              &
                sw_diag(2)%net_flux_trop

            IF (sw_diag(2)%l_up_flux_trop)                              &
                sw_diag(1)%up_flux_trop =                               &
                sw_diag(1)%up_flux_trop -                               &
                sw_diag(2)%up_flux_trop

            IF (sw_diag(2)%l_flux_direct)                               &
                sw_diag(1)%flux_direct =                                &
                sw_diag(1)%flux_direct -                                &
                sw_diag(2)%flux_direct

            IF (sw_diag(2)%l_flux_diffuse)                              &
                sw_diag(1)%flux_diffuse =                               &
                sw_diag(1)%flux_diffuse -                               &
                sw_diag(2)%flux_diffuse

            IF (sw_diag(2)%l_flux_direct_sph)                           &
                sw_diag(1)%flux_direct_sph =                            &
                sw_diag(1)%flux_direct_sph -                            &
                sw_diag(2)%flux_direct_sph

            IF (sw_diag(2)%l_flux_direct_div)                           &
                sw_diag(1)%flux_direct_div =                            &
                sw_diag(1)%flux_direct_div -                            &
                sw_diag(2)%flux_direct_div

            IF (sw_diag(2)%l_FlxSolBelow690nmSurf)                      &
                sw_diag(1)%FlxSolBelow690nmSurf =                       &
                sw_diag(1)%FlxSolBelow690nmSurf -                       &
                sw_diag(2)%FlxSolBelow690nmSurf

            IF (sw_diag(2)%l_FlxSeaBelow690nmSurf)                      &
                sw_diag(1)%FlxSeaBelow690nmSurf =                       &
                sw_diag(1)%FlxSeaBelow690nmSurf -                       &
                sw_diag(2)%FlxSeaBelow690nmSurf

          END IF

          IF (l_rad_step_diag) THEN

            IF (sw_diag(2)%l_flux_up) THEN
              sw_diag(2)%flux_up =                                      &
              sw_diag(2)%flux_up +                                      &
              sw_diag(1)%flux_up
              WHERE (sw_diag(2)%flux_up < 0.0)                          &
                     sw_diag(2)%flux_up = 0.0
            END IF

            IF (sw_diag(2)%l_flux_down) THEN
              sw_diag(2)%flux_down =                                    &
              sw_diag(2)%flux_down +                                    &
              sw_diag(1)%flux_down
              WHERE (sw_diag(2)%flux_down < 0.0)                        &
                     sw_diag(2)%flux_down = 0.0
            END IF

            IF (sw_diag(2)%l_solar_out_toa) THEN
              sw_diag(2)%solar_out_toa =                                &
              sw_diag(2)%solar_out_toa +                                &
              sw_diag(1)%solar_out_toa
              WHERE (sw_diag(2)%solar_out_toa < 0.0)                    &
                     sw_diag(2)%solar_out_toa = 0.0
            END IF

            IF (sw_diag(2)%l_surface_down_flux) THEN
              sw_diag(2)%surface_down_flux =                            &
              sw_diag(2)%surface_down_flux +                            &
              sw_diag(1)%surface_down_flux
              WHERE (sw_diag(2)%surface_down_flux < 0.0)                &
                     sw_diag(2)%surface_down_flux = 0.0
            END IF

            IF (sw_diag(2)%l_net_flux_trop) THEN
              sw_diag(2)%net_flux_trop =                                &
              sw_diag(2)%net_flux_trop +                                &
              sw_diag(1)%net_flux_trop
              WHERE (sw_diag(2)%net_flux_trop < 0.0)                    &
                     sw_diag(2)%net_flux_trop = 0.0
            END IF

            IF (sw_diag(2)%l_up_flux_trop) THEN
              sw_diag(2)%up_flux_trop =                                 &
              sw_diag(2)%up_flux_trop +                                 &
              sw_diag(1)%up_flux_trop
              WHERE (sw_diag(2)%up_flux_trop < 0.0)                     &
                     sw_diag(2)%up_flux_trop = 0.0
            END IF

            IF (sw_diag(2)%l_flux_direct) THEN
              sw_diag(2)%flux_direct =                                  &
              sw_diag(2)%flux_direct +                                  &
              sw_diag(1)%flux_direct
              WHERE (sw_diag(2)%flux_direct < 0.0)                      &
                     sw_diag(2)%flux_direct = 0.0
            END IF

            IF (sw_diag(2)%l_flux_diffuse) THEN
              sw_diag(2)%flux_diffuse =                                 &
              sw_diag(2)%flux_diffuse +                                 &
              sw_diag(1)%flux_diffuse
              WHERE (sw_diag(2)%flux_diffuse < 0.0)                     &
                     sw_diag(2)%flux_diffuse = 0.0
            END IF

            IF (sw_diag(2)%l_flux_direct_sph) THEN
              sw_diag(2)%flux_direct_sph =                              &
              sw_diag(2)%flux_direct_sph +                              &
              sw_diag(1)%flux_direct_sph
            END IF

            IF (sw_diag(2)%l_flux_direct_div) THEN
              sw_diag(2)%flux_direct_div =                              &
              sw_diag(2)%flux_direct_div +                              &
              sw_diag(1)%flux_direct_div
              WHERE (sw_diag(2)%flux_direct_div < 0.0)                  &
                     sw_diag(2)%flux_direct_div = 0.0
            END IF

            IF (sw_diag(2)%l_FlxSolBelow690nmSurf) THEN
              sw_diag(2)%FlxSolBelow690nmSurf =                         &
              sw_diag(2)%FlxSolBelow690nmSurf +                         &
              sw_diag(1)%FlxSolBelow690nmSurf
              WHERE (sw_diag(2)%FlxSolBelow690nmSurf < 0.0)             &
                     sw_diag(2)%FlxSolBelow690nmSurf = 0.0
            END IF

            IF (sw_diag(2)%l_FlxSeaBelow690nmSurf) THEN
              sw_diag(2)%FlxSeaBelow690nmSurf =                         &
              sw_diag(2)%FlxSeaBelow690nmSurf +                         &
              sw_diag(1)%FlxSeaBelow690nmSurf
              WHERE (sw_diag(2)%FlxSeaBelow690nmSurf < 0.0)             &
                     sw_diag(2)%FlxSeaBelow690nmSurf = 0.0
            END IF

          END IF ! l_rad_step_diag

        END IF ! j_sw == n_swcall


      ELSE IF (l_forcing .AND. l_rad_step_diag) THEN
        i_off=diagnostic_offset*(j_sw-1)
      ELSE IF (l_radiance) THEN
        IF (j_sw == 1) THEN
          i_off = 0
        ELSE
          i_off=90+(diagnostic_offset/20)*(j_sw-1)
        END IF
      ELSE
        i_off=0
      END IF

      IF (model_type /= mt_single_column) THEN
        CALL diagnostics_rad(                                                  &
          sw_sect, j_sw, sw_diag(j_sw), sw_spectrum(j_sw),                     &
          row_length, rows, ozone_levels, cloud_levels, ntiles,                &
          at_extremity, i_off,                                                 &
          ! Fields common to SW and LW radiation
          t_n, t_inc, q_n, qcl_n, cf_n, cfl_n,                                 &
          t_latest, q_latest, qcl_latest, cf_latest, cfl_latest,               &
          t_incr_diagnostic, surfsw, swsea, trans_solid_angle,                 &
          sw_down_sice_weighted_cat, sw_down_sice_weighted,                    &
          ! SW diagnostic fields
          itoasw, surfsw_cor, toasw_cor, surfdir_cor, surfdif_cor,             &
          flux_below_690nm_surf, photosynth_act_rad, flxdirparsurf,            &
          f_orog, sw_net_land, sw_net_sice, sea_salt_film, sea_salt_jet,       &
          sw_up_sice_weighted_cat, sw_up_sice_weighted,                        &
          albedo_sice_weighted_cat, albedo_sice_weighted,                      &
          salt_dim1, salt_dim2, salt_dim3,                                     &
          cos_zenith_angle, day_fraction,                                      &
          cos_zen_rts, day_frac_rts, sol_azm_rts,                              &
          ! LW diagnostic fields (dummys here)
          d_mass, density, layer_heat_capacity,                                &
          p_layer_boundaries, p_layer_centres, p_extra_layer,                  &
          t_layer_boundaries, t_extra_layer,                                   &
          olr, lw_down, ozone, o3_trop_level, o3_trop_height,                  &
          t_trop_level, t_trop_height,                                         &
          stashwork1)
      END IF ! model_type
    END IF ! TmpLogic
  END DO

  IF ( l_scmdiags(scmdiag_rad) .AND.                                    &
       model_type == mt_single_column) THEN

    ! Output some SCM diagnostics for SW radiation
    ! Equiv to stash 1,161
    CALL scmoutput(t_incr_diagnostic,'dt_sw',                           &
         'SW heating rate','K/timestep',                                &
         t_acc,d_all,default_streams,'',routinename)

    CALL scmoutput(t_incr_diagnostic, 'sw2',                            &
         'SW heating rate','K/day',                                     &
         t_mult,d_all,default_streams,'ntspday',routinename)

    CALL scmoutput(surfsw,'surfsw',                                     &
         'Net surface SW flux','W/m2',                                  &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    IF (l_surfsw_cor)                                                   &
      CALL scmoutput(surfsw_cor,'surfsw_cor',                           &
           'Corrected net surface SW flux','W/m2',                      &
           t_avg,d_sl,default_streams,'',routinename)

    IF (l_toasw_cor)                                                    &
         CALL scmoutput(toasw_cor,'toasw_cor',                          &
             'Corrected outgoing SW flux (TOA)','W/m2',                 &
             t_avg,d_sl,default_streams,'',routinename)

    IF (l_surfdir_cor)                                                  &
         CALL scmoutput(surfdir_cor,'surfdir_cor',                      &
             'Corrected direct surface SW flux','W/m2',                 &
             t_avg,d_sl,default_streams,'',routinename)

    IF (l_surfdif_cor)                                                  &
         CALL scmoutput(surfdif_cor,'surfdif_cor',                      &
             'Corrected diffuse surface SW flux','W/m2',                &
             t_avg,d_sl,default_streams,'',routinename)

    CALL scmoutput(itoasw,'is_toa',                                     &
         'Incoming solar radn (TOA)','W/m2',                            &
         t_avg,d_sl,default_streams,'',routinename)

    CALL scmoutput(flux_below_690nm_surf,'surf_sw_b1',                  &
         'Net SW surface flux in band 1','W/m2',                        &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    IF (l_timestep) j_sw=2
    IF (.NOT. l_timestep) j_sw=1

    IF ( ALLOCATED(sw_diag(j_sw)%solar_out_toa ) ) THEN
      TmpScm2d(:,:) = sw_diag(j_sw)%solar_out_toa(:,:)
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF

    ! stash 1,208
    CALL scmoutput(TmpScm2d,'os_toa',                                   &
         'Outgoing solar radn (TOA)','W/m2',                            &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    ! stash 1,182  Vapour increment
    TmpScm3d_1(:,:,:) = q_latest(:,:,:) - q_n(:,:,:)

    CALL scmoutput(TmpScm3d_1,'dq_sw',                                  &
         'Specific humidity increment swrad','kg/kg',                   &
         t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    ! stash 1,183  liquid water content increment
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_1(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_1,'dqcl_sw',                                &
         'qcl increment swrad','kg/kg',                                 &
         t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    ! stash 1,192  total cloud fraction increment
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_1(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_1,'dbcf_sw',                                &
        'bulk cloud fraction increment swrad','fraction',               &
        t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    ! stash 1,193  liquid cloud fraction increment
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_1(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_1,'dcfl_sw',                                &
        'liquid cloud fraction increment swrad','fraction',             &
        t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    ! stash 1,232
    ! SW heating = sw rad temp inc /timestep/timestep
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_2(i,j,k) = t_incr_diagnostic(i,j,k) / timestep
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_2,'sw1rate',                                &
         'SW heating rate all timesteps','K/ts/ts',                     &
        t_avg,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(sw_diag(j_sw)%surface_down_flux ) ) THEN
      TmpScm2d(:,:) = sw_diag(j_sw)%surface_down_flux(:,:)
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF
    ! stash 1,235
    CALL scmoutput(TmpScm2d,'surfsw_dn',                                &
         'Total down SW flux','W/m2',                                   &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)


    IF (l_rad_perturb) j_sw=1

    ! stash 1,209
    IF (ALLOCATED(sw_diag(j_sw)%solar_out_clear)) THEN
      TmpScm2d(:,:) = sw_diag(j_sw)%solar_out_clear
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm2d,'cs_os',               &
         'Clear-sky outgoing SW','W/m2',                                &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    ! stash 1,210
    IF (ALLOCATED(sw_diag(j_sw)%surf_down_clr)) THEN
      TmpScm2d(:,:) = sw_diag(j_sw)%surf_down_clr
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm2d,'cs_surf_dnsw',          &
         'Clear-sky down SW flux','W/m2',                               &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    ! stash 1,211
    IF (ALLOCATED(sw_diag(j_sw)%surf_up_clr)) THEN
      TmpScm2d(:,:) = sw_diag(j_sw)%surf_up_clr
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm2d,'cs_surf_upsw',            &
         'Clear-sky up SW flux','W/m2',                                 &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    ! stash 1,221
    IF (ALLOCATED(sw_diag(j_sw)%re_strat)) THEN
      TmpScm3d_1(:,:,:) = sw_diag(j_sw)%re_strat
    ELSE
      TmpScm3d_1(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_1,'re_strat',                   &
         'Layer cld liq effective radius * layer cld weight','-',       &
         t_avg+only_radsteps,d_cloud,default_streams,'',routinename)

    ! stash 1,223
    IF (ALLOCATED(sw_diag(j_sw)%wgt_strat)) THEN
      TmpScm3d_1(:,:,:) = sw_diag(j_sw)%wgt_strat
    ELSE
      TmpScm3d_1(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_1,'wgt_strat',                 &
         'Layer cloud weight for microphysics','-',                     &
         t_avg+only_radsteps,d_cloud,default_streams,'',routinename)

    ! stash 1,224
    IF (ALLOCATED(sw_diag(j_sw)%lwp_strat)) THEN
      TmpScm3d_1(:,:,:) = sw_diag(j_sw)%lwp_strat
    ELSE
      TmpScm3d_1(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_1,'lwp_strat',                 &
         'Layer cld liquid water path * weight','kg/m2',                &
         t_avg+only_radsteps,d_cloud,default_streams,'',routinename)

    ! stash 1,225
    IF (ALLOCATED(sw_diag(j_sw)%re_conv)) THEN
      TmpScm3d_1(:,:,:) = sw_diag(j_sw)%re_conv
    ELSE
      TmpScm3d_1(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_1,'re_conv',                     &
         'Conv cloud liq re * conv cld weight','-',                     &
         t_avg+only_radsteps,d_cloud,default_streams,'',routinename)

    ! stash 1,226
    IF (ALLOCATED(sw_diag(j_sw)%wgt_conv)) THEN
      TmpScm3d_1(:,:,:) = sw_diag(j_sw)%wgt_conv
    ELSE
      TmpScm3d_1(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_1,'wgt_conv',                   &
         'Conv cloud weight for microphysics','-',                      &
         t_avg+only_radsteps,d_cloud,default_streams,'',routinename)

    ! stash 1,233
    IF (ALLOCATED(sw_diag(j_sw)%clear_hr)) THEN
      TmpScm3d_1(:,:,:) = sw_diag(j_sw)%clear_hr
    ELSE
      TmpScm3d_1(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_1,'dt_cssw',                    &
         'Clear-sky SW heating rates','K/s',                            &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

  END IF ! scmdiag_rad / model_type


  DEALLOCATE ( albedo_sice_weighted_cat  )
  DEALLOCATE ( albedo_sice_weighted      )
  DEALLOCATE ( sw_down_sice_weighted_cat )
  DEALLOCATE ( sw_down_sice_weighted     )
  DEALLOCATE ( sw_up_sice_weighted_cat   )
  DEALLOCATE ( sw_up_sice_weighted       )

  DEALLOCATE ( t_incr_diagnostic )
  DEALLOCATE ( open_sea_albedo )
  DEALLOCATE ( f_orog )

  IF (MOD(timestep_number,a_sw_radstep_prog) == 0) THEN
    IF (ALLOCATED(surf_down_sw_rts)) DEALLOCATE(surf_down_sw_rts)
    IF (ALLOCATED(alb_tile))         DEALLOCATE(alb_tile)
    IF (l_orog) THEN
      DEALLOCATE(orog_corr)
    END IF
  END IF

  IF (l_timestep) THEN
    IF (l_rad_step_diag) THEN
      CALL deallocate_diag(sw_diag(2))
    END IF
    IF (timestep_number >  1) THEN
      IF (MOD(timestep_number,a_sw_radstep_prog) == 0) THEN

        DEALLOCATE(sw_incs_local)
        DEALLOCATE(swsea_local)
        DEALLOCATE(netsw_local)
        DEALLOCATE(top_abs_sw)
        DEALLOCATE(flux_b690nm_local)
        DEALLOCATE(surf_down_sw_local)

        CALL deallocate_diag(sw_diag(1))

      END IF
    END IF
  ELSE IF (l_forcing) THEN
    IF (l_rad_step_prog) THEN
      CALL deallocate_diag(sw_diag(1))
    END IF
    IF (l_rad_step_diag) THEN
      CALL deallocate_diag(sw_diag(2))
    END IF
  ELSE IF (l_radiance) THEN
    IF (l_rad_step_prog) THEN
      DO j_sw = 1, n_swcall
        CALL deallocate_diag(sw_diag(j_sw))
      END DO
    END IF
  ELSE
    IF (l_rad_step_prog) THEN
      CALL deallocate_diag(sw_diag(1))
    END IF
  END IF

  ! ----------------------------------------------------------------------
  ! Section RAD.2 Long Wave Radiation Code.
  ! -----------------------------------------------------------------------

  IF (l_rad_step_prog) THEN

    CALL allocate_diag(lw_diag(1), lw_spectrum(1),                      &
      row_length, rows, model_levels, cloud_levels, ntiles)

    IF (l_timestep) THEN
      ALLOCATE(lw_incs_local(row_length, rows, 0:model_levels, n_lwcall))
      ALLOCATE(lwsea_local(row_length, rows, n_lwcall))
      ALLOCATE(olr_local(row_length, rows, n_lwcall))
      ALLOCATE(lw_down_local(row_length, rows, n_lwcall))
      ALLOCATE(top_abs_lw(row_length, rows, n_lwcall))
    END IF

    IF (l_radiance) THEN
      DO j_lw = 2, n_lwcall
        CALL allocate_diag(lw_diag(j_lw), lw_spectrum(j_lw),            &
          row_length, rows, model_levels, cloud_levels, ntiles)
      END DO
    END IF

  END IF

  IF (l_rad_step_diag) THEN
    IF (l_timestep .OR. l_forcing) THEN
      CALL allocate_diag(lw_diag(2), lw_spectrum(2),                    &
        row_length, rows, model_levels, cloud_levels, ntiles)
    END IF
  END IF

  IF (ltimer) CALL timer ('AP1R LW Rad  ',5)

  IF (l_rad_step_prog .OR. l_rad_step_diag) THEN

    IF (l_planet_grey_surface) THEN
      t_rad_surf = t_surf
      surf_emission = planet_emissivity*t_surf**4
    ELSE

      ! Calculate JULES mean land emissivity for use in t_rad_surf and
      ! dOLR calculations and to pass down to R2_SET_SURFACE_FIELD
      ! Note that, unlike the albedo in the SW, the emissivity is not
      ! recalculated from frac_control on a diagnsotic call. Given the
      ! connection between surface temperatures and emissivity in the LW,
      ! it is less clear how such surface forcing might be calculated.
      ! Since also the impact would be much smaller, no provision is
      ! made for calculation of forcing due to land-use changes in the LW.

      ! For the following logic recall the definition of the mappings:
      ! [tile,type]_index goes to the appropriate index over land points
      ! and land_index then maps that to the grid-point.

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, surf_emission, t_rad_surf, t_rad_land,        &
!$OMP&         emis_land )                                                     &
!$OMP& PRIVATE( i, j )
      DO j = 1, rows
        DO i = 1, row_length
          surf_emission(i,j) = 0.0
          t_rad_surf(i,j)    = 0.0
          t_rad_land(i,j)    = 0.0
          emis_land(i,j)     = 0.0
        END DO
      END DO
!$OMP END PARALLEL DO

      IF (l_aggregate) THEN

        ! First aggregate the emissivity under snow-free conditions.
        DO n=1, npft
          DO point = 1, type_pts(n)
            l = type_index(point, n)
            i = land_index_i(l)
            j = land_index_j(l)
            emis_land(i,j) = emis_land(i,j) +                           &
              frac(l,n) * emis_pft(n)
          END DO
        END DO
        DO n=1, nnvg
          DO point = 1, type_pts(n+npft)
            l = type_index(point, n+npft)
            i = land_index_i(l)
            j = land_index_j(l)
            emis_land(i,j) = emis_land(i,j) +                           &
              frac(l,n+npft) * emis_nvg(n)
          END DO
        END DO

        ! Since JULES does not currently allow for the emissivity
        ! of snow, the surface emission is calculated with the
        ! snow-free emissivity.

        DO point = 1, tile_pts(1)
          l = tile_index(point, 1)
          i = land_index_i(l)
          j = land_index_j(l)
          IF (l_dolr_land_black) THEN
            surf_emission(i,j) = surf_emission(i,j) +                   &
              tstar_tile(l,1)**4
          ELSE
            surf_emission(i,j) = surf_emission(i,j) +                   &
              emis_land(i,j) * tstar_tile(l,1)**4
          END IF
          IF (l_t_land_nosnow) THEN
            emis_nosnow = emis_land(i,j)
          END IF
          IF (l_rad_snow_emis) THEN
            emis_land(i,j) = emis_land(i,j) +                           &
              (emis_nvg(ice-npft) - emis_land(i,j)) *                   &
              ( MAX(0.0, snow_tile(l,1)) /                              &
              ( MAX(0.0, snow_tile(l,1)) +                              &
              10.0 * z0_tile(l,1) * rho_snow_const ) )
          END IF
          ! Now temporarily over-write tstar_tile to replicate
          ! existing science.
          tstar_tile(l,1) = t_surf(i,j)
          IF (l_t_land_nosnow) THEN
            t_rad_land(i,j) = t_rad_land(i,j) +                         &
              emis_nosnow * tstar_tile(l,1)**4
          ELSE
            t_rad_land(i,j) = t_rad_land(i,j) +                         &
              emis_land(i,j) * tstar_tile(l,1)**4
          END IF
        END DO

      ELSE

        ! Types match tiles here.
        emis_tiles(1:npft) = emis_pft
        emis_tiles(npft+1:ntype) = emis_nvg
        DO n = 1, ntype
          DO point = 1, tile_pts(n)
            l = type_index(point, n)
            i = land_index_i(l)
            j = land_index_j(l)
            ! Calculate the surface emission without adjusting for
            ! the emissivity of snow to match what JULES will add
            ! back.
            emis_here = emis_tiles(n)
            IF (l_dolr_land_black) THEN
              surf_emission(i,j) = surf_emission(i,j) +                 &
                tile_frac(l,n) * tstar_tile(l,n)**4
            ELSE
              surf_emission(i,j) = surf_emission(i,j) +                 &
                tile_frac(l,n) * emis_here * tstar_tile(l,n)**4
            END IF
            ! Adjust for snow and accumulate into calculation
            ! of the radiative temperature.
            IF (l_t_land_nosnow) THEN
              emis_nosnow = emis_here
            END IF
            IF (l_rad_snow_emis) THEN
              emis_here = emis_here +                                   &
                (emis_nvg(ice-npft) - emis_here) *                      &
                ( MAX(0.0, snow_tile(l,n)) /                            &
                ( MAX(0.0, snow_tile(l,n)) +                            &
                10.0 * z0_tile(l,n) * rho_snow_const ) )
            END IF
            emis_land(i,j) = emis_land(i,j) +                           &
              frac(l, n) * emis_here
            IF (l_t_land_nosnow) THEN
              t_rad_land(i,j) = t_rad_land(i,j) +                       &
                tile_frac(l,n) * emis_nosnow * tstar_tile(l,n)**4
            ELSE
              t_rad_land(i,j) = t_rad_land(i,j) +                       &
                tile_frac(l,n) * emis_here * tstar_tile(l,n)**4
            END IF
          END DO
        END DO

      END IF

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, flandg, surf_emission, t_rad_surf,            &
!$OMP&         t_rad_land, emis_land, t_rad_sice )                             &
!$OMP& PRIVATE( i, j )
      DO j = 1, rows
        DO i = 1, row_length
          IF (flandg(i,j) > 0.0) THEN
            surf_emission(i,j) = flandg(i,j) * surf_emission(i,j)
            t_rad_surf(i,j) =  t_rad_surf(i,j) + flandg(i,j) * t_rad_land(i,j)
            t_rad_land(i,j) = SQRT( SQRT( t_rad_land(i,j) / emis_land(i,j) ) )
          END IF
          t_rad_sice(i,j) = 0.0
        END DO
      END DO
!$OMP END PARALLEL DO

      ! Effective surface radiative temperature over sea ice
      t_rad_sice(:,:) = 0.0

      DO n = 1, nice_use
        DO point = 1, sice_pts_ncat(n)
          l = sice_index_ncat(point,n)
          i = ssi_index_i(l)
          j = ssi_index_j(l)
          t_rad_sice(i,j) = t_rad_sice(i,j) +                           &
                 ice_fract_cat(i,j,n) * emis_sice *                     &
                 tstar_sice_cat(i,j,n)**4
        END DO
      END DO

      ! The test for the presence of sea ice in atmos_physics1 is
      ! whether ice_fract exceeds 0.
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( rows, row_length, ice_fract, t_rad_surf, flandg,                &
!$OMP&         t_rad_sice, surf_emission, emis_sice, emis_sea,                 &
!$OMP&         tstar_sea, emis_land, t_rad_solid )                             &
!$OMP& PRIVATE( i, j )
      DO j = 1, rows
        DO i = 1, row_length
          IF (ice_fract(i,j) > 0.0) THEN
            t_rad_surf(i,j) =  t_rad_surf(i,j) + (1.0 - flandg(i,j)) *         &
                               t_rad_sice(i,j)
            surf_emission(i,j) =  surf_emission(i,j) + (1.0 -                  &
                                  flandg(i,j)) * t_rad_sice(i,j)
            t_rad_sice(i,j) = SQRT( SQRT( t_rad_sice(i,j) /                    &
                                 (emis_sice * ice_fract(i,j) ) ) )
          END IF

      ! Contributions from open sea
          IF (flandg(i,j) < 1.0) THEN
            t_rad_surf(i,j) =  t_rad_surf(i,j) +                               &
                         (1.0 - flandg(i,j)) * (1.0 - ice_fract(i,j)) *        &
                         emis_sea * tstar_sea(i,j)**4
            surf_emission(i,j) =  surf_emission(i,j) +                         &
                         (1.0 - flandg(i,j)) * (1.0 - ice_fract(i,j)) *        &
                         emis_sea * tstar_sea(i,j)**4
          END IF

          t_rad_surf(i,j) = SQRT( SQRT( t_rad_surf(i,j) / (flandg(i,j)*        &
                                        emis_land(i,j) +                       &
                       (1.0 - flandg(i,j)) * ( ice_fract(i,j) *                &
                                        emis_sice +                            &
                       (1.0 - ice_fract(i,j)) * emis_sea) ) ) )
          t_rad_solid(i,j) = 0.0
        END DO
      END DO
!$OMP END PARALLEL DO

      IF ( l_ctile .AND. .NOT. l_quad_t_coast ) THEN
        ! Use linear averaging at coastal points for consistency.
        WHERE ( (flandg > 0) .AND. (.NOT. land0p5) )
          t_rad_surf = (1.0 - flandg) * (1.0 - ice_fract) * tstar_sea
        END WHERE
        DO n = 1, nice_use
          DO point = 1, sice_pts_ncat(n)
            l = sice_index_ncat(point,n)
            i = ssi_index_i(l)
            j = ssi_index_j(l)
            IF ( (flandg(i,j) > 0) .AND. (.NOT. land0p5(i,j)) )         &
              t_rad_surf(i,j) = t_rad_surf(i,j) +                       &
                (1.0 - flandg(i,j)) *                                   &
                ice_fract_cat(i,j,n) * tstar_sice_cat(i,j,n)
          END DO
        END DO
        DO n = 1, ntype
          DO point = 1, tile_pts(n)
            l = type_index(point, n)
            i = land_index_i(l)
            j = land_index_j(l)
            IF ( (flandg(i,j) > 0) .AND. (.NOT. land0p5(i,j)) )         &
              t_rad_surf(i,j) = t_rad_surf(i,j) +                       &
                flandg(i,j) * tile_frac(l,n) * tstar_tile(l,n)
          END DO
        END DO
      END IF

      IF ( l_t_rad_solid ) THEN
        ! Use deprecated common solid temperature for land and sea ice
        ! Only coastal points require modification.
        DO n = 1, nice_use
          DO point = 1, sice_pts_ncat(n)
            l = sice_index_ncat(point,n)
            i = ssi_index_i(l)
            j = ssi_index_j(l)
            IF (flandg(i,j) > 0)                                        &
              t_rad_solid(i,j) = t_rad_solid(i,j) +                     &
                (1.0 - flandg(i,j)) * ice_fract_cat(i,j,n) *            &
                emis_sice *tstar_sice_cat(i,j,n)**4
          END DO
        END DO

!$OMP  PARALLEL DO DEFAULT(NONE)                                               &
!$OMP& SHARED( rows, row_length, flandg, t_rad_solid, emis_land,               &
!$OMP&         t_rad_land, ice_fract, emis_sice, t_rad_sice )                  &
!$OMP& PRIVATE( i, j )
        DO j = 1, rows
          DO i = 1, row_length
            IF ( (flandg(i,j) > 0) .AND. (flandg(i,j) < 1) ) THEN
              t_rad_solid(i,j) = t_rad_solid(i,j) +                     &
                flandg(i,j) * emis_land(i,j) * t_rad_land(i,j)**4
              t_rad_solid(i,j) = SQRT( SQRT( t_rad_solid(i,j) /         &
                ( (1.0 - flandg(i,j)) * ice_fract(i,j) * emis_sice +    &
                flandg(i,j) * emis_land(i,j) ) ) )
              ! Over-write.
              t_rad_land(i,j) = t_rad_solid(i,j)
              t_rad_sice(i,j) = t_rad_solid(i,j)
            END IF
          END DO
        END DO
!$OMP END PARALLEL DO
      END IF

    END IF ! l_planet_grey_surface

    IF (l_t_bdy_surf) THEN
      ! Take the temperature of the air just above the surface as
      ! the temperature at the surface.
      DO j = 1, rows
        DO i = 1, row_length
          t_layer_boundaries(i,j,0) = t_rad_surf(i,j)
        END DO
      END DO
    END IF

    DO j_lw = n_lwcall, 1, -1
      ! Decide whether it is necessary to call radiation
      l_call_lwrad=.FALSE.
      IF (j_lw==1 .AND. l_rad_step_prog) THEN
        l_call_lwrad=.TRUE.
      ELSE IF (l_timestep .AND. j_lw == 2) THEN
        l_call_lwrad=.TRUE.
      ELSE IF (l_forcing .AND. j_lw > 1) THEN
        IF (lw_diag(j_lw)%l_diag_call) THEN
          l_call_lwrad=.TRUE.
        END IF
      END IF

      IF (l_call_lwrad) THEN

        ! Set COSP flag only on prognostic radiation steps
        IF (l_cosp_in .AND. (j_lw == 1) .AND. l_rad_step_prog) THEN
          l_cosp=.TRUE.
        ELSE
          l_cosp=.FALSE.
        END IF

        ! For each LW call set the correct mass mixing ratios and switch the
        ! aerosols on/off, &c.  Note that as for the SW, in the case of the
        ! aerosols the reference  state is 'no aerosols' and hence they are
        ! switched off in the diagnostic call.

        ! As in the SW, since the "default" for each variable is to have it
        ! set in the diagnostic call to the same as in the prognostic, we set
        ! them all unconditionally to the prognostic & then re-set as needed.

        ! Set default values for greenhouse gases
        j_co2_mmr    = co2_mmr
        j_n2o_mmr    = n2ommr
        j_ch4_mmr    = ch4mmr
        j_so2_mmr    = so2mmr
        j_cfc11_mmr  = c11mmr
        j_cfc12_mmr  = c12mmr
        j_c113_mmr   = c113mmr
        j_c114_mmr   = c114mmr
        j_hcfc22_mmr = hcfc22mmr
        j_hfc125_mmr = hfc125mmr
        j_hfc134_mmr = hfc134ammr
        j_ozone      = ozone

        ! Set default flags for aerosols
        j_l_sulpc_so2          = l_sulpc_so2
        j_l_use_sulpc_direct   = l_use_sulpc_direct
        j_l_use_seasalt_direct = l_use_seasalt_direct
        j_l_soot               = l_soot
        j_l_use_soot_direct    = l_use_soot_direct
        j_l_biomass            = l_biomass
        j_l_use_bmass_direct   = l_use_bmass_direct
        j_l_ocff               = l_ocff
        j_l_use_ocff_direct    = l_use_ocff_direct
        j_l_nitrate            = l_nitrate
        j_l_use_nitrate_direct = l_use_nitrate_direct
        j_l_dust               = l_dust
        j_l_use_dust           = l_use_dust
        j_l_use_biogenic       = l_use_biogenic
        j_l_climat_aerosol     = l_climat_aerosol
        j_l_murk_rad           = l_murk_rad
        j_l_use_arcl           = l_use_arcl
        j_n_arcl_species       = n_arcl_species
        j_l_use_ukca_radaer    = l_ukca_radaer
        j_l_use_glomap_clim_radaer = l_glomap_clim_radaer
        j_l_use_easyaerosol    = l_easyaerosol_lw

        ! If timestepping required
        IF (l_timestep .AND. j_lw == 2) THEN

          ! Turn aerosols off for the "cloud only" call in order to minimise
          ! computational cost. Clear-sky and other diagnostics have already
          ! been turned off in set_lwdiag_logic.
          IF (l_rad_perturb) THEN
            j_l_sulpc_so2          = .FALSE.
            j_l_use_sulpc_direct   = .FALSE.
            j_l_use_seasalt_direct = .FALSE.
            j_l_soot               = .FALSE.
            j_l_use_soot_direct    = .FALSE.
            j_l_biomass            = .FALSE.
            j_l_use_bmass_direct   = .FALSE.
            j_l_ocff               = .FALSE.
            j_l_use_ocff_direct    = .FALSE.
            j_l_nitrate            = .FALSE.
            j_l_use_nitrate_direct = .FALSE.
            j_l_dust               = .FALSE.
            j_l_use_dust           = .FALSE.
            j_l_use_biogenic       = .FALSE.
            j_l_climat_aerosol     = .FALSE.
            j_l_murk_rad           = .FALSE.
            j_l_use_arcl           = .FALSE.
            j_n_arcl_species       = 0
            j_l_use_ukca_radaer    = .FALSE.
            j_l_use_glomap_clim_radaer = .FALSE.
            j_l_use_easyaerosol    = .FALSE.
          END IF

        ! If diagnostic call required
        ELSE IF (l_forcing .AND. j_lw > 1) THEN

          ! Greenhouse Gases
          IF ( c2c_co2 ) THEN
            j_co2_mmr = co2_mmr*co2_mmr_scl + co2_mmr_add
          END IF
          IF ( c2c_n2o ) THEN
            j_n2o_mmr = n2ommr*n2o_mmr_scl + n2o_mmr_add
            IF (ngrgas >= p_n2o) THEN
              ALLOCATE(n2o_mmr_tmp(row_length, rows, model_levels))
              n2o_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_n2o)
              grgas_field(:,:,:,p_n2o) = &
                grgas_field(:,:,:,p_n2o)*n2o_mmr_scl + n2o_mmr_add
            END IF
          END IF
          IF ( c2c_ch4 ) THEN
            j_ch4_mmr = ch4mmr*ch4_mmr_scl + ch4_mmr_add
            IF (ngrgas >= p_ch4) THEN
              ALLOCATE(ch4_mmr_tmp(row_length, rows, model_levels))
              ch4_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_ch4)
              grgas_field(:,:,:,p_ch4) = &
                grgas_field(:,:,:,p_ch4)*ch4_mmr_scl + ch4_mmr_add
            END IF
          END IF
          IF ( c2c_cfc11 ) THEN
            j_cfc11_mmr = c11mmr*cfc11_mmr_scl + cfc11_mmr_add
            IF (ngrgas >= p_f11) THEN
              ALLOCATE(cfc11_mmr_tmp(row_length, rows, model_levels))
              cfc11_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_f11)
              grgas_field(:,:,:,p_f11) = &
                grgas_field(:,:,:,p_f11)*cfc11_mmr_scl + cfc11_mmr_add
            END IF
          END IF
          IF ( c2c_cfc12 ) THEN
            j_cfc12_mmr = c12mmr*cfc12_mmr_scl + cfc12_mmr_add
            IF (ngrgas >= p_f12) THEN
              ALLOCATE(cfc12_mmr_tmp(row_length, rows, model_levels))
              cfc12_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_f12)
              grgas_field(:,:,:,p_f12) = &
                grgas_field(:,:,:,p_f12)*cfc12_mmr_scl + cfc12_mmr_add
            END IF
          END IF
          IF ( c2c_c113 ) THEN
            j_c113_mmr = c113mmr*cfc113_mmr_scl + cfc113_mmr_add
            IF (ngrgas >= p_f113) THEN
              ALLOCATE(cfc113_mmr_tmp(row_length, rows, model_levels))
              cfc113_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_f113)
              grgas_field(:,:,:,p_f113) = &
                grgas_field(:,:,:,p_f113)*cfc113_mmr_scl + cfc113_mmr_add
            END IF
          END IF
          IF ( c2c_hcfc22 ) THEN
            j_hcfc22_mmr = hcfc22mmr*hcfc22_mmr_scl + hcfc22_mmr_add
            IF (ngrgas >= p_f22) THEN
              ALLOCATE(hcfc22_mmr_tmp(row_length, rows, model_levels))
              hcfc22_mmr_tmp(:,:,:) = grgas_field(:,:,:,p_f22)
              grgas_field(:,:,:,p_f22) = &
                grgas_field(:,:,:,p_f22)*hcfc22_mmr_scl + hcfc22_mmr_add
            END IF
          END IF
          IF ( c2c_hfc125 ) THEN
            j_hfc125_mmr = hfc125mmr*hfc125_mmr_scl + hfc125_mmr_add
          END IF
          IF ( c2c_hfc134 ) THEN
            j_hfc134_mmr = hfc134ammr*hfc134a_mmr_scl + hfc134a_mmr_add
          END IF

          ! Ozone: Note that aerosol contains the MURK array. Care needs to be
          ! taken so that MURK has the correct dimensions before it arrives in
          ! RAD_CTL.
          IF ( c2c_o3 ) j_ozone = aerosol

          ! Aerosols
          IF ( c2c_sulpc_d ) j_l_sulpc_so2          = .FALSE.
          IF ( c2c_sulpc_d ) j_l_use_sulpc_direct   = .FALSE.
          IF ( c2c_seas_d )  j_l_use_seasalt_direct = .FALSE.
          IF ( c2c_soot_d )  j_l_soot               = .FALSE.
          IF ( c2c_soot_d )  j_l_use_soot_direct    = .FALSE.
          IF ( c2c_bmb_d )   j_l_biomass            = .FALSE.
          IF ( c2c_bmb_d )   j_l_use_bmass_direct   = .FALSE.
          IF ( c2c_ocff_d )  j_l_ocff               = .FALSE.
          IF ( c2c_ocff_d )  j_l_use_ocff_direct    = .FALSE.
          IF ( c2c_nitr_d )  j_l_nitrate            = .FALSE.
          IF ( c2c_nitr_d )  j_l_use_nitrate_direct = .FALSE.
          IF ( c2c_dust_d )  j_l_dust               = .FALSE.
          IF ( c2c_dust_d )  j_l_use_dust           = .FALSE.
          IF ( c2c_biog_d )  j_l_use_biogenic       = .FALSE.
          IF ( c2c_ukca_d )  j_l_use_ukca_radaer    = .FALSE.
          IF ( c2c_ukca_d )  j_l_use_glomap_clim_radaer = .FALSE.
          IF ( c2c_easy_d )  j_l_use_easyaerosol    = .FALSE.
        END IF

        CALL prelim_lwrad(                                              &
        ! Parallel variables
          at_extremity, n_proc,                                         &
        ! Model Dimensions
          row_length,rows,model_levels,                                 &
        ! Model Switches
          l_rad_deg, l_extra_top,                                       &
          lw_control(j_lw)%l_subsample,                                 &
          lw_control(j_lw)%l_geostationary,                             &
        ! Time stepping Information
          timestep_number,a_lw_radstep_prog,                            &
        ! ancillary fields and fields needed to be kept from timestep to
        ! timestep
          lw_incs,                                                      &
        ! Satellite Geometry
          lw_control(j_lw)%min_view_lon,                                &
          lw_control(j_lw)%max_view_lon,                                &
          lw_control(j_lw)%min_view_lat,                                &
          lw_control(j_lw)%max_view_lat,                                &
        ! Number of Call
          j_lw,                                                         &
        ! Other variables
          true_latitude,true_longitude,                                 &
          seconds_since_midnight,                                       &
          rad_mask, list_lw_points, first_data_interp,                  &
          first_row,last_row,diag_row_list(1,j_lw),                     &
          diag_col_list(1,j_lw),                                        &
          lw_points,                                                    &
          olr, lw_down, lwsea, top_absorption )

        CALL set_control(lw_control(j_lw), lw_spectrum(j_lw),                  &
        ! Control flags
          j_l_climat_aerosol, j_l_use_sulpc_direct, j_l_use_soot_direct,       &
          j_l_use_biogenic, j_l_use_dust, j_l_use_bmass_direct,                &
          j_l_use_ocff_direct, j_l_use_nitrate_direct, j_l_use_seasalt_direct, &
          j_l_murk_rad, l_use_sulpc_indirect_sw, j_n_arcl_species,             &
          j_l_use_ukca_radaer, j_l_use_glomap_clim_radaer,                     &
          j_l_use_easyaerosol,                                                 &
        ! Diagnostic options
          lw_diag(j_lw), l_cosp)

        !Number of threads is 1 at this point: we are not load-balancing between
        !threads, then segmenting here.  Rather, set up segmentation on one
        !thread, then distribute the segments among threads.
        num_parallel_lw=1
        ipar=1

        !Set up automatic segment tuning for LW radiation
        IF (l_autotune_segments .AND. .NOT. ALLOCATED(lw_autotune_state)) THEN
          ALLOCATE(lw_autotune_state)
          CALL autotune_init(                   &
            lw_autotune_state,                  &
            name           = 'LW_Radiation',    &
            tag            = 'LW-RAD',          &
            start_size     = a_lw_seg_size,     &
            calls_per_step = n_lwcall)
        END IF

        !If autotuning is active, retrieves a new trial segment size. If not,
        !the returned segment size will stay the same.
        segment_size = a_lw_seg_size
        IF (l_autotune_segments) THEN
          segment_size = autotune_get_trial_size(lw_autotune_state)
        END IF

        ! Set up segmentation information
        CALL segments_mod_seg_meta(meta_segments, ipar, num_parallel_lw,&
          lw_points, segment_size, a_lw_segments)

        ! Allocate space for segmentation arrays
        ALLOCATE( segments( meta_segments%num_segments ) )

        ! Fill the segments array
        CALL segments_mod_segments( segments, meta_segments,            &
               row_length, rows,                                        &
               list_points=list_lw_points)

        IF (l_autotune_segments) THEN
          CALL autotune_start_region(lw_autotune_state)
        END IF

!$OMP  PARALLEL  DEFAULT(SHARED)                                             &
!$OMP& PRIVATE(i, j, first_point, first_point_dust_a,                        &
!$OMP& first_point_dust_b, first_point_sulpc,                                &
!$OMP& first_point_soot,first_point_biomass,                                 &
!$OMP& error_code,                                                           &
!$OMP& first_point_biogenic, first_point_ocff, first_point_arcl,             &
!$OMP& first_point_nitrate, dimen,                                           &
!$OMP& nd_field_flux_diag, nd_field_rad_diag, ptr_local, first_point_ukca,   &
!$OMP& ii, jj)

!$OMP DO SCHEDULE(DYNAMIC)
        DO i = 1, meta_segments%num_segments

          first_point = segments(i)%start_index

          ! Additional work on list_points
          DO j = first_point, segments(i)%end_index
            list_lw_points(j) = list_lw_points(j)-segments(i)%start+1
          END DO

          ! Set the pointer to the beginning of the current segment.
          ! This is done solely to reduce the number of continuation
          ! lines required by the call.
          ptr_local=segments(i)%fp
          ii=segments(i)%first_x
          jj=segments(i)%first_y

          ! Set the first point of the dust arrays to be used.
          IF (l_dust) THEN
            IF (l_twobin_dust) THEN
              first_point_dust_a = segments(i)%fp
              first_point_dust_b = 1
            ELSE
              first_point_dust_a = segments(i)%fp
              first_point_dust_b = segments(i)%fp
            END IF
          ELSE
            first_point_dust_a = 1
            first_point_dust_b = 1
          END IF

          ! Set the first point of the biogenic array
          IF (l_use_biogenic) THEN
            first_point_biogenic=segments(i)%fp
          ELSE
            first_point_biogenic=1
          END IF

          ! Set the first points of the arrays of sulphates to be used.

          ! A separate assignment is necessary since
          ! not be of the full size unless the sulphur cycle is on.
          IF (l_sulpc_so2 .OR. l_use_sulpc_indirect_lw) THEN
            first_point_sulpc=segments(i)%fp
          ELSE
            first_point_sulpc=1
          END IF

          IF (l_soot) THEN
            first_point_soot=segments(i)%fp
          ELSE
            first_point_soot=1
          END IF

          IF (l_biomass .OR. l_use_bmass_indirect) THEN
            first_point_biomass=segments(i)%fp
          ELSE
            first_point_biomass=1
          END IF

          IF (l_ocff .OR. l_use_ocff_indirect) THEN
            first_point_ocff=segments(i)%fp
          ELSE
            first_point_ocff=1
          END IF

          IF (n_arcl_species > 0) THEN
            first_point_arcl = segments(i)%fp
          ELSE
            first_point_arcl = 1
          END IF

          IF (l_nitrate .OR. l_use_nitrate_indirect) THEN
            first_point_nitrate=segments(i)%fp
          ELSE
            first_point_nitrate=1
          END IF

          IF (l_ukca_radaer .OR. l_glomap_clim_radaer) THEN
            first_point_ukca=segments(i)%fp
          ELSE
            first_point_ukca=1
          END IF

          ! Set dynamically determined array sizes
          CALL set_dimen(                                               &
            lw_control(j_lw), dimen, lw_spectrum(j_lw),                 &
            row_length*rows, segments(i)%use_points,                    &
            n_rad_layers, global_cloud_top, n_ukca_mode,                &
            j_l_use_easyaerosol,                                        &
            nd_field_flux_diag, nd_field_rad_diag)


          CALL lw_rad(error_code,                                       &
          ! Input data
            q_n(ii,jj,1), j_co2_mmr, j_ozone(ii,jj,1),                  &
            co2_dim1, co2_dim2, co2_3d(ii,jj,1), l_co2_3d,              &

          ! chemical greenhouse gas fields
            ngrgas, grgas_field(ii,jj,1,1),                             &
            j_n2o_mmr, j_ch4_mmr, j_so2_mmr, j_cfc11_mmr, j_cfc12_mmr,  &
            j_c113_mmr, j_c114_mmr, j_hcfc22_mmr, j_hfc125_mmr,         &
            j_hfc134_mmr, t_n(ii,jj,1),t_rad_surf(ii,jj),               &
            t_rad_land(ii,jj), t_rad_sice(ii,jj),                       &
            tstar_sea(ii,jj),                                           &
            p_star(ii,jj),p_layer_boundaries(ii,jj,0),                  &
            p_layer_centres(ii,jj,0), t_layer_boundaries(ii,jj,0),      &
            p_extra_layer(ii,jj), t_extra_layer(ii,jj),                 &
            d_mass(ii,jj,1), density(ii,jj,1),                          &
            layer_heat_capacity(ii,jj,1),                               &
            r_layer_centres(ii,jj,1), r_layer_boundaries(ii,jj,0),      &

          ! Options for COSP
            l_cosp,                                                     &

          ! Options for treating clouds
            l_inhom_cloud, inhom_cloud_lw, dp_corr_strat, dp_corr_conv, &

          ! Stratiform Cloud Fields
            area_cloud_fraction(ii,jj,1), cf_n(ii,jj,1),                &
            qcl_n(ii,jj,1),qcf_total(ii,jj,1),                          &
            n_drop_pot(ii,jj,1),                                        &

          ! Convective Cloud Fields
            cca(ii,jj,1),cclwp(ii,jj),                                  &
            ccw(ii,jj,1), lcbase(ii,jj),                                &
            ccb(ii,jj), cct(ii,jj),                                     &

          ! Surface Fields. Only want the 0.5 threshold LAND mask
          ! and fractional land:
            land0p5(ii,jj),flandg(ii,jj),                               &
            ice_fract(ii,jj),snow_depth(ii,jj),                         &
            emis_land(ii,jj),                                           &

          ! Solar Fields
            cos_zen_rts(ii,jj), day_frac_rts(ii,jj), solcon_rts, scs, sindec,  &

          ! Aerosol Fields
            j_l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,   &
            zh(ii,jj), aero_bl_levels,                                  &
            j_l_dust, j_l_use_dust, dust_dim1, dust_dim2,               &
            dust_1(first_point_dust_a,1), dust_2(first_point_dust_a,1), &
            dust_3(first_point_dust_b,1), dust_4(first_point_dust_b,1), &
            dust_5(first_point_dust_b,1), dust_6(first_point_dust_b,1), &
            j_l_use_biogenic, biogenic_dim1, biogenic_dim2,             &
            local_biogenic(first_point_biogenic, 1),                    &
            j_l_sulpc_so2, j_l_use_sulpc_direct, l_use_sulpc_indirect_lw,&
            sulp_dim1, sulp_dim2,                                       &
            accum_sulphate(first_point_sulpc, 1),                       &
            aitken_sulphate(first_point_sulpc, 1),                      &
            diss_sulphate(first_point_sulpc, 1),                        &
            sea_salt_film, sea_salt_jet,                                &
            l_use_seasalt_indirect, j_l_use_seasalt_direct,             &
            salt_dim1, salt_dim2, salt_dim3, j_l_soot, j_l_use_soot_direct,&
            soot_dim1, soot_dim2, fresh_soot(first_point_soot, 1),      &
            aged_soot(first_point_soot, 1), j_l_biomass, j_l_use_bmass_direct,&
            bmass_dim1, bmass_dim2, fresh_bmass(first_point_biomass, 1),&
            aged_bmass(first_point_biomass, 1),                         &
            cloud_bmass(first_point_biomass, 1), l_use_bmass_indirect,  &
            j_l_ocff, j_l_use_ocff_direct, ocff_dim1, ocff_dim2,        &
            fresh_ocff(first_point_ocff, 1), aged_ocff(first_point_ocff, 1),&
            cloud_ocff(first_point_ocff, 1), l_use_ocff_indirect,       &
            j_l_nitrate, j_l_use_nitrate_direct, nitrate_dim1, nitrate_dim2,&
            accum_nitrate(first_point_nitrate, 1),                      &
            diss_nitrate(first_point_nitrate, 1), l_use_nitrate_indirect,&
            j_l_use_arcl, arcl_dim1, arcl_dim2, j_n_arcl_species,       &
            n_arcl_compnts, i_arcl_compnts,local_arcl(first_point_arcl,1,1),&
            aerosol(ii,jj,1), j_l_murk_rad,                             &
            j_l_use_ukca_radaer, j_l_use_glomap_clim_radaer,            &
            ukca_radaer, ukca_dim1, ukca_dim2,                          &
            local_ukca_mmr(first_point_ukca, 1, 1),                     &
            local_ukca_cvl(first_point_ukca, 1, 1),                     &
            local_ukca_dry(first_point_ukca, 1, 1),                     &
            local_ukca_wet(first_point_ukca, 1, 1),                     &
            local_ukca_rho(first_point_ukca, 1, 1),                     &
            local_ukca_vol(first_point_ukca, 1, 1),                     &
            local_ukca_wtv(first_point_ukca, 1, 1),                     &
            local_ukca_nbr(first_point_ukca, 1, 1),                     &
            j_l_use_easyaerosol, easyaerosol_lw,                        &
            l_easyaerosol_cdnc, easyaerosol_cdnc,                       &

          ! time
            previous_time, seconds_since_midnight,                      &

          ! grid-dependent arrays
            true_latitude(ii,jj), true_longitude(ii,jj),                &
            obs_solid_angle(ii,jj),                                     &

          ! Level of tropopause
            trindx(ii, jj),                                             &

          ! Spectral data
            lw_spectrum(j_lw),                                          &

          ! Algorithmic options
            lw_control(j_lw), timestep,                                 &
            list_lw_points(first_point),                                &

          ! All diagnostics
            lw_diag(j_lw), diag_row_list(first_point,j_lw),             &
            diag_col_list(first_point,j_lw),                            &

          ! Physical Dimensions
            dimen,                                                      &
            segments(i)%use_points, n_rad_layers, global_cloud_top,     &
            ozone_levels, row_length, rows, row_length*rows,            &
            nd_field_flux_diag, nd_field_rad_diag,                      &
            n_cca_levels, n_ukca_mode, n_ukca_cpnt,                     &

          ! Output data
            olr(ii, jj), lw_down(ii, jj),                               &
            top_absorption(ii, jj), lwsea(ii, jj),                      &
            lw_incs(ii,jj,0),                                           &

          ! COSP input arguments
            cosp_gbx, cosp_sgx, cosp_sgh)

        END DO ! end loop over long-wave segments
!$OMP END DO nowait
!$OMP END PARALLEL

        !If autotuning is active, decide what to do with the
        !segment size and report the current status.
        IF (l_autotune_segments) THEN
          CALL autotune_stop_region(lw_autotune_state, lw_points)
          CALL autotune_advance(lw_autotune_state)
          CALL autotune_report(lw_autotune_state, quiet=(j_lw > 1))
        END IF

        ! Deallocate the segmentation arrays
        DEALLOCATE(segments)

        ! Deallocate band-by-band control options
        CALL deallocate_control(lw_control(j_lw))

        IF (model_type /= mt_single_column) THEN

          ! Radiative fluxes may not have been calculated at all
          ! points. We now fill in as required.

          l_complete_north=.FALSE.
          l_complete_south=.FALSE.

          ! When spatial degradation is performed fields must be
          ! filled in at alternate points.
          l_complete_deg = ( l_rad_deg )

          ! Set addressing limits for spatial degradation.
          IF ( l_complete_deg ) THEN
            first_row=1
            last_row=rows
          END IF

          CALL fill_missing_data_lw(                                    &
            off_x, off_y, row_length, rows, model_levels,               &
            cloud_levels, lw_spectrum(j_lw)%aerosol%n_aod_wavel,        &
            first_row,last_row,                                         &
            first_data_interp, es_space_interp,                         &
            l_complete_north, l_complete_south,                         &
            l_complete_deg, lw_control(j_lw)%n_channel, j_lw,           &
            l_extra_top, lw_incs, olr, lw_down, lwsea, top_absorption )

        END IF ! model_type

        IF (l_rad_perturb .AND. (j_lw==1)) THEN
          lw_incs_local(:,:,:,1) = lw_incs                              &
                                 - lw_incs_local(:,:,:,2)
          olr_local(:,:,1) = olr                                        &
                           - olr_local(:,:,2)
          lw_down_local(:,:,1) = lw_down                                &
                               - lw_down_local(:,:,2)
          lwsea_local(:,:,1) = lwsea                                    &
                             - lwsea_local(:,:,2)
          IF (l_extra_top) THEN
            top_abs_lw(:,:,1) = top_absorption                          &
                              - top_abs_lw(:,:,2)
          END IF
        ELSE IF (l_timestep) THEN
          lw_incs_local(:,:,:,j_lw)=lw_incs
          olr_local(:,:,j_lw)=olr
          lw_down_local(:,:,j_lw)=lw_down
          lwsea_local(:,:,j_lw)=lwsea
          IF (l_extra_top) THEN
            top_abs_lw(:,:,j_lw)=top_absorption
          END IF
        END IF

        IF (l_forcing .AND. j_lw > 1) THEN
          ! Reset grgas_field array to original values
          IF (c2c_n2o .AND. ngrgas >= p_n2o) THEN
            grgas_field(:,:,:,p_n2o) = n2o_mmr_tmp(:,:,:)
            DEALLOCATE(n2o_mmr_tmp)
          END IF
          IF (c2c_ch4 .AND. ngrgas >= p_ch4) THEN
            grgas_field(:,:,:,p_ch4) = ch4_mmr_tmp(:,:,:)
            DEALLOCATE(ch4_mmr_tmp)
          END IF
          IF (c2c_cfc11 .AND. ngrgas >= p_f11) THEN
            grgas_field(:,:,:,p_f11) = cfc11_mmr_tmp(:,:,:)
            DEALLOCATE(cfc11_mmr_tmp)
          END IF
          IF (c2c_cfc12 .AND. ngrgas >= p_f12) THEN
            grgas_field(:,:,:,p_f12) = cfc12_mmr_tmp(:,:,:)
            DEALLOCATE(cfc12_mmr_tmp)
          END IF
          IF (c2c_c113 .AND. ngrgas >= p_f113) THEN
            grgas_field(:,:,:,p_f113) = cfc113_mmr_tmp(:,:,:)
            DEALLOCATE(cfc113_mmr_tmp)
          END IF
          IF (c2c_hcfc22 .AND. ngrgas >= p_f22) THEN
            grgas_field(:,:,:,p_f22) = hcfc22_mmr_tmp(:,:,:)
            DEALLOCATE(hcfc22_mmr_tmp)
          END IF
        END IF

      END IF ! if radiation call required
    END DO ! loop over radiation calls


    IF (l_rad_perturb .AND. l_rad_step_prog) THEN

      ! For this case the full fluxes are already set.

    ELSE IF (l_timestep .AND. l_rad_step_diag) THEN
      lw_incs(:,:,:) = lw_incs_local(:,:,:,1)                           &
                     + lw_incs_local(:,:,:,2)
      olr(:,:) = olr_local(:,:,1)                                       &
               + olr_local(:,:,2)
      lw_down(:,:) = lw_down_local(:,:,1)                               &
                   + lw_down_local(:,:,2)
      lwsea(:,:) = lwsea_local(:,:,1)                                   &
                 + lwsea_local(:,:,2)

      IF (l_extra_top) THEN
        top_absorption(:,:) = top_abs_lw(:,:,1)                         &
                            + top_abs_lw(:,:,2)
      END IF
    END IF

   ! Downward LW flux diagnostics weighted by sea ice fraction:
    IF (l_lwdn_sice_wt_cat .OR. l_lwdn_sice_wt ) THEN
      ALLOCATE ( lw_down_sice_weighted(row_length,rows) )
      lw_down_sice_weighted(:,:) = 0.0
      ALLOCATE ( lw_down_sice_weighted_cat(row_length,rows,nice_use) )
      lw_down_sice_weighted_cat(:,:,:) = 0.0
      DO n = 1, nice_use
        DO point = 1, sice_pts_ncat(n)
          l = sice_index_ncat(point,n)
          i = ssi_index_i(l)
          j = ssi_index_j(l)
          IF (l_lwdn_sice_wt_cat) lw_down_sice_weighted_cat(i,j,n) =    &
                                    lw_down(i,j) * ice_fract_cat(i,j,n)
          IF (l_lwdn_sice_wt)  lw_down_sice_weighted(i,j) =             &
                        lw_down_sice_weighted(i,j)                      &
                              + lw_down(i,j) * ice_fract_cat(i,j,n)
        END DO
      END DO
    END IF

!$OMP  PARALLEL DEFAULT(NONE)                                                  &
!$OMP& SHARED( l_ctile, rows, row_length, surflw, lw_incs, flandg,             &
!$OMP&         lwsea, land_sea_mask, dOLR_rts, olr, surf_emission )            &
!$OMP& PRIVATE( i, j )
    IF ( l_ctile ) THEN
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          surflw(i,j) = lw_incs(i,j,0)                                  &
                  + (1.0-flandg(i,j))*lwsea(i,j)
          IF (flandg(i,j) == 1.0) lwsea(i,j)=rmdi
        END DO
      END DO
!$OMP END DO NOWAIT
    ELSE
!$OMP DO SCHEDULE(STATIC)
      DO j = 1, rows
        DO i = 1, row_length
          surflw(i,j) = lw_incs(i,j,0) + lwsea(i,j)
          IF (land_sea_mask(i,j)) lwsea(i,j)=rmdi
        END DO
      END DO
!$OMP END DO NOWAIT
    END IF

    ! To allow for changes in surface temperature on non-radiation
    ! timesteps, dOLR is now defined as the grid-box mean OLR
    ! less the contributions to the upward LW flux at the surface
    ! from each portions of the surface, so conceptually
    ! emissivity * areal fraction * sigma T^4 for each portion.
    ! This is based on the effective temperatures seen by radiation,
    ! contributions haveing been accumulated in surf_emission.

!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length
        dOLR_rts(i,j) = olr(i,j) - sbcon * surf_emission(i,j)
      END DO
    END DO
!$OMP END DO
!$OMP END PARALLEL

  END IF ! end conditional on being a radiation timestep


  ! Is the PC2 cloud scheme being used?
  IF (i_cld_vn == i_cld_pc2) THEN

    ! Reset _latest values to _n values
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                              &
!$OMP& SHARED( model_levels, rows, row_length, t_latest, t_n, q_latest, q_n,   &
!$OMP&         qcl_latest, qcl_n, cf_latest, cf_n, cfl_latest, cfl_n )         &
!$OMP& PRIVATE( i, j, k )
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_latest(i,j,k)   = t_n(i,j,k)

          q_latest(i,j,k)   = q_n(i,j,k)
          qcl_latest(i,j,k) = qcl_n(i,j,k)
          cf_latest(i,j,k)  = cf_n(i,j,k)
          cfl_latest(i,j,k) = cfl_n(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

    ! ----------------------------------------------------------------------
    ! Homogeneous forcing. Note the temperature increment from longwave
    ! is added in this routine
    ! ----------------------------------------------------------------------

    CALL pc2_homog_plus_turb(p_layer_centres(1,1,1),                    &
      model_levels,                                                     &
      timestep, t_latest, cf_latest, cfl_latest,                        &
      cff_latest, q_latest, qcl_latest, lw_incs(1,1,1),                 &
      zeros, zeros, zeros, 0.0, 0.0,                                    &
      l_mr_physics)

    ! Add increments from the homogeneous forcing to the increment variables
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE ( i, j, k )
!$OMP DO SCHEDULE(STATIC)
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_inc(i,j,k) = t_inc(i,j,k) + t_latest(i,j,k)-t_n(i,j,k)

          q_inc(i,j,k) = q_inc(i,j,k) + q_latest(i,j,k)-q_n(i,j,k)
          qcl_inc(i,j,k) = qcl_inc(i,j,k)                               &
                           + qcl_latest(i,j,k)-qcl_n(i,j,k)
          cf_inc(i,j,k) = cf_inc(i,j,k)                                 &
                           + cf_latest(i,j,k)-cf_n(i,j,k)
          cfl_inc(i,j,k) = cfl_inc(i,j,k)                               &
                           + cfl_latest(i,j,k)-cfl_n(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT

   ! Compute and save LW Tendencies for physics_tendencies_mod
    IF (l_retain_rad_tendencies) THEN
!$OMP DO SCHEDULE(STATIC)
      DO k = 1, model_levels
        DO j = 1, rows
          DO i =  1, row_length
            dq_lw(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
          END DO
        END DO
      END DO
!$OMP END DO
    END IF !end if over l_retain_rad_tendencies
!$OMP END PARALLEL

  ELSE  ! i_cld_pc2

    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_inc(i,j,k) = t_inc(i,j,k) + lw_incs(i,j,k)
          t_latest(i,j,k) = t_n(i,j,k) + lw_incs(i,j,k)

        END DO
      END DO
    END DO

  END IF  ! i_cld_pc2

  ! Compute and save T LW Tendency for physics_tendencies_mod
  IF (l_retain_rad_tendencies) THEN
    DO k = 1, model_levels
      DO j = 1, rows
        DO i =  1, row_length
          dt_lw(i,j,k) = T_latest(i,j,k) - T_n(i,j,k)
        END DO
      END DO
    END DO
  END IF !end if over l_retain_rad_tendencies


  ! Get T_incr for output as STASH diagnostic
  IF ( l_t_incr_lw ) THEN  ! STASHflag set
    ALLOCATE ( t_incr_diagnostic(row_length,rows,model_levels) )
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          t_incr_diagnostic(i,j,k) =  lw_incs(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k


    IF ( l_scmdiags(scmdiag_pc2) .AND.                                  &
         model_type == mt_single_column) THEN

      ! Calculate equivalent to stash code 2,181 - the total
      ! temperature increment including the pc2 scheme
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            TmpScm3d_2(i,j,k) = t_latest(i,j,k) - t_n(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

      ! stash item = 2,181  ! temperature increment minus pc2
      CALL scmoutput(TmpScm3d_2,'dt_lwpc2',                             &
           'LW heating rate incl PC2','K/timestep',                     &
           t_acc,d_all,default_streams,'',routinename)

      CALL scmoutput(TmpScm3d_2,'lw2pc2',                               &
           'LW heating rate incl PC2','K/day',                          &
           t_mult,d_all,default_streams,'ntspday',routinename)

    END IF ! scmdiag_pc2 /model_type

  ELSE
    ALLOCATE ( t_incr_diagnostic(1,1,1) )
  END IF ! on STASHflag

  ! Set up radiative heating rates for 6A boundary layer code
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE( i, j, k )
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = 1, rows
      DO i = 1, row_length
        rad_hr(i,j,1,k) = lw_incs(i,j,k) * recip_timestep
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

  ! Calculate surface temperature from the radiative balance of the
  ! surface fluxes and a set intrinsic planet flux/temperature
  IF (l_planet_intrinsic_flux) THEN
!$OMP DO SCHEDULE(STATIC)
    DO j = 1, rows
      DO i = 1, row_length

        t_surf(i,j) = ( planet_t_intrinsic**4                           &
             + (lw_down(i,j)*planet_emissivity + surfsw_cor(i,j))       &
             / (planet_emissivity*sbcon)                                &
             )**0.25
      END DO
    END DO
!$OMP END DO NOWAIT

  END IF

  ! Copy dOLR from last LW timestep
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      dOLR(i,j) = dOLR_rts(i,j)
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  IF (ltimer) CALL timer ('AP1R LW Rad  ',6)

  ! ----------------------------------------------------------------------
  ! Section RAD.2.2 Long Wave Radiation Energy correction code
  ! -----------------------------------------------------------------------

  IF ((l_rad_step_prog .OR. l_rad_step_diag) .AND. l_emcorr) THEN

    ! Sum long wave fluxes into the atmosphere and
    ! add into the net diabatic fluxes into the
    ! atmosphere for use in the energy correction
    ! procedure

    IF (l_extra_top) THEN
      ! The energy absorbed above the top of the model in
      ! the radiation scheme does not contribute to the
      ! energy absorbed, but the diagnostics are calculated
      ! at the top of the atmosphere, so the net atmospheric
      ! flux must be adjusted.
      DO j = 1, rows
        DO i = 1, row_length
          net_atm_flux(i,j) = -olr(i,j) - surflw(i,j)                   &
                              -top_absorption(i,j)
        END DO
      END DO
    ELSE
      DO j = 1, rows
        DO i = 1, row_length
          net_atm_flux(i,j) = -olr(i,j) - surflw(i,j)
        END DO
      END DO
    END IF

    CALL flux_diag(net_atm_flux, cos_theta_latitude,                    &
      row_length, rows ,off_x,off_y,1.0,                                &
      sum_eng_fluxes,radiation_tstep_diag)

  END IF

  ! ----------------------------------------------------------------------
  ! Section RAD.2.3 Long Wave Radiation diagnostics
  ! -----------------------------------------------------------------------

  ! Minimal allocation for unused diagnostics
  IF (.NOT. ALLOCATED( lw_down_sice_weighted) ) &
            ALLOCATE ( lw_down_sice_weighted(1,1) )
  IF (.NOT. ALLOCATED( lw_down_sice_weighted_cat) ) &
            ALLOCATE ( lw_down_sice_weighted_cat(1,1,1) )

  DO j_lw = 1, n_lwcall

    ! Check that lw diagnostics requested this timestep

    SELECT CASE (model_type)
    CASE (mt_single_column)
      TmpLogic = ( error_code == 0 )
    CASE DEFAULT
      TmpLogic = ( ( error_code == 0 ) .AND. sf_calc(0,2) )
    END SELECT

    IF (TmpLogic) THEN

      IF (l_timestep) THEN

        i_off=0

        IF (j_lw == n_lwcall) THEN

          IF (l_rad_perturb .AND. l_rad_step_prog) THEN

            IF (lw_diag(2)%l_flux_up)                                   &
              lw_diag(1)%flux_up =                                      &
              lw_diag(1)%flux_up -                                      &
              lw_diag(2)%flux_up

            IF (lw_diag(2)%l_flux_down)                                 &
              lw_diag(1)%flux_down =                                    &
              lw_diag(1)%flux_down -                                    &
              lw_diag(2)%flux_down

            IF (lw_diag(2)%l_net_flux_trop)                             &
              lw_diag(1)%net_flux_trop =                                &
              lw_diag(1)%net_flux_trop -                                &
              lw_diag(2)%net_flux_trop

            IF (lw_diag(2)%l_down_flux_trop)                            &
              lw_diag(1)%down_flux_trop =                               &
              lw_diag(1)%down_flux_trop -                               &
              lw_diag(2)%down_flux_trop

          END IF

          IF (l_rad_step_diag) THEN

            IF (lw_diag(2)%l_flux_up)                                   &
              lw_diag(2)%flux_up =                                      &
              lw_diag(2)%flux_up +                                      &
              lw_diag(1)%flux_up

            IF (lw_diag(2)%l_flux_down)                                 &
              lw_diag(2)%flux_down =                                    &
              lw_diag(2)%flux_down +                                    &
              lw_diag(1)%flux_down

            IF (lw_diag(2)%l_net_flux_trop)                             &
              lw_diag(2)%net_flux_trop =                                &
              lw_diag(2)%net_flux_trop +                                &
              lw_diag(1)%net_flux_trop

            IF (lw_diag(2)%l_down_flux_trop)                            &
              lw_diag(2)%down_flux_trop =                               &
              lw_diag(2)%down_flux_trop +                               &
              lw_diag(1)%down_flux_trop

          END IF ! l_rad_step_diag
        END IF ! j_lw == n_lwcall

      ELSE IF (l_forcing .AND. l_rad_step_diag) THEN
        i_off=diagnostic_offset*(j_lw-1)
      ELSE IF (l_radiance) THEN
        IF (j_lw == 1) THEN
          i_off=0
        ELSE
          i_off=90+(diagnostic_offset/20)*(j_lw-1)
        END IF
      ELSE
        i_off=0
      END IF

      IF (model_type /= mt_single_column) THEN
        CALL diagnostics_rad(                                                  &
          lw_sect, j_lw, lw_diag(j_lw), lw_spectrum(j_lw),                     &
          row_length, rows, ozone_levels, cloud_levels, ntiles,                &
          at_extremity, i_off,                                                 &
          ! Fields common to SW and LW radiation
          t_n, t_inc, q_n, qcl_n, cf_n, cfl_n,                                 &
          t_latest, q_latest, qcl_latest, cf_latest, cfl_latest,               &
          t_incr_diagnostic, surflw, lwsea, obs_solid_angle,                   &
          lw_down_sice_weighted_cat, lw_down_sice_weighted,                    &
          ! SW diagnostic fields (dummys here)
          itoasw, surfsw_cor, toasw_cor, surfdir_cor, surfdif_cor,             &
          flux_below_690nm_surf, photosynth_act_rad, flxdirparsurf,            &
          dummy2d, sw_net_land, sw_net_sice, sea_salt_film, sea_salt_jet,      &
          dummy3d, dummy2d, dummy3d, dummy2d,                                  &
          salt_dim1, salt_dim2, salt_dim3,                                     &
          cos_zenith_angle, day_fraction,                                      &
          cos_zen_rts, day_frac_rts, sol_azm_rts,                              &
          ! LW diagnostic fields
          d_mass, density, layer_heat_capacity,                                &
          p_layer_boundaries, p_layer_centres, p_extra_layer,                  &
          t_layer_boundaries, t_extra_layer,                                   &
          olr, lw_down, ozone, o3_trop_level, o3_trop_height,                  &
          t_trop_level, t_trop_height,                                         &
          stashwork2)
      END IF ! model_type
    END IF ! Tmplogic
  END DO

  IF ( l_scmdiags(scmdiag_rad) .AND.                                    &
       model_type == mt_single_column) THEN

    ! Output some SCM diagnostics for LW radiation
    ! stash 2,161
    CALL scmoutput(t_incr_diagnostic,'dt_lw',                           &
         'LW heating rate minus PC2','K/timestep',                      &
         t_acc,d_all,default_streams,'',routinename)

    CALL scmoutput(t_incr_diagnostic,'lw2',                             &
         'LW heating rate minus PC2','K/day',                           &
         t_mult,d_all,default_streams,'ntspday',routinename)

    CALL scmoutput(surflw,'surf_lw',                                    &
         'Net surface LW flux','W/m2',                                  &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    ! stash 2,205
    CALL scmoutput(olr,'olr_toa',                                       &
         'Outgoing LW','W/m2',                                          &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    ! stash 2,182  Vapour increment
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_1(i,j,k) = q_latest(i,j,k) - q_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_1,'dq_lw',                                  &
         'Specific humidity increment lwrad','kg/kg',                   &
         t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    ! stash 2,183  liquid water content increment
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_1(i,j,k) = qcl_latest(i,j,k) - qcl_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_1,'dqcl_lw',                                &
         'qcl increment lwrad','kg/kg',                                 &
         t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    ! stash 2,192  total cloud fraction increment
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_1(i,j,k) = cf_latest(i,j,k) - cf_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_1,'dbcf_lw',                                &
         'bulk cloud fraction increment lwrad','fraction',              &
         t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    ! stash 2,193  liquid cloud fraction increment
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          TmpScm3d_1(i,j,k) = cfl_latest(i,j,k) - cfl_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k
    CALL scmoutput(TmpScm3d_1,'dcfl_lw',                                &
         'liquid cloud fraction increment lwrad','fraction',            &
         t_avg+only_radsteps,d_wet,default_streams,'',routinename)

    IF (l_timestep) j_lw=2
    IF (.NOT. l_timestep) j_lw=1

    ! stash 2,204
    IF ( ALLOCATED(lw_diag(j_lw)%total_cloud_cover) ) THEN
      TmpScm2d(:,:) = lw_diag(j_lw)%total_cloud_cover(:,:)
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm2d,'tca_lw',                                   &
         'Total cloud amount in LW rad','fraction',                     &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)

    CALL scmoutput(lw_down,'surf_dnlw',                                 &
         'Downward LW surface flux','W/m2',                             &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)


    ! In-radiation cloud diagnostics --------------------------------------
    IF ( ALLOCATED(lw_diag(j_lw)%ls_qcl_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%ls_qcl_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'ls_qcl_rad',                             &
         'Stratiform cloud liquid water','kg/kg',                       &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(lw_diag(j_lw)%ls_qcf_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%ls_qcf_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'ls_qcf_rad',                             &
         'Stratiform cloud ice water','kg/kg',                          &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(lw_diag(j_lw)%cc_qcl_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%cc_qcl_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'cc_qcl_rad',                             &
         'Convective cloud liquid water','kg/kg',                       &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(lw_diag(j_lw)%cc_qcf_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%cc_qcf_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'cc_qcf_rad',                             &
         'Convective cloud ice water','kg/kg',                          &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(lw_diag(j_lw)%ls_cl_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%ls_cl_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'ls_cl_rad',                              &
         'Stratiform liquid cloud fraction','',                         &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(lw_diag(j_lw)%ls_cf_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%ls_cf_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'ls_cf_rad',                              &
         'Stratiform ice cloud fraction','',                            &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(lw_diag(j_lw)%cc_cl_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%cc_cl_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'ccal_rad',                               &
         'Convective liquid cloud fraction','',                         &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF ( ALLOCATED(lw_diag(j_lw)%cc_cf_rad) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%cc_cf_rad(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'ccaf_rad',                               &
         'Convective ice cloud fraction','',                            &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

    IF (l_rad_perturb) j_lw=1

    ! stash 2,206
    IF ( ALLOCATED(lw_diag(j_lw)%clear_olr) ) THEN
      TmpScm2d(:,:) = lw_diag(j_lw)%clear_olr(:,:)
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm2d,                                            &
         'cs_olr','Clear-sky outgoing LW','W/m2',                       &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)
    ! stash 2,208
    IF ( ALLOCATED(lw_diag(j_lw)%surf_down_clr) ) THEN
      TmpScm2d(:,:) = lw_diag(j_lw)%surf_down_clr(:,:)
    ELSE
      TmpScm2d(:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm2d,'cs_surf_dnlw',                             &
         'Clear-sky down LW flux','W/m2',                               &
         t_avg+only_radsteps,d_sl,default_streams,'',routinename)
    ! stash 2,233
    IF ( ALLOCATED(lw_diag(j_lw)%clear_hr) ) THEN
      TmpScm3d_2(:,:,:) = lw_diag(j_lw)%clear_hr(:,:,:)
    ELSE
      TmpScm3d_2(:,:,:) = 0.0
    END IF
    CALL scmoutput(TmpScm3d_2,'dt_cslw',                                &
         'Clear-sky LW heating rates','K/s',                            &
         t_avg+only_radsteps,d_all,default_streams,'',routinename)

  END IF ! scmdiag_rad / model_type


  DEALLOCATE ( t_incr_diagnostic )

  DEALLOCATE ( lw_down_sice_weighted_cat )
  DEALLOCATE ( lw_down_sice_weighted     )

  ! Deallocate the diagnostic space that is no longer required.
  IF (l_timestep) THEN
    IF (l_rad_step_diag) THEN
      CALL deallocate_diag(lw_diag(2))
    END IF
    IF (timestep_number >  1) THEN
      IF (MOD(timestep_number,a_lw_radstep_prog) == 0) THEN

        DEALLOCATE(lw_incs_local)
        DEALLOCATE(lwsea_local)
        DEALLOCATE(olr_local)
        DEALLOCATE(lw_down_local)
        DEALLOCATE(top_abs_lw)

        CALL deallocate_diag(lw_diag(1))

      END IF
    END IF
  ELSE IF (l_forcing) THEN
    IF (l_rad_step_prog) THEN
      CALL deallocate_diag(lw_diag(1))
    END IF
    IF (l_rad_step_diag) THEN
      CALL deallocate_diag(lw_diag(2))
    END IF
  ELSE IF (l_radiance) THEN
    IF (l_rad_step_prog) THEN
      DO j_lw = 1, n_lwcall
        CALL deallocate_diag(lw_diag(j_lw))
      END DO
    END IF
  ELSE
    IF (l_rad_step_prog) THEN
      CALL deallocate_diag(lw_diag(1))
    END IF
  END IF

  ! Deallocate astronomy fields prior to next full radiation timestep
  IF (MOD(timestep_number,a_sw_radstep_prog) == 0) THEN
    DEALLOCATE(day_frac_sph_rts)
    DEALLOCATE(cos_zen_sph_rts)
    DEALLOCATE(day_frac_rts)
    DEALLOCATE(sol_azm_rts)
    DEALLOCATE(cos_zen_rts)
  END IF

END IF ! on error_code

! Restore the true tiled temperature in aggregate runs.
IF (l_aggregate) tstar_tile(:,1)=tstar_tile_tmp

9999 CONTINUE
! Check error condition
IF (error_code /= 0) THEN
  CALL ereport(routinename, error_code, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE rad_ctl
END MODULE rad_ctl_mod
