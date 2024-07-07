! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to declare a structure for radiation diagnostic fields.
!
! Description:
!   This module contains the declaration of the structure
!   used to store diagnostic fields for the radiation code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE def_diag

IMPLICIT NONE

TYPE StrDiag

! Have diagnostics been requested that require a second diagnostic call
  LOGICAL :: l_diag_call                        = .FALSE.

! SW and LW
  LOGICAL :: l_flux_up                          = .FALSE.
  LOGICAL :: l_flux_down                        = .FALSE.
  LOGICAL :: l_flux_up_clear                    = .FALSE.
  LOGICAL :: l_flux_down_clear                  = .FALSE.
  LOGICAL :: l_flux_up_clean                    = .FALSE.
  LOGICAL :: l_flux_down_clean                  = .FALSE.
  LOGICAL :: l_flux_up_clear_clean              = .FALSE.
  LOGICAL :: l_flux_down_clear_clean            = .FALSE.
  LOGICAL :: l_flux_up_band                     = .FALSE.
  LOGICAL :: l_flux_down_band                   = .FALSE.
  LOGICAL :: l_flux_up_clear_band               = .FALSE.
  LOGICAL :: l_flux_down_clear_band             = .FALSE.
  LOGICAL :: l_flux_up_clean_band               = .FALSE.
  LOGICAL :: l_flux_down_clean_band             = .FALSE.
  LOGICAL :: l_flux_up_clear_clean_band         = .FALSE.
  LOGICAL :: l_flux_down_clear_clean_band       = .FALSE.
  LOGICAL :: l_flux_up_forc                     = .FALSE.
  LOGICAL :: l_flux_down_forc                   = .FALSE.
  LOGICAL :: l_flux_up_clear_forc               = .FALSE.
  LOGICAL :: l_flux_down_clear_forc             = .FALSE.
  LOGICAL :: l_flux_up_forc_band                = .FALSE.
  LOGICAL :: l_flux_down_forc_band              = .FALSE.
  LOGICAL :: l_flux_up_clear_forc_band          = .FALSE.
  LOGICAL :: l_flux_down_clear_forc_band        = .FALSE.
  LOGICAL :: l_surf_down_clr                    = .FALSE.
  LOGICAL :: l_clear_hr                         = .FALSE.
  LOGICAL :: l_net_flux_trop                    = .FALSE.
  LOGICAL :: l_toa_radiance                     = .FALSE.
  LOGICAL :: l_aerosol_optical_depth            = .FALSE.
  LOGICAL :: l_aerosol_scat_optical_depth       = .FALSE.
  LOGICAL :: l_aerosol_asymmetry_scat           = .FALSE.
  LOGICAL :: l_rad_mask                         = .FALSE.
  LOGICAL :: l_emission_spectrum                = .FALSE.
  LOGICAL :: l_emission_spectrum_clear          = .FALSE.
  LOGICAL :: l_emission_spectrum_clean          = .FALSE.
  LOGICAL :: l_emission_spectrum_clear_clean    = .FALSE.
  LOGICAL :: l_toa_clear_weighted               = .FALSE.
  LOGICAL :: l_total_clear_area                 = .FALSE.
  LOGICAL :: l_easyaerosol_extinction           = .FALSE.
  LOGICAL :: l_easyaerosol_absorption           = .FALSE.
  LOGICAL :: l_easyaerosol_scattering           = .FALSE.
  LOGICAL :: l_easyaerosol_asytimscat           = .FALSE.

  REAL, ALLOCATABLE :: flux_up(:, :, :)
  !                      Upward fluxes on model levels
  REAL, ALLOCATABLE :: flux_down(:, :, :)
  !                      Downward fluxes on model levels
  REAL, ALLOCATABLE :: flux_up_clear(:, :, :)
  !                      Clear-sky upward fluxes on model levels
  REAL, ALLOCATABLE :: flux_down_clear(:, :, :)
  !                      Clear-sky downward fluxes on model levels
  REAL, ALLOCATABLE :: flux_up_clean(:, :, :)
  !                      Clean-air upward fluxes on model levels
  REAL, ALLOCATABLE :: flux_down_clean(:, :, :)
  !                      Clean-air downward fluxes on model levels
  REAL, ALLOCATABLE :: flux_up_clear_clean(:, :, :)
  !                      Clear-sky, clean-air upward fluxes on model levels
  REAL, ALLOCATABLE :: flux_down_clear_clean(:, :, :)
  !                      Clear-sky, clean-air downward fluxes on model levels
  REAL, ALLOCATABLE :: flux_up_band(:, :, :, :)
  !                      Upward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_down_band(:, :, :, :)
  !                      Downward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_up_clear_band(:, :, :, :)
  !                      Clear-sky upward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_down_clear_band(:, :, :, :)
  !                      Clear-sky downward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_up_clean_band(:, :, :, :)
  !                      Clean-air upward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_down_clean_band(:, :, :, :)
  !                      Clean-air downward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_up_clear_clean_band(:, :, :, :)
  !                      Clear, clean upward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_down_clear_clean_band(:, :, :, :)
  !                      Clear, clean downward fluxes on model levels and bands
  REAL, ALLOCATABLE :: flux_up_forc(:, :, :)
  !                      Upward fluxes on model levels with forcing
  REAL, ALLOCATABLE :: flux_down_forc(:, :, :)
  !                      Downward fluxes on model levels with forcing
  REAL, ALLOCATABLE :: flux_up_clear_forc(:, :, :)
  !                      Clear-sky upward fluxes on model levels with forcing
  REAL, ALLOCATABLE :: flux_down_clear_forc(:, :, :)
  !                      Clear-sky downward fluxes on model levels with forcing
  REAL, ALLOCATABLE :: flux_up_forc_band(:, :, :, :)
  !                      Upward fluxes on model levels and bands with forcing
  REAL, ALLOCATABLE :: flux_down_forc_band(:, :, :, :)
  !                      Downward fluxes on model levels and bands with forcing
  REAL, ALLOCATABLE :: flux_up_clear_forc_band(:, :, :, :)
  !                      Clear-sky upward fluxes on model levels and bands
  !                      with forcing
  REAL, ALLOCATABLE :: flux_down_clear_forc_band(:, :, :, :)
  !                      Clear-sky downward fluxes on model levels and bands
  !                      with forcing
  REAL, ALLOCATABLE :: surf_down_clr(:, :)
  !                      Clear-sky downward flux at the
  !                      surface (Method II)
  REAL, ALLOCATABLE :: clear_hr(:, :, :)
  !                      Clear-sky heating rates, calculated
  !                      by ignoring all clouds (Method II):
  !                      these are not necessarily the same
  !                      as the heating rates in the cloud-
  !                      free parts of a grid-box
  REAL, ALLOCATABLE :: net_flux_trop(:, :)
  !                      Net downward flux at the tropopause
  REAL, ALLOCATABLE :: toa_radiance(:, :, :)
  !                      The radiance observed at the top of
  !                      the atmosphere
  REAL, ALLOCATABLE :: aerosol_optical_depth(:, :, :, :)
  !                      Total aerosol optical depth per model level and band
  REAL, ALLOCATABLE :: aerosol_scat_optical_depth(:, :, :, :)
  !                      Total aerosol scattering optical depth
  !                      per model level and band
  REAL, ALLOCATABLE :: aerosol_asymmetry_scat(:, :, :, :)
  !                      Total aerosol asymmetry per model level and band
  !                      weighted by scattering optical depth
  REAL, ALLOCATABLE :: rad_mask(:, :)
  !                      Mask for radiation calculations:
  !                        1.0 if fluxes have been calculated
  !                        0.0 if no calculations have been done
  REAL, ALLOCATABLE :: emission_spectrum(:, :, :)
  !                      Outgoing TOA flux towards a distant observer
  !                      in each band normalised to units of Wm-2 at 1 AU
  REAL, ALLOCATABLE :: emission_spectrum_clear(:, :, :)
  !                      Outgoing TOA clear-sky flux towards a distant observer
  !                      in each band normalised to units of Wm-2 at 1 AU
  REAL, ALLOCATABLE :: emission_spectrum_clean(:, :, :)
  !                      Outgoing TOA clean-air flux towards a distant observer
  !                      in each band normalised to units of Wm-2 at 1 AU
  REAL, ALLOCATABLE :: emission_spectrum_clear_clean(:, :, :)
  !                      Outgoing TOA clear-sky, clean-air flux towards a
  !                      distant observer in each band (Wm-2 at 1 AU)
  REAL, ALLOCATABLE :: toa_clear_weighted(:, :)
  !                      Clear-sky (method I) upward flux at top-of-atmosphere
  !                      weighted by the clear-sky fraction (total_clear_area)
  REAL, ALLOCATABLE :: total_clear_area(:, :)
  !                      Total clear-sky fraction in the column
  !                      (1 - total_cloud_cover)
  REAL, ALLOCATABLE :: easyaerosol_extinction(:, :, :, :)
  !                      Extinction of EasyAerosol per model level and
  !                      waveband.
  REAL, ALLOCATABLE :: easyaerosol_absorption(:, :, :, :)
  !                      Absorption of EasyAerosol per model level and
  !                      waveband.
  REAL, ALLOCATABLE :: easyaerosol_scattering(:, :, :, :)
  !                      Scattering of EasyAerosol per model level and
  !                      waveband.
  REAL, ALLOCATABLE :: easyaerosol_asytimscat(:, :, :, :)
  !                      Scattering-weighted asymmetry of EasyAerosol
  !                      per model level and waveband.

! SW
  LOGICAL :: l_flux_direct_sph                  = .FALSE.
  LOGICAL :: l_flux_direct_div                  = .FALSE.
  LOGICAL :: l_flux_direct_clear_sph            = .FALSE.
  LOGICAL :: l_flux_direct_clear_div            = .FALSE.
  LOGICAL :: l_flux_direct_clean_sph            = .FALSE.
  LOGICAL :: l_flux_direct_clean_div            = .FALSE.
  LOGICAL :: l_flux_direct_clear_clean_sph      = .FALSE.
  LOGICAL :: l_flux_direct_clear_clean_div      = .FALSE.
  LOGICAL :: l_flux_direct_sph_band             = .FALSE.
  LOGICAL :: l_flux_direct_div_band             = .FALSE.
  LOGICAL :: l_flux_direct_clear_sph_band       = .FALSE.
  LOGICAL :: l_flux_direct_clear_div_band       = .FALSE.
  LOGICAL :: l_flux_direct_clean_sph_band       = .FALSE.
  LOGICAL :: l_flux_direct_clean_div_band       = .FALSE.
  LOGICAL :: l_flux_direct_clear_clean_sph_band = .FALSE.
  LOGICAL :: l_flux_direct_clear_clean_div_band = .FALSE.
  LOGICAL :: l_transmission_spectrum            = .FALSE.
  LOGICAL :: l_transmission_spectrum_clear      = .FALSE.
  LOGICAL :: l_transmission_spectrum_clean      = .FALSE.
  LOGICAL :: l_transmission_spectrum_clear_clean= .FALSE.
  LOGICAL :: l_flux_direct_sph_forc             = .FALSE.
  LOGICAL :: l_flux_direct_div_forc             = .FALSE.
  LOGICAL :: l_flux_direct_clear_sph_forc       = .FALSE.
  LOGICAL :: l_flux_direct_clear_div_forc       = .FALSE.
  LOGICAL :: l_flux_direct_sph_forc_band        = .FALSE.
  LOGICAL :: l_flux_direct_div_forc_band        = .FALSE.
  LOGICAL :: l_flux_direct_clear_sph_forc_band  = .FALSE.
  LOGICAL :: l_flux_direct_clear_div_forc_band  = .FALSE.
  LOGICAL :: l_spherical_path                   = .FALSE.
  LOGICAL :: l_solar_out_toa                    = .FALSE.
  LOGICAL :: l_solar_out_clear                  = .FALSE.
  LOGICAL :: l_surface_down_flux                = .FALSE.
  LOGICAL :: l_surf_up_clr                      = .FALSE.
  LOGICAL :: l_up_flux_trop                     = .FALSE.
  LOGICAL :: l_flux_direct                      = .FALSE.
  LOGICAL :: l_flux_diffuse                     = .FALSE.
  LOGICAL :: re_conv_flag                       = .FALSE.
  LOGICAL :: re_strat_flag                      = .FALSE.
  LOGICAL :: wgt_conv_flag                      = .FALSE.
  LOGICAL :: wgt_strat_flag                     = .FALSE.
  LOGICAL :: lwp_strat_flag                     = .FALSE.
  LOGICAL :: weighted_re_flag                   = .FALSE.
  LOGICAL :: sum_weight_re_flag                 = .FALSE.
  LOGICAL :: wgtd_warm_re_flag                  = .FALSE.
  LOGICAL :: sum_wgt_warm_re_flag               = .FALSE.
  LOGICAL :: cdnc_ct_diag_flag                  = .FALSE.
  LOGICAL :: cdnc_ct_weight_flag                = .FALSE. 
  LOGICAL :: ntot_diag_flag                     = .FALSE.
  LOGICAL :: strat_lwc_diag_flag                = .FALSE.
  LOGICAL :: so4_ccn_diag_flag                  = .FALSE.
  LOGICAL :: cond_samp_wgt_flag                 = .FALSE.
  LOGICAL :: seasalt_film_flag                  = .FALSE.
  LOGICAL :: seasalt_jet_flag                   = .FALSE.
  LOGICAL :: nc_diag_flag                       = .FALSE.
  LOGICAL :: nc_weight_flag                     = .FALSE.
  LOGICAL :: l_FlxSolBelow690nmSurf             = .FALSE.
  LOGICAL :: l_FlxSeaBelow690nmSurf             = .FALSE.
  LOGICAL :: l_cloud_extinction                 = .FALSE.
  LOGICAL :: l_cloud_weight_extinction          = .FALSE.
  LOGICAL :: l_ls_cloud_extinction              = .FALSE.
  LOGICAL :: l_ls_cloud_weight_extinction       = .FALSE.
  LOGICAL :: l_cnv_cloud_extinction             = .FALSE.
  LOGICAL :: l_cnv_cloud_weight_extinction      = .FALSE.
  LOGICAL :: l_orog_corr                        = .FALSE.
  LOGICAL :: l_uvflux_direct                    = .FALSE.
  LOGICAL :: l_uvflux_up                        = .FALSE.
  LOGICAL :: l_uvflux_down                      = .FALSE.
  LOGICAL :: l_surf_uv                          = .FALSE.
  LOGICAL :: l_surf_uv_clr                      = .FALSE.
  LOGICAL :: l_direct_albedo                    = .FALSE.
  LOGICAL :: l_diffuse_albedo                   = .FALSE.
  LOGICAL :: l_vis_albedo_sc                    = .FALSE.
  LOGICAL :: l_nir_albedo_sc                    = .FALSE.

  REAL, ALLOCATABLE :: flux_direct_sph(:, :, :)
  !                      Direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_div(:, :, :)
  !                      Direct flux divergence across layer
  REAL, ALLOCATABLE :: flux_direct_clear_sph(:, :, :)
  !                      Clear-sky direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_clear_div(:, :, :)
  !                      Clear-sky direct flux divergence across layer
  REAL, ALLOCATABLE :: flux_direct_clean_sph(:, :, :)
  !                      Clean-air direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_clean_div(:, :, :)
  !                      Clean-air direct flux divergence across layer
  REAL, ALLOCATABLE :: flux_direct_clear_clean_sph(:, :, :)
  !                      Clear, clean direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_clear_clean_div(:, :, :)
  !                      Clear, clean direct flux divergence across layer
  REAL, ALLOCATABLE :: flux_direct_sph_band(:, :, :, :)
  !                      Direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_div_band(:, :, :, :)
  !                      Direct flux divergence across layer
  REAL, ALLOCATABLE :: flux_direct_clear_sph_band(:, :, :, :)
  !                      Clear-sky direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_clear_div_band(:, :, :, :)
  !                      Clear-sky direct flux divergence across layer
  REAL, ALLOCATABLE :: flux_direct_clean_sph_band(:, :, :, :)
  !                      Clean-air direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_clean_div_band(:, :, :, :)
  !                      Clean-air direct flux divergence across layer
  REAL, ALLOCATABLE :: flux_direct_clear_clean_sph_band(:, :, :, :)
  !                      Clear, clean direct flux for spherical geometry
  REAL, ALLOCATABLE :: flux_direct_clear_clean_div_band(:, :, :, :)
  !                      Clear, clean direct flux divergence across layer
  REAL, ALLOCATABLE :: transmission_spectrum(:, :, :)
  !                      Outgoing TOA direct flux towards a distant observer
  !                      in each band normalised to units of Wm-2 at 1 AU
  REAL, ALLOCATABLE :: transmission_spectrum_clear(:, :, :)
  !                      Outgoing TOA clear-sky direct flux towards a distant
  !                      observer in each band (Wm-2 at 1 AU)
  REAL, ALLOCATABLE :: transmission_spectrum_clean(:, :, :)
  !                      Outgoing TOA clean-air direct flux towards a distant
  !                      observer in each band (Wm-2 at 1 AU)
  REAL, ALLOCATABLE :: transmission_spectrum_clear_clean(:, :, :)
  !                      Outgoing TOA clear-sky, clean-air direct flux towards
  !                      a distant observer in each band (Wm-2 at 1 AU)
  REAL, ALLOCATABLE :: flux_direct_sph_forc(:, :, :)
  !                      Direct flux for spherical geometry with forcing
  REAL, ALLOCATABLE :: flux_direct_div_forc(:, :, :)
  !                      Direct flux divergence across layer with forcing
  REAL, ALLOCATABLE :: flux_direct_clear_sph_forc(:, :, :)
  !                      Clear-sky direct flux for spherical geometry
  !                      with forcing
  REAL, ALLOCATABLE :: flux_direct_clear_div_forc(:, :, :)
  !                      Clear-sky direct flux divergence across layer
  !                      with forcing
  REAL, ALLOCATABLE :: flux_direct_sph_forc_band(:, :, :, :)
  !                      Direct flux for spherical geometry with forcing
  REAL, ALLOCATABLE :: flux_direct_div_forc_band(:, :, :, :)
  !                      Direct flux divergence across layer with forcing
  REAL, ALLOCATABLE :: flux_direct_clear_sph_forc_band(:, :, :, :)
  !                      Clear-sky direct flux for spherical geometry
  !                      with forcing
  REAL, ALLOCATABLE :: flux_direct_clear_div_forc_band(:, :, :, :)
  !                      Clear-sky direct flux divergence across layer
  !                      with forcing
  REAL, ALLOCATABLE :: spherical_path(:, :, :, :)
  !                      Path length for direct beam through spherical layers
  !                      as a multiple of vertical height
  REAL, ALLOCATABLE :: solar_out_toa(:, :)
  !                      Reflected SW flux at the top
  !                      of the atmosphere
  REAL, ALLOCATABLE :: solar_out_clear(:, :)
  !                      Reflected clear-sky SW flux
  !                      at the top of the atmosphere:
  !                      this is calculated at all points
  !                      omitting cloud (Method II)
  REAL, ALLOCATABLE :: surface_down_flux(:, :)
  !                      Downward SW flux at the surface
  !                      (not net)
  REAL, ALLOCATABLE :: surf_up_clr(:, :)
  !                      Clear-sky upward SW flux at the
  !                      surface (Method II)
  REAL, ALLOCATABLE :: up_flux_trop(:, :)
  !                      Actual upward flux at the tropopause
  REAL, ALLOCATABLE :: flux_direct(:, :, :)
  !                      Direct downward SW flux
  REAL, ALLOCATABLE :: flux_diffuse(:, :, :)
  !                      Diffuse downward SW flux
  REAL, ALLOCATABLE :: re_strat(:, :, :)
  !                      The weighted effective radius in
  !                      stratiform clouds multiplied by 10^6
  !                      to avoid packing problems
  REAL, ALLOCATABLE :: wgt_strat(:, :, :)
  !                      The weighting factor for re_strat
  !                      and lwp_strat
  REAL, ALLOCATABLE :: lwp_strat(:, :, :)
  !                      The liquid water path in stratiform
  !                      cloud weighted by wgt_strat
  REAL, ALLOCATABLE :: re_conv(:, :, :)
  !                      The weighted effective radius in
  !                      Convective clouds multiplied by 10^6
  !                      to avoid packing problems
  REAL, ALLOCATABLE :: wgt_conv(:, :, :)
  !                      The weighting factor for re_conv
  REAL, ALLOCATABLE :: ntot_diag(:, :, :)
  !                      The number concentration of droplets
  !                      multiplied by stratiform weighting factor
  REAL, ALLOCATABLE :: strat_lwc_diag(:, :, :)
  !                      The liquid water content of stratiform
  !                      clouds multiplied by stratiform weighting
  !                      factor
  REAL, ALLOCATABLE :: so4_ccn_diag(:, :, :)
  !                      The mass concentration of SO4 CCN multiplied
  !                      by the conditional sampling weight
  REAL, ALLOCATABLE :: cond_samp_wgt(:, :, :)
  !                      The conditional sampling weight for
  !                      so4_ccn_diag
  REAL, ALLOCATABLE :: weighted_re(:, :)
  !                      The effective radius as seen from space
  !                      multiplied by an appropriate weight
  REAL, ALLOCATABLE :: sum_weight_re(:, :)
  !                      The weighting factor for the effective
  !                      radius as viewed from space
  REAL, ALLOCATABLE :: weighted_warm_re(:, :)
  !                      The effective radius as seen from space
  !                      for warm clouds only (T>273K) multiplied
  !                      by an appropriate weight
  REAL, ALLOCATABLE :: sum_weight_warm_re(:, :)
  !                      The weighting factor for the warm-cloud-
  !                      only effective radius as viewed from space
  REAL, ALLOCATABLE :: cdnc_ct_diag(:, :)
  !                      CDNC at cloud top multiplied by
  !                      an appropriate weight
  REAL, ALLOCATABLE :: cdnc_ct_weight(:, :)
  !                      The weighting factor for the CDNC
  !                      at cloud top
  REAL, ALLOCATABLE :: nc_diag(:, :)
  !                      The column-integrated cloud droplet number
  !                      multiplied by an appropriate weight
  REAL, ALLOCATABLE :: nc_weight(:, :)
  !                      The weighting factor for the column-
  !                      integrated cloud droplet number
  REAL, ALLOCATABLE :: FlxSolBelow690nmSurf(:, :)
  !                      The grid-box mean flux below 690 nm
  !                      into the solid surface
  REAL, ALLOCATABLE :: FlxSeaBelow690nmSurf(:, :)
  !                      The grid-box mean flux below 690 nm
  !                      into the sea surface
  REAL, ALLOCATABLE :: cloud_extinction(:, :, :)
  !                      Mean extinction coefficient in clouds,
  !                      weighted by the cloud amount and the
  !                      clear sky flux
  REAL, ALLOCATABLE :: cloud_weight_extinction(:, :, :)
  !                      Weighting factor for extinction in clouds :
  !                      the product of the cloud amount and the
  !                      clear-sky direct flux.
  REAL, ALLOCATABLE :: ls_cloud_extinction(:, :, :)
  !                      Mean extinction coefficient in layer clouds,
  !                      weighted by the cloud amount and the
  !                      clear sky flux
  REAL, ALLOCATABLE :: ls_cloud_weight_extinction(:, :, :)
  !                      Weighting factor for extinction in layer
  !                      clouds : the product of the cloud amount
  !                      and the clear-sky direct flux.
  REAL, ALLOCATABLE :: cnv_cloud_extinction(:, :, :)
  !                      Mean extinction coefficient in conv. clouds,
  !                      weighted by the cloud amount and the
  !                      clear sky flux
  REAL, ALLOCATABLE :: cnv_cloud_weight_extinction(:, :, :)
  !                      Weighting factor for extinction in conv.
  !                      clouds : the product of the cloud amount
  !                      and the clear-sky direct flux.
  REAL, ALLOCATABLE :: orog_corr(:, :)
  !                      Correction factor for the direct solar flux
  !                      reaching the surface for sloping terrain.
  REAL, ALLOCATABLE :: uvflux_direct(:, :, :)
  !                      direct UV-flux
  REAL, ALLOCATABLE :: uvflux_down(:, :, :)
  !                      downward UV-flux
  REAL, ALLOCATABLE :: uvflux_up(:, :, :)
  !                      upward UV-flux
  REAL, ALLOCATABLE :: surf_uv(:, :)
  !                      Surface down UV flux
  REAL, ALLOCATABLE :: surf_uv_clr(:, :)
  !                      Clear-sky surface down UV flux
  REAL, ALLOCATABLE :: direct_albedo(:, :, :)
  !                      direct albedo on SW bands
  REAL, ALLOCATABLE :: diffuse_albedo(:, :, :)
  !                      diffuse albedo on SW bands
  REAL, ALLOCATABLE :: vis_albedo_sc(:, :, :)
  !                      sclaing to VIS albedo obs on land tiles
  REAL, ALLOCATABLE :: nir_albedo_sc(:, :, :)
  !                      sclaing to NIR albedo obs on land tiles


! LW
  LOGICAL :: l_total_cloud_cover                = .FALSE.
  LOGICAL :: l_clear_olr                        = .FALSE.
  LOGICAL :: l_down_flux_trop                   = .FALSE.
  LOGICAL :: l_total_cloud_on_levels            = .FALSE.
  LOGICAL :: l_cloud_absorptivity               = .FALSE.
  LOGICAL :: l_cloud_weight_absorptivity        = .FALSE.
  LOGICAL :: l_ls_cloud_absorptivity            = .FALSE.
  LOGICAL :: l_ls_cloud_weight_absorptivity     = .FALSE.
  LOGICAL :: l_cnv_cloud_absorptivity           = .FALSE.
  LOGICAL :: l_cnv_cloud_weight_absorptivity    = .FALSE.
  LOGICAL :: l_ls_qcl_rad                       = .FALSE.
  LOGICAL :: l_ls_qcf_rad                       = .FALSE.
  LOGICAL :: l_cc_qcl_rad                       = .FALSE.
  LOGICAL :: l_cc_qcf_rad                       = .FALSE.
  LOGICAL :: l_ls_cl_rad                        = .FALSE.
  LOGICAL :: l_ls_cf_rad                        = .FALSE.
  LOGICAL :: l_cc_cl_rad                        = .FALSE.
  LOGICAL :: l_cc_cf_rad                        = .FALSE.
  LOGICAL :: l_ccore_clt_rad                    = .FALSE.
  LOGICAL :: l_ccore_qcl_rad                    = .FALSE.
  LOGICAL :: l_ccore_qcf_rad                    = .FALSE.
  LOGICAL :: l_ls_del_rad                       = .FALSE.
  LOGICAL :: l_ls_def_rad                       = .FALSE.
  LOGICAL :: l_cc_del_rad                       = .FALSE.
  LOGICAL :: l_cc_def_rad                       = .FALSE.
  LOGICAL :: l_aod_sulphate                     = .FALSE.
  LOGICAL :: l_aod_dust                         = .FALSE.
  LOGICAL :: l_aod_seasalt                      = .FALSE.
  LOGICAL :: l_aod_soot                         = .FALSE.
  LOGICAL :: l_aod_biomass                      = .FALSE.
  LOGICAL :: l_aod_biogenic                     = .FALSE.
  LOGICAL :: l_aod_ocff                         = .FALSE.
  LOGICAL :: l_aod_delta                        = .FALSE.
  LOGICAL :: l_aod_nitrate                      = .FALSE.
  LOGICAL :: l_aod_total_radn                   = .FALSE.
  LOGICAL :: l_angst_total_radn                 = .FALSE.
  LOGICAL :: l_aod_prog_sulphate                = .FALSE.
  LOGICAL :: l_aod_prog_dust                    = .FALSE.
  LOGICAL :: l_aod_prog_seasalt                 = .FALSE.
  LOGICAL :: l_aod_prog_soot                    = .FALSE.
  LOGICAL :: l_aod_prog_biomass                 = .FALSE.
  LOGICAL :: l_aod_prog_ocff                    = .FALSE.
  LOGICAL :: l_aod_prog_nitrate                 = .FALSE.
  LOGICAL :: l_aaod_sulphate                    = .FALSE.
  LOGICAL :: l_aaod_dust                        = .FALSE.
  LOGICAL :: l_aaod_seasalt                     = .FALSE.
  LOGICAL :: l_aaod_soot                        = .FALSE.
  LOGICAL :: l_aaod_biomass                     = .FALSE.
  LOGICAL :: l_aaod_biogenic                    = .FALSE.
  LOGICAL :: l_aaod_ocff                        = .FALSE.
  LOGICAL :: l_aaod_nitrate                     = .FALSE.
  LOGICAL :: l_aod_ukca_ait_sol                 = .FALSE.
  LOGICAL :: l_aod_ukca_acc_sol                 = .FALSE.
  LOGICAL :: l_aod_ukca_cor_sol                 = .FALSE.
  LOGICAL :: l_aod_ukca_ait_ins                 = .FALSE.
  LOGICAL :: l_aod_ukca_acc_ins                 = .FALSE.
  LOGICAL :: l_aod_ukca_cor_ins                 = .FALSE.
  LOGICAL :: l_sod_ukca_ait_sol                 = .FALSE.
  LOGICAL :: l_sod_ukca_acc_sol                 = .FALSE.
  LOGICAL :: l_sod_ukca_cor_sol                 = .FALSE.
  LOGICAL :: l_sod_ukca_ait_ins                 = .FALSE.
  LOGICAL :: l_sod_ukca_acc_ins                 = .FALSE.
  LOGICAL :: l_sod_ukca_cor_ins                 = .FALSE.
  LOGICAL :: l_ls_qcl_rad_path                  = .FALSE.
  LOGICAL :: l_ls_qcf_rad_path                  = .FALSE.
  LOGICAL :: l_cc_qcl_rad_path                  = .FALSE.
  LOGICAL :: l_cc_qcf_rad_path                  = .FALSE.
  LOGICAL :: l_ccore_qcl_rad_path               = .FALSE.
  LOGICAL :: l_ccore_qcf_rad_path               = .FALSE.
  LOGICAL :: l_aaod_ukca_ait_sol                = .FALSE.
  LOGICAL :: l_aaod_ukca_acc_sol                = .FALSE.
  LOGICAL :: l_aaod_ukca_cor_sol                = .FALSE.
  LOGICAL :: l_aaod_ukca_ait_ins                = .FALSE.
  LOGICAL :: l_aaod_ukca_acc_ins                = .FALSE.
  LOGICAL :: l_aaod_ukca_cor_ins                = .FALSE.
  LOGICAL, ALLOCATABLE :: l_clas_aerosol_ext(:)
  LOGICAL, ALLOCATABLE :: l_clas_aerosol_abs(:)
  LOGICAL, ALLOCATABLE :: l_clas_aerosol_sca(:)
  LOGICAL, ALLOCATABLE :: l_ukca_aerosol_ext(:)
  LOGICAL, ALLOCATABLE :: l_ukca_aerosol_abs(:)
  LOGICAL, ALLOCATABLE :: l_ukca_aerosol_sca(:)
  LOGICAL, ALLOCATABLE :: l_ukca_aerosol_gsca(:)

  INTEGER :: n_clas_aerosol_ext
  INTEGER :: n_clas_aerosol_abs
  INTEGER :: n_clas_aerosol_sca
  INTEGER :: n_ukca_aerosol_ext
  INTEGER :: n_ukca_aerosol_abs
  INTEGER :: n_ukca_aerosol_sca
  INTEGER :: n_ukca_aerosol_gsca


  REAL, ALLOCATABLE :: total_cloud_cover(:, :)
  !                      Total cloud cover at all grid-points
  REAL, ALLOCATABLE :: clear_olr(:, :)
  !                      Clear-sky outgoing LW radiation calculated
  !                      at all grid-points omitting cloud
  !                      (This is known as Method II)
  REAL, ALLOCATABLE :: down_flux_trop(:, :)
  !                      The actual downward flux at the tropopause
  REAL, ALLOCATABLE :: total_cloud_on_levels(:, :, :)
  !                      Total cloud on model layers
  REAL, ALLOCATABLE :: cloud_absorptivity(:, :, :)
  !                      Mean absorption coefficient in clouds,
  !                      weighted by the cloud amount and the
  !                      clear sky flux
  REAL, ALLOCATABLE :: cloud_weight_absorptivity(:, :, :)
  !                      Weighting factor for absorption in clouds :
  !                      the product of the cloud amount and the
  !                      clear-sky direct flux.
  REAL, ALLOCATABLE :: ls_cloud_absorptivity(:, :, :)
  !                      Mean absorption coefficient in layer clouds,
  !                      weighted by the cloud amount and the
  !                      clear sky flux
  REAL, ALLOCATABLE :: ls_cloud_weight_absorptivity(:, :, :)
  !                      Weighting factor for absorption in layer
  !                      clouds : the product of the cloud amount
  !                      and the clear-sky direct flux.
  REAL, ALLOCATABLE :: cnv_cloud_absorptivity(:, :, :)
  !                      Mean absorption coefficient in conv. clouds,
  !                      weighted by the cloud amount and the
  !                      clear sky flux
  REAL, ALLOCATABLE :: cnv_cloud_weight_absorptivity(:, :, :)
  !                      Weighting factor for absorption in conv.
  !                      clouds : the product of the cloud amount
  !                      and the clear-sky direct flux.
  REAL, ALLOCATABLE :: ls_qcl_rad(:, :, :)
  !                      Grid-box mean LS liquid water mixing ratio
  REAL, ALLOCATABLE :: ls_qcf_rad(:, :, :)
  !                      Grid-box mean LS ice water mixing ratio
  REAL, ALLOCATABLE :: cc_qcl_rad(:, :, :)
  !                      Grid-box mean CONV liquid water mixing ratio
  REAL, ALLOCATABLE :: cc_qcf_rad(:, :, :)
  !                      Grid-box mean CONV ice water mixing ratio
  REAL, ALLOCATABLE :: ls_cl_rad(:, :, :)
  !                      Grid-box mean LS liquid cloud fraction
  REAL, ALLOCATABLE :: ls_cf_rad(:, :, :)
  !                      Grid-box mean LS ice cloud fraction
  REAL, ALLOCATABLE :: cc_cl_rad(:, :, :)
  !                      Grid-box mean CONV liquid cloud fraction
  REAL, ALLOCATABLE :: cc_cf_rad(:, :, :)
  !                      Grid-box mean CONV ice cloud fraction
  REAL, ALLOCATABLE :: ccore_clt_rad(:, :, :)
  !                      Convective core cloud fraction
  REAL, ALLOCATABLE :: ccore_qcl_rad(:, :, :)
  !                      Grid-box mean convective core liquid water mixing ratio
  REAL, ALLOCATABLE :: ccore_qcf_rad(:, :, :)
  !                      Grid-box mean convective core ice water mixing ratio
  REAL, ALLOCATABLE :: ls_del_rad(:, :, :)
  !                      Large-scale liquid weighted effective dimension
  REAL, ALLOCATABLE :: ls_def_rad(:, :, :)
  !                      Large-scale ice weighted effective dimension
  REAL, ALLOCATABLE :: cc_del_rad(:, :, :)
  !                      Convective liquid weighted effective dimension
  REAL, ALLOCATABLE :: cc_def_rad(:, :, :)
  !                      Convective ice weighted effective dimension
  REAL, ALLOCATABLE :: aod_sulphate(:, :, :)
  !                      Sulphate aerosol optical depth
  REAL, ALLOCATABLE :: aod_dust(:, :, :)
  !                      Mineral dust aerosol optical depth
  REAL, ALLOCATABLE :: aod_seasalt(:, :, :)
  !                      Sea salt aerosol optical depth
  REAL, ALLOCATABLE :: aod_soot(:, :, :)
  !                      Soot (black-carbon) aerosol optical depth
  REAL, ALLOCATABLE :: aod_biomass(:, :, :)
  !                      Biomass-burning aerosol optical depth
  REAL, ALLOCATABLE :: aod_biogenic(:, :, :)
  !                      Biogenic aerosol optical depth
  REAL, ALLOCATABLE :: aod_ocff(:, :, :)
  !                      Fossil-fuel org carb aerosol optical depth
  REAL, ALLOCATABLE :: aod_delta(:, :, :)
  !                      Delta aerosol optical depth
  REAL, ALLOCATABLE :: aod_nitrate(:, :, :)
  !                      Nitrate aerosol optical depth
  REAL, ALLOCATABLE :: aod_total_radn(:, :, :)
  !                      Total aerosol optical depth in radiation
  REAL, ALLOCATABLE :: angst_total_radn(:, :, :)
  !                      Angstom Exp from Total AOD in radiation
  REAL, ALLOCATABLE :: aod_prog_sulphate(:, :, :)
  !                      Prognostic Sulphate aerosol optical depth (radn or not)
  REAL, ALLOCATABLE :: aod_prog_dust(:, :, :)
  !                      Prognostic Mineral dust aerosol optical depth
  REAL, ALLOCATABLE :: aod_prog_seasalt(:, :, :)
  !                      Seasalt AOD diagnosed from prognostics
  REAL, ALLOCATABLE :: aod_prog_soot(:, :, :)
  !                      Prognostic Soot (black-carbon) aerosol optical depth
  REAL, ALLOCATABLE :: aod_prog_biomass(:, :, :)
  !                      Prognostic Biomass-burning aerosol optical depth
  REAL, ALLOCATABLE :: aod_prog_ocff(:, :, :)
  !                      Prognostic Fossil-fuel org carb aerosol optical depth
  REAL, ALLOCATABLE :: aod_prog_nitrate(:, :, :)
  !                      Prognostic Nitrate aerosol optical depth
  REAL, ALLOCATABLE :: aaod_sulphate(:, :, :)
  !                      Sulphate absorption aerosol optical depth
  REAL, ALLOCATABLE :: aaod_dust(:, :, :)
  !                      Mineral dust absorption aerosol optical depth
  REAL, ALLOCATABLE :: aaod_seasalt(:, :, :)
  !                      Sea salt absorption aerosol optical depth
  REAL, ALLOCATABLE :: aaod_soot(:, :, :)
  !                      Soot (black-carbon) absorption aerosol optical depth
  REAL, ALLOCATABLE :: aaod_biomass(:, :, :)
  !                      Biomass-burning absorption aerosol optical depth
  REAL, ALLOCATABLE :: aaod_biogenic(:, :, :)
  !                      Biogenic aerosol absorption optical depth
  REAL, ALLOCATABLE :: aaod_ocff(:, :, :)
  !                      Fossil-fuel org carb absorption aerosol optical depth
  REAL, ALLOCATABLE :: aaod_nitrate(:, :, :)
  !                      Nitrate absorption aerosol optical depth
  REAL, ALLOCATABLE :: aod_ukca_ait_sol(:, :, :)
  !                      UKCA Aitken soluble mode optical depth
  REAL, ALLOCATABLE :: aod_ukca_acc_sol(:, :, :)
  !                      UKCA accum. soluble mode optical depth
  REAL, ALLOCATABLE :: aod_ukca_cor_sol(:, :, :)
  !                      UKCA coarse soluble mode optical depth
  REAL, ALLOCATABLE :: aod_ukca_ait_ins(:, :, :)
  !                      UKCA Aitken insoluble mode optical depth
  REAL, ALLOCATABLE :: aod_ukca_acc_ins(:, :, :)
  !                      UKCA accum. insoluble mode optical depth
  REAL, ALLOCATABLE :: aod_ukca_cor_ins(:, :, :)
  !                      UKCA coarse insoluble mode optical depth
  REAL, ALLOCATABLE :: sod_ukca_ait_sol(:, :, :)
  !                      UKCA Aitken soluble mode stratospheric optical depth
  REAL, ALLOCATABLE :: sod_ukca_acc_sol(:, :, :)
  !                      UKCA accum. soluble mode stratospheric optical depth
  REAL, ALLOCATABLE :: sod_ukca_cor_sol(:, :, :)
  !                      UKCA coarse soluble mode stratospheric optical depth
  REAL, ALLOCATABLE :: sod_ukca_ait_ins(:, :, :)
  !                      UKCA Aitken insoluble mode stratospheric optical depth
  REAL, ALLOCATABLE :: sod_ukca_acc_ins(:, :, :)
  !                      UKCA accum. insoluble mode stratospheric optical depth
  REAL, ALLOCATABLE :: sod_ukca_cor_ins(:, :, :)
  !                      UKCA coarse insoluble mode stratospheric optical depth
  REAL, ALLOCATABLE :: ls_qcl_rad_path(:, :)
  !                      Large-scale liquid water path
  REAL, ALLOCATABLE :: ls_qcf_rad_path(:, :)
  !                      Large-scale ice water path
  REAL, ALLOCATABLE :: cc_qcl_rad_path(:, :)
  !                      Convective liquid water path
  REAL, ALLOCATABLE :: cc_qcf_rad_path(:, :)
  !                      Convective ice water path
  REAL, ALLOCATABLE :: ccore_qcl_rad_path(:, :)
  !                      Convective core liquid water path
  REAL, ALLOCATABLE :: ccore_qcf_rad_path(:, :)
  !                      Convective core ice water path
  REAL, ALLOCATABLE :: aaod_ukca_ait_sol(:, :, :)
  !                      UKCA Aitken soluble mode optical depth
  REAL, ALLOCATABLE :: aaod_ukca_acc_sol(:, :, :)
  !                      UKCA accum. soluble mode optical depth
  REAL, ALLOCATABLE :: aaod_ukca_cor_sol(:, :, :)
  !                      UKCA coarse soluble mode optical depth
  REAL, ALLOCATABLE :: aaod_ukca_ait_ins(:, :, :)
  !                      UKCA Aitken insoluble mode optical depth
  REAL, ALLOCATABLE :: aaod_ukca_acc_ins(:, :, :)
  !                      UKCA accum. insoluble mode optical depth
  REAL, ALLOCATABLE :: aaod_ukca_cor_ins(:, :, :)
  !                      UKCA coarse insoluble mode optical depth
  REAL, ALLOCATABLE :: clas_aerosol_ext(:,:,:,:)
  !                      CLASSIC 3D aerosol extinction
  REAL, ALLOCATABLE :: clas_aerosol_abs(:,:,:,:)
  !                      CLASSIC 3D aerosol absorption
  REAL, ALLOCATABLE :: clas_aerosol_sca(:,:,:,:)
  !                      CLASSIC 3D aerosol scattering
  REAL, ALLOCATABLE :: ukca_aerosol_ext(:,:,:,:)
  !                      UKCA 3D aerosol extinction
  REAL, ALLOCATABLE :: ukca_aerosol_abs(:,:,:,:)
  !                      UKCA 3D aerosol absorption
  REAL, ALLOCATABLE :: ukca_aerosol_sca(:,:,:,:)
  !                      UKCA 3D aerosol scattering
  REAL, ALLOCATABLE :: ukca_aerosol_gsca(:,:,:,:)
  !                      UKCA 3D aerosol asymmetry * scattering

END TYPE StrDiag


  CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEF_DIAG'

CONTAINS
!------------------------------------------------------------------------------
SUBROUTINE init_diag_logic(diag, spectrum)

USE def_spectrum, ONLY: StrSpecData

IMPLICIT NONE

TYPE (StrDiag)    , INTENT(INOUT) :: diag
TYPE (StrSpecData), INTENT(IN)    :: spectrum

! SW and LW
  diag%l_flux_up                          = .FALSE.
  diag%l_flux_down                        = .FALSE.
  diag%l_flux_up_clear                    = .FALSE.
  diag%l_flux_down_clear                  = .FALSE.
  diag%l_flux_up_clean                    = .FALSE.
  diag%l_flux_down_clean                  = .FALSE.
  diag%l_flux_up_clear_clean              = .FALSE.
  diag%l_flux_down_clear_clean            = .FALSE.
  diag%l_flux_up_band                     = .FALSE.
  diag%l_flux_down_band                   = .FALSE.
  diag%l_flux_up_clear_band               = .FALSE.
  diag%l_flux_down_clear_band             = .FALSE.
  diag%l_flux_up_clean_band               = .FALSE.
  diag%l_flux_down_clean_band             = .FALSE.
  diag%l_flux_up_clear_clean_band         = .FALSE.
  diag%l_flux_down_clear_clean_band       = .FALSE.
  diag%l_flux_up_forc                     = .FALSE.
  diag%l_flux_down_forc                   = .FALSE.
  diag%l_flux_up_clear_forc               = .FALSE.
  diag%l_flux_down_clear_forc             = .FALSE.
  diag%l_flux_up_forc_band                = .FALSE.
  diag%l_flux_down_forc_band              = .FALSE.
  diag%l_flux_up_clear_forc_band          = .FALSE.
  diag%l_flux_down_clear_forc_band        = .FALSE.
  diag%l_surf_down_clr                    = .FALSE.
  diag%l_clear_hr                         = .FALSE.
  diag%l_net_flux_trop                    = .FALSE.
  diag%l_toa_radiance                     = .FALSE.
  diag%l_aerosol_optical_depth            = .FALSE.
  diag%l_aerosol_scat_optical_depth       = .FALSE.
  diag%l_aerosol_asymmetry_scat           = .FALSE.
  diag%l_rad_mask                         = .FALSE.
  diag%l_emission_spectrum                = .FALSE.
  diag%l_emission_spectrum_clear          = .FALSE.
  diag%l_emission_spectrum_clean          = .FALSE.
  diag%l_emission_spectrum_clear_clean    = .FALSE.
  diag%l_toa_clear_weighted               = .FALSE.
  diag%l_total_clear_area                 = .FALSE.
  diag%l_easyaerosol_extinction           = .FALSE.
  diag%l_easyaerosol_absorption           = .FALSE.
  diag%l_easyaerosol_scattering           = .FALSE.
  diag%l_easyaerosol_asytimscat           = .FALSE.

! SW
  diag%l_flux_direct_sph                  = .FALSE.
  diag%l_flux_direct_div                  = .FALSE.
  diag%l_flux_direct_clear_sph            = .FALSE.
  diag%l_flux_direct_clear_div            = .FALSE.
  diag%l_flux_direct_clean_sph            = .FALSE.
  diag%l_flux_direct_clean_div            = .FALSE.
  diag%l_flux_direct_clear_clean_sph      = .FALSE.
  diag%l_flux_direct_clear_clean_div      = .FALSE.
  diag%l_flux_direct_sph_band             = .FALSE.
  diag%l_flux_direct_div_band             = .FALSE.
  diag%l_flux_direct_clear_sph_band       = .FALSE.
  diag%l_flux_direct_clear_div_band       = .FALSE.
  diag%l_flux_direct_clean_sph_band       = .FALSE.
  diag%l_flux_direct_clean_div_band       = .FALSE.
  diag%l_flux_direct_clear_clean_sph_band = .FALSE.
  diag%l_flux_direct_clear_clean_div_band = .FALSE.
  diag%l_transmission_spectrum            = .FALSE.
  diag%l_transmission_spectrum_clear      = .FALSE.
  diag%l_transmission_spectrum_clean      = .FALSE.
  diag%l_transmission_spectrum_clear_clean= .FALSE.
  diag%l_flux_direct_sph_forc             = .FALSE.
  diag%l_flux_direct_div_forc             = .FALSE.
  diag%l_flux_direct_clear_sph_forc       = .FALSE.
  diag%l_flux_direct_clear_div_forc       = .FALSE.
  diag%l_flux_direct_sph_forc_band        = .FALSE.
  diag%l_flux_direct_div_forc_band        = .FALSE.
  diag%l_flux_direct_clear_sph_forc_band  = .FALSE.
  diag%l_flux_direct_clear_div_forc_band  = .FALSE.
  diag%l_spherical_path                   = .FALSE.
  diag%l_solar_out_toa                    = .FALSE.
  diag%l_solar_out_clear                  = .FALSE.
  diag%l_surface_down_flux                = .FALSE.
  diag%l_surf_up_clr                      = .FALSE.
  diag%l_up_flux_trop                     = .FALSE.
  diag%l_flux_direct                      = .FALSE.
  diag%l_flux_diffuse                     = .FALSE.
  diag%re_conv_flag                       = .FALSE.
  diag%re_strat_flag                      = .FALSE.
  diag%wgt_conv_flag                      = .FALSE.
  diag%wgt_strat_flag                     = .FALSE.
  diag%lwp_strat_flag                     = .FALSE.
  diag%weighted_re_flag                   = .FALSE.
  diag%sum_weight_re_flag                 = .FALSE.
  diag%wgtd_warm_re_flag                  = .FALSE.
  diag%sum_wgt_warm_re_flag               = .FALSE.
  diag%cdnc_ct_diag_flag                  = .FALSE.
  diag%cdnc_ct_weight_flag                = .FALSE.
  diag%ntot_diag_flag                     = .FALSE.
  diag%strat_lwc_diag_flag                = .FALSE.
  diag%so4_ccn_diag_flag                  = .FALSE.
  diag%cond_samp_wgt_flag                 = .FALSE.
  diag%seasalt_film_flag                  = .FALSE.
  diag%seasalt_jet_flag                   = .FALSE.
  diag%nc_diag_flag                       = .FALSE.
  diag%nc_weight_flag                     = .FALSE.
  diag%l_FlxSolBelow690nmSurf             = .FALSE.
  diag%l_FlxSeaBelow690nmSurf             = .FALSE.
  diag%l_cloud_extinction                 = .FALSE.
  diag%l_cloud_weight_extinction          = .FALSE.
  diag%l_ls_cloud_extinction              = .FALSE.
  diag%l_ls_cloud_weight_extinction       = .FALSE.
  diag%l_cnv_cloud_extinction             = .FALSE.
  diag%l_cnv_cloud_weight_extinction      = .FALSE.
  diag%l_orog_corr                        = .FALSE.
  diag%l_uvflux_direct                    = .FALSE.
  diag%l_uvflux_up                        = .FALSE.
  diag%l_uvflux_down                      = .FALSE.
  diag%l_surf_uv                          = .FALSE.
  diag%l_surf_uv_clr                      = .FALSE.
  diag%l_direct_albedo                    = .FALSE.
  diag%l_diffuse_albedo                   = .FALSE.
  diag%l_vis_albedo_sc                    = .FALSE.
  diag%l_nir_albedo_sc                    = .FALSE.

! LW
  diag%l_total_cloud_cover                = .FALSE.
  diag%l_clear_olr                        = .FALSE.
  diag%l_down_flux_trop                   = .FALSE.
  diag%l_total_cloud_on_levels            = .FALSE.
  diag%l_cloud_absorptivity               = .FALSE.
  diag%l_cloud_weight_absorptivity        = .FALSE.
  diag%l_ls_cloud_absorptivity            = .FALSE.
  diag%l_ls_cloud_weight_absorptivity     = .FALSE.
  diag%l_cnv_cloud_absorptivity           = .FALSE.
  diag%l_cnv_cloud_weight_absorptivity    = .FALSE.
  diag%l_ls_qcl_rad                       = .FALSE.
  diag%l_ls_qcf_rad                       = .FALSE.
  diag%l_cc_qcl_rad                       = .FALSE.
  diag%l_cc_qcf_rad                       = .FALSE.
  diag%l_ls_cl_rad                        = .FALSE.
  diag%l_ls_cf_rad                        = .FALSE.
  diag%l_cc_cl_rad                        = .FALSE.
  diag%l_cc_cf_rad                        = .FALSE.
  diag%l_ccore_clt_rad                    = .FALSE.
  diag%l_ccore_qcl_rad                    = .FALSE.
  diag%l_ccore_qcf_rad                    = .FALSE.
  diag%l_ls_qcl_rad_path                  = .FALSE.
  diag%l_ls_qcf_rad_path                  = .FALSE.
  diag%l_cc_qcl_rad_path                  = .FALSE.
  diag%l_cc_qcf_rad_path                  = .FALSE.
  diag%l_ccore_qcl_rad_path               = .FALSE.
  diag%l_ccore_qcf_rad_path               = .FALSE.
  diag%l_ls_del_rad                       = .FALSE.
  diag%l_ls_def_rad                       = .FALSE.
  diag%l_cc_del_rad                       = .FALSE.
  diag%l_cc_def_rad                       = .FALSE.
  diag%l_aod_sulphate                     = .FALSE.
  diag%l_aod_dust                         = .FALSE.
  diag%l_aod_seasalt                      = .FALSE.
  diag%l_aod_soot                         = .FALSE.
  diag%l_aod_biomass                      = .FALSE.
  diag%l_aod_biogenic                     = .FALSE.
  diag%l_aod_ocff                         = .FALSE.
  diag%l_aod_delta                        = .FALSE.
  diag%l_aod_nitrate                      = .FALSE.
  diag%l_aod_total_radn                   = .FALSE.
  diag%l_angst_total_radn                 = .FALSE.
  diag%l_aod_prog_sulphate                = .FALSE.
  diag%l_aod_prog_dust                    = .FALSE.
  diag%l_aod_prog_seasalt                 = .FALSE.
  diag%l_aod_prog_soot                    = .FALSE.
  diag%l_aod_prog_biomass                 = .FALSE.
  diag%l_aod_prog_ocff                    = .FALSE.
  diag%l_aod_prog_nitrate                 = .FALSE.
  diag%l_aaod_sulphate                    = .FALSE.
  diag%l_aaod_dust                        = .FALSE.
  diag%l_aaod_seasalt                     = .FALSE.
  diag%l_aaod_soot                        = .FALSE.
  diag%l_aaod_biomass                     = .FALSE.
  diag%l_aaod_biogenic                    = .FALSE.
  diag%l_aaod_ocff                        = .FALSE.
  diag%l_aaod_nitrate                     = .FALSE.
  diag%l_aod_ukca_ait_sol                 = .FALSE.
  diag%l_aod_ukca_acc_sol                 = .FALSE.
  diag%l_aod_ukca_cor_sol                 = .FALSE.
  diag%l_aod_ukca_ait_ins                 = .FALSE.
  diag%l_aod_ukca_acc_ins                 = .FALSE.
  diag%l_aod_ukca_cor_ins                 = .FALSE.
  diag%l_sod_ukca_ait_sol                 = .FALSE.
  diag%l_sod_ukca_acc_sol                 = .FALSE.
  diag%l_sod_ukca_cor_sol                 = .FALSE.
  diag%l_sod_ukca_ait_ins                 = .FALSE.
  diag%l_sod_ukca_acc_ins                 = .FALSE.
  diag%l_sod_ukca_cor_ins                 = .FALSE.
  diag%l_aaod_ukca_ait_sol                = .FALSE.
  diag%l_aaod_ukca_acc_sol                = .FALSE.
  diag%l_aaod_ukca_cor_sol                = .FALSE.
  diag%l_aaod_ukca_ait_ins                = .FALSE.
  diag%l_aaod_ukca_acc_ins                = .FALSE.
  diag%l_aaod_ukca_cor_ins                = .FALSE.
IF (.NOT. ALLOCATED(diag%l_clas_aerosol_ext)) &
  ALLOCATE( diag%l_clas_aerosol_ext  (spectrum%dim%nd_aod_wavel))
IF (.NOT. ALLOCATED(diag%l_clas_aerosol_abs)) &
  ALLOCATE( diag%l_clas_aerosol_abs  (spectrum%dim%nd_aod_wavel))
IF (.NOT. ALLOCATED(diag%l_clas_aerosol_sca)) &
  ALLOCATE( diag%l_clas_aerosol_sca  (spectrum%dim%nd_aod_wavel))

IF (.NOT. ALLOCATED(diag%l_ukca_aerosol_ext)) &
  ALLOCATE( diag%l_ukca_aerosol_ext  (spectrum%dim%nd_aod_wavel))
IF (.NOT. ALLOCATED(diag%l_ukca_aerosol_abs)) &
  ALLOCATE( diag%l_ukca_aerosol_abs  (spectrum%dim%nd_aod_wavel))
IF (.NOT. ALLOCATED(diag%l_ukca_aerosol_sca)) &
  ALLOCATE( diag%l_ukca_aerosol_sca  (spectrum%dim%nd_aod_wavel))
IF (.NOT. ALLOCATED(diag%l_ukca_aerosol_gsca)) &
  ALLOCATE( diag%l_ukca_aerosol_gsca (spectrum%dim%nd_aod_wavel))

  diag%l_clas_aerosol_ext(:)              = .FALSE.
  diag%l_clas_aerosol_abs(:)              = .FALSE.
  diag%l_clas_aerosol_sca(:)              = .FALSE.
  diag%l_ukca_aerosol_ext(:)              = .FALSE.
  diag%l_ukca_aerosol_abs(:)              = .FALSE.
  diag%l_ukca_aerosol_sca(:)              = .FALSE.
  diag%l_ukca_aerosol_gsca(:)             = .FALSE.
  diag%n_clas_aerosol_ext                 = 0
  diag%n_clas_aerosol_abs                 = 0
  diag%n_clas_aerosol_sca                 = 0
  diag%n_ukca_aerosol_ext                 = 0
  diag%n_ukca_aerosol_abs                 = 0
  diag%n_ukca_aerosol_sca                 = 0
  diag%n_ukca_aerosol_gsca                = 0

END SUBROUTINE init_diag_logic
!------------------------------------------------------------------------------
SUBROUTINE allocate_diag(diag, spectrum,                                       &
                         row_length, rows, model_levels, cloud_levels, ntiles)

! If a radiation diagnostic is required the correct amount of memory space is
! allocated, and the diagnostic initialised to zero.

USE def_spectrum, ONLY: StrSpecData
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE missing_data_mod, ONLY: rmdi

IMPLICIT NONE

TYPE (StrDiag),     INTENT(INOUT) :: diag
TYPE (StrSpecData), INTENT(IN)    :: spectrum

INTEGER, INTENT(IN) :: row_length       ! Length of rows
INTEGER, INTENT(IN) :: rows             ! Number of rows
INTEGER, INTENT(IN) :: model_levels     ! Number of model levels
INTEGER, INTENT(IN) :: cloud_levels     ! Number of cloud levels
INTEGER, INTENT(IN) :: ntiles           ! Number of land surface tiles

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ALLOCATE_DIAG'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! SW and LW
IF (diag%l_flux_up) THEN
  ALLOCATE(diag%flux_up(row_length, rows, model_levels+1))
  diag%flux_up = 0.0
END IF

IF (diag%l_flux_down) THEN
  ALLOCATE(diag%flux_down(row_length, rows, model_levels+1))
  diag%flux_down = 0.0
END IF

IF (diag%l_flux_up_clear) THEN
  ALLOCATE(diag%flux_up_clear(row_length, rows, model_levels+1))
  diag%flux_up_clear = 0.0
END IF

IF (diag%l_flux_down_clear) THEN
  ALLOCATE(diag%flux_down_clear(row_length, rows, model_levels+1))
  diag%flux_down_clear = 0.0
END IF

IF (diag%l_flux_up_clean) THEN
  ALLOCATE(diag%flux_up_clean(row_length, rows, model_levels+1))
  diag%flux_up_clean = 0.0
END IF

IF (diag%l_flux_down_clean) THEN
  ALLOCATE(diag%flux_down_clean(row_length, rows, model_levels+1))
  diag%flux_down_clean = 0.0
END IF

IF (diag%l_flux_up_clear_clean) THEN
  ALLOCATE(diag%flux_up_clear_clean(row_length, rows, model_levels+1))
  diag%flux_up_clear_clean = 0.0
END IF

IF (diag%l_flux_down_clear_clean) THEN
  ALLOCATE(diag%flux_down_clear_clean(row_length, rows, model_levels+1))
  diag%flux_down_clear_clean = 0.0
END IF

IF (diag%l_flux_up_band) THEN
  ALLOCATE(diag%flux_up_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_up_band = 0.0
END IF

IF (diag%l_flux_down_band) THEN
  ALLOCATE(diag%flux_down_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_down_band = 0.0
END IF

IF (diag%l_flux_up_clear_band) THEN
  ALLOCATE(diag%flux_up_clear_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_up_clear_band = 0.0
END IF

IF (diag%l_flux_down_clear_band) THEN
  ALLOCATE(diag%flux_down_clear_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_down_clear_band = 0.0
END IF

IF (diag%l_flux_up_clean_band) THEN
  ALLOCATE(diag%flux_up_clean_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_up_clean_band = 0.0
END IF

IF (diag%l_flux_down_clean_band) THEN
  ALLOCATE(diag%flux_down_clean_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_down_clean_band = 0.0
END IF

IF (diag%l_flux_up_clear_clean_band) THEN
  ALLOCATE(diag%flux_up_clear_clean_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_up_clear_clean_band = 0.0
END IF

IF (diag%l_flux_down_clear_clean_band) THEN
  ALLOCATE(diag%flux_down_clear_clean_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_down_clear_clean_band = 0.0
END IF

IF (diag%l_flux_up_forc) THEN
  ALLOCATE(diag%flux_up_forc(row_length, rows, model_levels+1))
  diag%flux_up_forc = 0.0
END IF

IF (diag%l_flux_down_forc) THEN
  ALLOCATE(diag%flux_down_forc(row_length, rows, model_levels+1))
  diag%flux_down_forc = 0.0
END IF

IF (diag%l_flux_up_clear_forc) THEN
  ALLOCATE(diag%flux_up_clear_forc(row_length, rows, model_levels+1))
  diag%flux_up_clear_forc = 0.0
END IF

IF (diag%l_flux_down_clear_forc) THEN
  ALLOCATE(diag%flux_down_clear_forc(row_length, rows, model_levels+1))
  diag%flux_down_clear_forc = 0.0
END IF

IF (diag%l_flux_up_forc_band) THEN
  ALLOCATE(diag%flux_up_forc_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_up_forc_band = 0.0
END IF

IF (diag%l_flux_down_forc_band) THEN
  ALLOCATE(diag%flux_down_forc_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_down_forc_band = 0.0
END IF

IF (diag%l_flux_up_clear_forc_band) THEN
  ALLOCATE(diag%flux_up_clear_forc_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_up_clear_forc_band = 0.0
END IF

IF (diag%l_flux_down_clear_forc_band) THEN
  ALLOCATE(diag%flux_down_clear_forc_band(row_length, rows, model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_down_clear_forc_band = 0.0
END IF

IF (diag%l_surf_down_clr) THEN
  ALLOCATE(diag%surf_down_clr(row_length, rows))
  diag%surf_down_clr = 0.0
END IF

IF (diag%l_clear_hr) THEN
  ALLOCATE(diag%clear_hr(row_length, rows, model_levels))
  diag%clear_hr = 0.0
END IF

IF (diag%l_net_flux_trop) THEN
  ALLOCATE(diag%net_flux_trop(row_length, rows))
  diag%net_flux_trop = 0.0
END IF

IF (diag%l_toa_radiance) THEN
  ALLOCATE(diag%toa_radiance(row_length, rows, spectrum%dim%nd_band))
  diag%toa_radiance = 0.0
END IF

IF (diag%l_aerosol_optical_depth) THEN
  ALLOCATE(diag%aerosol_optical_depth(row_length, rows, model_levels, &
                                                  spectrum%dim%nd_band))
  diag%aerosol_optical_depth = 0.0
END IF

IF (diag%l_aerosol_scat_optical_depth) THEN
  ALLOCATE(diag%aerosol_scat_optical_depth(row_length, rows, model_levels, &
                                                  spectrum%dim%nd_band))
  diag%aerosol_scat_optical_depth = 0.0
END IF

IF (diag%l_aerosol_asymmetry_scat) THEN
  ALLOCATE(diag%aerosol_asymmetry_scat(row_length, rows, model_levels, &
                                                  spectrum%dim%nd_band))
  diag%aerosol_asymmetry_scat = 0.0
END IF

IF (diag%l_rad_mask) THEN
  ALLOCATE(diag%rad_mask(row_length, rows))
  diag%rad_mask = 0.0
END IF

IF (diag%l_emission_spectrum) THEN
  ALLOCATE(diag%emission_spectrum(row_length, rows, spectrum%dim%nd_band))
  diag%emission_spectrum = 0.0
END IF

IF (diag%l_emission_spectrum_clear) THEN
  ALLOCATE(diag%emission_spectrum_clear(row_length, rows, &
                                                  spectrum%dim%nd_band))
  diag%emission_spectrum_clear = 0.0
END IF

IF (diag%l_emission_spectrum_clean) THEN
  ALLOCATE(diag%emission_spectrum_clean(row_length, rows, &
                                                  spectrum%dim%nd_band))
  diag%emission_spectrum_clean = 0.0
END IF

IF (diag%l_emission_spectrum_clear_clean) THEN
  ALLOCATE(diag%emission_spectrum_clear_clean(row_length, rows, &
                                                  spectrum%dim%nd_band))
  diag%emission_spectrum_clear_clean = 0.0
END IF
IF (diag%l_toa_clear_weighted) THEN
  ALLOCATE(diag%toa_clear_weighted(row_length, rows))
  diag%toa_clear_weighted = 0.0
END IF

IF (diag%l_total_clear_area) THEN
  ALLOCATE(diag%total_clear_area(row_length, rows))
  diag%total_clear_area = 0.0
END IF

IF (diag%l_easyaerosol_extinction) THEN
  ALLOCATE(diag%easyaerosol_extinction(row_length, rows, model_levels, &
                                       spectrum%dim%nd_band))
  diag%easyaerosol_extinction = 0.0
END IF

IF (diag%l_easyaerosol_absorption) THEN
  ALLOCATE(diag%easyaerosol_absorption(row_length, rows, model_levels, &
                                       spectrum%dim%nd_band))
  diag%easyaerosol_absorption = 0.0
END IF

IF (diag%l_easyaerosol_scattering) THEN
  ALLOCATE(diag%easyaerosol_scattering(row_length, rows, model_levels, &
                                       spectrum%dim%nd_band))
  diag%easyaerosol_scattering = 0.0
END IF

IF (diag%l_easyaerosol_asytimscat) THEN
  ALLOCATE(diag%easyaerosol_asytimscat(row_length, rows, model_levels, &
                                       spectrum%dim%nd_band))
  diag%easyaerosol_asytimscat = 0.0
END IF

! SW
IF (diag%l_flux_direct_sph) THEN
  ALLOCATE(diag%flux_direct_sph(row_length, rows, 0:model_levels+1))
  diag%flux_direct_sph = 0.0
END IF

IF (diag%l_flux_direct_div) THEN
  ALLOCATE(diag%flux_direct_div(row_length, rows, model_levels))
  diag%flux_direct_div = 0.0
END IF

IF (diag%l_flux_direct_clear_sph) THEN
  ALLOCATE(diag%flux_direct_clear_sph(row_length, rows, 0:model_levels+1))
  diag%flux_direct_clear_sph = 0.0
END IF

IF (diag%l_flux_direct_clear_div) THEN
  ALLOCATE(diag%flux_direct_clear_div(row_length, rows, model_levels))
  diag%flux_direct_clear_div = 0.0
END IF

IF (diag%l_flux_direct_clean_sph) THEN
  ALLOCATE(diag%flux_direct_clean_sph(row_length, rows, 0:model_levels+1))
  diag%flux_direct_clean_sph = 0.0
END IF

IF (diag%l_flux_direct_clean_div) THEN
  ALLOCATE(diag%flux_direct_clean_div(row_length, rows, model_levels))
  diag%flux_direct_clean_div = 0.0
END IF

IF (diag%l_flux_direct_clear_clean_sph) THEN
  ALLOCATE(diag%flux_direct_clear_clean_sph(row_length, rows, 0:model_levels+1))
  diag%flux_direct_clear_clean_sph = 0.0
END IF

IF (diag%l_flux_direct_clear_clean_div) THEN
  ALLOCATE(diag%flux_direct_clear_clean_div(row_length, rows, model_levels))
  diag%flux_direct_clear_clean_div = 0.0
END IF

IF (diag%l_flux_direct_sph_band) THEN
  ALLOCATE(diag%flux_direct_sph_band(row_length, rows, 0:model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_sph_band = 0.0
END IF

IF (diag%l_flux_direct_div_band) THEN
  ALLOCATE(diag%flux_direct_div_band(row_length, rows, model_levels, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_div_band = 0.0
END IF

IF (diag%l_flux_direct_clear_sph_band) THEN
  ALLOCATE(diag%flux_direct_clear_sph_band(row_length, rows, 0:model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clear_sph_band = 0.0
END IF

IF (diag%l_flux_direct_clear_div_band) THEN
  ALLOCATE(diag%flux_direct_clear_div_band(row_length, rows, model_levels, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clear_div_band = 0.0
END IF

IF (diag%l_flux_direct_clean_sph_band) THEN
  ALLOCATE(diag%flux_direct_clean_sph_band(row_length, rows, 0:model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clean_sph_band = 0.0
END IF

IF (diag%l_flux_direct_clean_div_band) THEN
  ALLOCATE(diag%flux_direct_clean_div_band(row_length, rows, model_levels, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clean_div_band = 0.0
END IF

IF (diag%l_flux_direct_clear_clean_sph_band) THEN
  ALLOCATE(diag%flux_direct_clear_clean_sph_band( row_length, rows, &
                                                  0:model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clear_clean_sph_band = 0.0
END IF

IF (diag%l_flux_direct_clear_clean_div_band) THEN
  ALLOCATE(diag%flux_direct_clear_clean_div_band( row_length, rows, &
                                                  model_levels, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clear_clean_div_band = 0.0
END IF

IF (diag%l_transmission_spectrum) THEN
  ALLOCATE(diag%transmission_spectrum(row_length, rows, &
                                                  spectrum%dim%nd_band))
  diag%transmission_spectrum = 0.0
END IF

IF (diag%l_transmission_spectrum_clear) THEN
  ALLOCATE(diag%transmission_spectrum_clear(row_length, rows, &
                                                  spectrum%dim%nd_band))
  diag%transmission_spectrum_clear = 0.0
END IF

IF (diag%l_transmission_spectrum_clean) THEN
  ALLOCATE(diag%transmission_spectrum_clean(row_length, rows, &
                                                  spectrum%dim%nd_band))
  diag%transmission_spectrum_clean = 0.0
END IF

IF (diag%l_transmission_spectrum_clear_clean) THEN
  ALLOCATE(diag%transmission_spectrum_clear_clean(row_length, rows, &
                                                  spectrum%dim%nd_band))
  diag%transmission_spectrum_clear_clean = 0.0
END IF

IF (diag%l_flux_direct_sph_forc) THEN
  ALLOCATE(diag%flux_direct_sph_forc(row_length, rows, 0:model_levels+1))
  diag%flux_direct_sph_forc = 0.0
END IF

IF (diag%l_flux_direct_div_forc) THEN
  ALLOCATE(diag%flux_direct_div_forc(row_length, rows, model_levels))
  diag%flux_direct_div_forc = 0.0
END IF

IF (diag%l_flux_direct_clear_sph_forc) THEN
  ALLOCATE(diag%flux_direct_clear_sph_forc(row_length, rows, 0:model_levels+1))
  diag%flux_direct_clear_sph_forc = 0.0
END IF

IF (diag%l_flux_direct_clear_div_forc) THEN
  ALLOCATE(diag%flux_direct_clear_div_forc(row_length, rows, model_levels))
  diag%flux_direct_clear_div_forc = 0.0
END IF

IF (diag%l_flux_direct_sph_forc_band) THEN
  ALLOCATE(diag%flux_direct_sph_forc_band(row_length, rows, 0:model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_sph_forc_band = 0.0
END IF

IF (diag%l_flux_direct_div_forc_band) THEN
  ALLOCATE(diag%flux_direct_div_forc_band(row_length, rows, model_levels, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_div_forc_band = 0.0
END IF

IF (diag%l_flux_direct_clear_sph_forc_band) THEN
  ALLOCATE(diag%flux_direct_clear_sph_forc_band(row_length, rows, &
                                                  0:model_levels+1, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clear_sph_forc_band = 0.0
END IF

IF (diag%l_flux_direct_clear_div_forc_band) THEN
  ALLOCATE(diag%flux_direct_clear_div_forc_band(row_length, rows, &
                                                  model_levels, &
                                                  spectrum%dim%nd_band))
  diag%flux_direct_clear_div_forc_band = 0.0
END IF

IF (diag%l_spherical_path) THEN
  ALLOCATE(diag%spherical_path(row_length, rows, model_levels, &
                                                  0:model_levels+1))
  diag%spherical_path = 0.0
END IF

IF (diag%l_solar_out_toa) THEN
  ALLOCATE(diag%solar_out_toa(row_length, rows))
  diag%solar_out_toa = 0.0
END IF

IF (diag%l_solar_out_clear) THEN
  ALLOCATE(diag%solar_out_clear(row_length, rows))
  diag%solar_out_clear = 0.0
END IF

IF (diag%l_surface_down_flux) THEN
  ALLOCATE(diag%surface_down_flux(row_length, rows))
  diag%surface_down_flux = 0.0
END IF

IF (diag%l_surf_up_clr) THEN
  ALLOCATE(diag%surf_up_clr(row_length,rows))
  diag%surf_up_clr = 0.0
END IF

IF (diag%l_up_flux_trop) THEN
  ALLOCATE(diag%up_flux_trop(row_length, rows))
  diag%up_flux_trop = 0.0
END IF

IF (diag%l_flux_direct) THEN
  ALLOCATE(diag%flux_direct(row_length, rows, model_levels+1))
  diag%flux_direct = 0.0
END IF

IF (diag%l_flux_diffuse) THEN
  ALLOCATE(diag%flux_diffuse(row_length, rows, model_levels+1))
  diag%flux_diffuse = 0.0
END IF

IF (diag%l_uvflux_direct) THEN
  ALLOCATE(diag%uvflux_direct(row_length, rows, model_levels+1))
  diag%uvflux_direct = 0.0
END IF

IF (diag%l_uvflux_up) THEN
  ALLOCATE(diag%uvflux_up(row_length, rows, model_levels+1))
  diag%uvflux_up = 0.0
END IF

IF (diag%l_uvflux_down) THEN
  ALLOCATE(diag%uvflux_down(row_length, rows, model_levels+1))
  diag%uvflux_down = 0.0
END IF

IF (diag%l_surf_uv) THEN
  ALLOCATE(diag%surf_uv(row_length, rows))
  diag%surf_uv = 0.0
END IF

IF (diag%l_surf_uv_clr) THEN
  ALLOCATE(diag%surf_uv_clr(row_length, rows))
  diag%surf_uv_clr = 0.0
END IF

IF (diag%l_direct_albedo) THEN
  ALLOCATE(diag%direct_albedo(row_length, rows, spectrum%dim%nd_band))
  diag%direct_albedo = rmdi
END IF

IF (diag%l_diffuse_albedo) THEN
  ALLOCATE(diag%diffuse_albedo(row_length, rows, spectrum%dim%nd_band))
  diag%diffuse_albedo = rmdi
END IF

IF (diag%l_vis_albedo_sc) THEN
  ALLOCATE(diag%vis_albedo_sc(row_length, rows, ntiles))
  diag%vis_albedo_sc = rmdi
END IF

IF (diag%l_nir_albedo_sc) THEN
  ALLOCATE(diag%nir_albedo_sc(row_length, rows, ntiles))
  diag%nir_albedo_sc = rmdi
END IF

IF (diag%re_strat_flag) THEN
  ALLOCATE(diag%re_strat(row_length, rows, cloud_levels))
  diag%re_strat = 0.0
END IF

IF (diag%wgt_strat_flag) THEN
  ALLOCATE(diag%wgt_strat(row_length, rows, cloud_levels))
  diag%wgt_strat = 0.0
END IF

IF (diag%lwp_strat_flag) THEN
  ALLOCATE(diag%lwp_strat(row_length, rows, cloud_levels))
  diag%lwp_strat = 0.0
END IF

IF (diag%re_conv_flag) THEN
  ALLOCATE(diag%re_conv(row_length, rows, cloud_levels))
  diag%re_conv = 0.0
END IF

IF (diag%wgt_conv_flag) THEN
  ALLOCATE(diag%wgt_conv(row_length, rows,cloud_levels))
  diag%wgt_conv = 0.0
END IF

IF (diag%ntot_diag_flag) THEN
  ALLOCATE(diag%ntot_diag(row_length, rows, cloud_levels))
  diag%ntot_diag = 0.0
END IF

IF (diag%strat_lwc_diag_flag) THEN
  ALLOCATE(diag%strat_lwc_diag(row_length, rows, cloud_levels))
  diag%strat_lwc_diag = 0.0
END IF

IF (diag%so4_ccn_diag_flag) THEN
  ALLOCATE(diag%so4_ccn_diag(row_length, rows, cloud_levels))
  diag%so4_ccn_diag = 0.0
END IF

IF (diag%cond_samp_wgt_flag) THEN
  ALLOCATE(diag%cond_samp_wgt(row_length, rows, cloud_levels))
  diag%cond_samp_wgt = 0.0
END IF

IF (diag%weighted_re_flag) THEN
  ALLOCATE(diag%weighted_re(row_length, rows))
  diag%weighted_re = 0.0
END IF

IF (diag%sum_weight_re_flag) THEN
  ALLOCATE(diag%sum_weight_re(row_length, rows))
  diag%sum_weight_re = 0.0
END IF

IF (diag%wgtd_warm_re_flag) THEN
  ALLOCATE(diag%weighted_warm_re(row_length, rows))
  diag%weighted_warm_re = 0.0
END IF

IF (diag%sum_wgt_warm_re_flag) THEN
  ALLOCATE(diag%sum_weight_warm_re(row_length, rows))
  diag%sum_weight_warm_re = 0.0
END IF

IF (diag%cdnc_ct_diag_flag) THEN
  ALLOCATE(diag%cdnc_ct_diag(row_length, rows))
  diag%cdnc_ct_diag = 0.0
END IF

IF (diag%cdnc_ct_weight_flag) THEN
  ALLOCATE(diag%cdnc_ct_weight(row_length, rows))
  diag%cdnc_ct_weight = 0.0
END IF

IF (diag%nc_diag_flag) THEN
  ALLOCATE(diag%nc_diag(row_length, rows))
  diag%nc_diag = 0.0
END IF

IF (diag%nc_weight_flag) THEN
  ALLOCATE(diag%nc_weight(row_length, rows))
  diag%nc_weight = 0.0
END IF

IF (diag%l_FlxSolBelow690nmSurf) THEN
  ALLOCATE(diag%FlxSolBelow690nmSurf(row_length, rows))
  diag%FlxSolBelow690nmSurf = 0.0
END IF

IF (diag%l_FlxSeaBelow690nmSurf) THEN
  ALLOCATE(diag%FlxSeaBelow690nmSurf(row_length, rows))
  diag%FlxSeaBelow690nmSurf = 0.0
END IF

IF (diag%l_orog_corr) THEN
  ALLOCATE(diag%orog_corr(row_length, rows))
  diag%orog_corr = 1.0
END IF

IF (diag%l_cloud_extinction) THEN
  ALLOCATE(diag%cloud_extinction(row_length, rows, cloud_levels))
  diag%cloud_extinction = 0.0
END IF

IF (diag%l_cloud_weight_extinction) THEN
  ALLOCATE(diag%cloud_weight_extinction(row_length, rows, cloud_levels))
  diag%cloud_weight_extinction = 0.0
END IF

IF (diag%l_ls_cloud_extinction) THEN
  ALLOCATE(diag%ls_cloud_extinction(row_length, rows, cloud_levels))
  diag%ls_cloud_extinction = 0.0
END IF

IF (diag%l_ls_cloud_weight_extinction) THEN
  ALLOCATE(diag%ls_cloud_weight_extinction(row_length, rows, cloud_levels))
  diag%ls_cloud_weight_extinction = 0.0
END IF

IF (diag%l_cnv_cloud_extinction) THEN
  ALLOCATE(diag%cnv_cloud_extinction(row_length, rows, cloud_levels))
  diag%cnv_cloud_extinction = 0.0
END IF

IF (diag%l_cnv_cloud_weight_extinction) THEN
  ALLOCATE(diag%cnv_cloud_weight_extinction(row_length, rows, cloud_levels))
  diag%cnv_cloud_weight_extinction = 0.0
END IF


! LW
IF (diag%l_total_cloud_cover) THEN
  ALLOCATE(diag%total_cloud_cover(row_length, rows))
  diag%total_cloud_cover(:,:)=0.0
END IF

IF (diag%l_clear_olr) THEN
  ALLOCATE(diag%clear_olr(row_length, rows))
  diag%clear_olr(:,:)=0.0
END IF

IF (diag%l_down_flux_trop) THEN
  ALLOCATE(diag%down_flux_trop(row_length, rows))
  diag%down_flux_trop(:,:)=0.0
END IF

IF (diag%l_cloud_absorptivity) THEN
  ALLOCATE(diag%cloud_absorptivity(row_length, rows, cloud_levels))
  diag%cloud_absorptivity(:,:,:)=0.0
END IF

IF (diag%l_cloud_weight_absorptivity) THEN
  ALLOCATE(diag%cloud_weight_absorptivity(row_length, rows, cloud_levels))
  diag%cloud_weight_absorptivity(:,:,:)=0.0
END IF

IF (diag%l_ls_cloud_absorptivity) THEN
  ALLOCATE(diag%ls_cloud_absorptivity(row_length, rows, cloud_levels))
  diag%ls_cloud_absorptivity(:,:,:)=0.0
END IF

IF (diag%l_ls_cloud_weight_absorptivity) THEN
  ALLOCATE(diag%ls_cloud_weight_absorptivity(row_length, rows, cloud_levels))
  diag%ls_cloud_weight_absorptivity(:,:,:)=0.0
END IF

IF (diag%l_cnv_cloud_absorptivity) THEN
  ALLOCATE(diag%cnv_cloud_absorptivity(row_length, rows, cloud_levels))
  diag%cnv_cloud_absorptivity(:,:,:)=0.0
END IF

IF (diag%l_cnv_cloud_weight_absorptivity) THEN
  ALLOCATE(diag%cnv_cloud_weight_absorptivity(row_length, rows, cloud_levels))
  diag%cnv_cloud_weight_absorptivity(:,:,:)=0.0
END IF

IF (diag%l_total_cloud_on_levels) THEN
  ALLOCATE(diag%total_cloud_on_levels(row_length, rows, cloud_levels))
  diag%total_cloud_on_levels(:,:,:)=0.0
END IF

IF (diag%l_aod_sulphate) THEN
  ALLOCATE(diag%aod_sulphate(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_sulphate(:,:,:) = 0.0
END IF

IF (diag%l_aod_dust) THEN
  ALLOCATE(diag%aod_dust(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_dust(:,:,:) = 0.0
END IF

IF (diag%l_aod_seasalt) THEN
  ALLOCATE(diag%aod_seasalt(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_seasalt(:,:,:) = 0.0
END IF

IF (diag%l_aod_soot) THEN
  ALLOCATE(diag%aod_soot(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_soot(:,:,:) = 0.0
END IF

IF (diag%l_aod_biomass) THEN
  ALLOCATE(diag%aod_biomass(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_biomass(:,:,:) = 0.0
END IF

IF (diag%l_aod_biogenic) THEN
  ALLOCATE(diag%aod_biogenic(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_biogenic(:,:,:) = 0.0
END IF

IF (diag%l_aod_ocff) THEN
  ALLOCATE(diag%aod_ocff(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_ocff(:,:,:) = 0.0
END IF

IF (diag%l_aod_delta) THEN
  ALLOCATE(diag%aod_delta(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_delta(:,:,:) = 0.0
END IF

IF (diag%l_aod_nitrate) THEN
  ALLOCATE(diag%aod_nitrate(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_nitrate(:,:,:) = 0.0
END IF

IF (diag%l_aod_total_radn) THEN
  ALLOCATE(diag%aod_total_radn(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_total_radn(:,:,:) = 0.0
END IF

IF (diag%l_angst_total_radn) THEN
  ALLOCATE(diag%angst_total_radn(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%angst_total_radn(:,:,:) = 0.0
END IF

IF (diag%l_aod_prog_sulphate) THEN
  ALLOCATE(diag%aod_prog_sulphate(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_prog_sulphate(:,:,:) = 0.0
END IF

IF (diag%l_aod_prog_dust) THEN
  ALLOCATE(diag%aod_prog_dust(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_prog_dust(:,:,:) = 0.0
END IF

IF (diag%l_aod_prog_seasalt) THEN
  ALLOCATE(diag%aod_prog_seasalt(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_prog_seasalt(:,:,:) = 0.0
END IF

IF (diag%l_aod_prog_soot) THEN
  ALLOCATE(diag%aod_prog_soot(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_prog_soot(:,:,:) = 0.0
END IF

IF (diag%l_aod_prog_biomass) THEN
  ALLOCATE(diag%aod_prog_biomass(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_prog_biomass(:,:,:) = 0.0
END IF

IF (diag%l_aod_prog_ocff) THEN
  ALLOCATE(diag%aod_prog_ocff(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_prog_ocff(:,:,:) = 0.0
END IF

IF (diag%l_aod_prog_nitrate) THEN
  ALLOCATE(diag%aod_prog_nitrate(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_prog_nitrate(:,:,:) = 0.0
END IF

IF (diag%l_aaod_sulphate) THEN
  ALLOCATE(diag%aaod_sulphate(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_sulphate(:,:,:) = 0.0
END IF

IF (diag%l_aaod_dust) THEN
  ALLOCATE(diag%aaod_dust(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_dust(:,:,:) = 0.0
END IF

IF (diag%l_aaod_seasalt) THEN
  ALLOCATE(diag%aaod_seasalt(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_seasalt(:,:,:) = 0.0
END IF

IF (diag%l_aaod_soot) THEN
  ALLOCATE(diag%aaod_soot(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_soot(:,:,:) = 0.0
END IF

IF (diag%l_aaod_biomass) THEN
  ALLOCATE(diag%aaod_biomass(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_biomass(:,:,:) = 0.0
END IF

IF (diag%l_aaod_biogenic) THEN
  ALLOCATE(diag%aaod_biogenic(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_biogenic(:,:,:) = 0.0
END IF

IF (diag%l_aaod_ocff) THEN
  ALLOCATE(diag%aaod_ocff(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_ocff(:,:,:) = 0.0
END IF

IF (diag%l_aaod_nitrate) THEN
  ALLOCATE(diag%aaod_nitrate(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_nitrate(:,:,:) = 0.0
END IF

IF (diag%l_aod_ukca_ait_sol) THEN
  ALLOCATE(diag%aod_ukca_ait_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_ukca_ait_sol(:,:,:) = 0.0
END IF

IF (diag%l_aod_ukca_acc_sol) THEN
  ALLOCATE(diag%aod_ukca_acc_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_ukca_acc_sol(:,:,:) = 0.0
END IF

IF (diag%l_aod_ukca_cor_sol) THEN
  ALLOCATE(diag%aod_ukca_cor_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_ukca_cor_sol(:,:,:) = 0.0
END IF

IF (diag%l_aod_ukca_ait_ins) THEN
  ALLOCATE(diag%aod_ukca_ait_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_ukca_ait_ins(:,:,:) = 0.0
END IF

IF (diag%l_aod_ukca_acc_ins) THEN
  ALLOCATE(diag%aod_ukca_acc_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_ukca_acc_ins(:,:,:) = 0.0
END IF

IF (diag%l_aod_ukca_cor_ins) THEN
  ALLOCATE(diag%aod_ukca_cor_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aod_ukca_cor_ins(:,:,:) = 0.0
END IF

IF (diag%l_sod_ukca_ait_sol) THEN
  ALLOCATE(diag%sod_ukca_ait_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%sod_ukca_ait_sol(:,:,:) = 0.0
END IF

IF (diag%l_sod_ukca_acc_sol) THEN
  ALLOCATE(diag%sod_ukca_acc_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%sod_ukca_acc_sol(:,:,:) = 0.0
END IF

IF (diag%l_sod_ukca_cor_sol) THEN
  ALLOCATE(diag%sod_ukca_cor_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%sod_ukca_cor_sol(:,:,:) = 0.0
END IF

IF (diag%l_sod_ukca_ait_ins) THEN
  ALLOCATE(diag%sod_ukca_ait_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%sod_ukca_ait_ins(:,:,:) = 0.0
END IF

IF (diag%l_sod_ukca_acc_ins) THEN
  ALLOCATE(diag%sod_ukca_acc_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%sod_ukca_acc_ins(:,:,:) = 0.0
END IF

IF (diag%l_sod_ukca_cor_ins) THEN
  ALLOCATE(diag%sod_ukca_cor_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%sod_ukca_cor_ins(:,:,:) = 0.0
END IF

IF (diag%l_aaod_ukca_ait_sol) THEN
  ALLOCATE(diag%aaod_ukca_ait_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_ukca_ait_sol(:,:,:) = 0.0
END IF

IF (diag%l_aaod_ukca_acc_sol) THEN
  ALLOCATE(diag%aaod_ukca_acc_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_ukca_acc_sol(:,:,:) = 0.0
END IF

IF (diag%l_aaod_ukca_cor_sol) THEN
  ALLOCATE(diag%aaod_ukca_cor_sol(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_ukca_cor_sol(:,:,:) = 0.0
END IF

IF (diag%l_aaod_ukca_ait_ins) THEN
  ALLOCATE(diag%aaod_ukca_ait_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_ukca_ait_ins(:,:,:) = 0.0
END IF

IF (diag%l_aaod_ukca_acc_ins) THEN
  ALLOCATE(diag%aaod_ukca_acc_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_ukca_acc_ins(:,:,:) = 0.0
END IF

IF (diag%l_aaod_ukca_cor_ins) THEN
  ALLOCATE(diag%aaod_ukca_cor_ins(row_length, rows, spectrum%dim%nd_aod_wavel))
  diag%aaod_ukca_cor_ins(:,:,:) = 0.0
END IF

IF (diag%n_clas_aerosol_ext > 0) THEN   
  ALLOCATE(diag%clas_aerosol_ext(row_length, rows, model_levels,         &
       diag%n_clas_aerosol_ext))
  diag%clas_aerosol_ext(:,:,:,:) = 0.0
END IF

IF (diag%n_clas_aerosol_abs > 0) THEN   
  ALLOCATE(diag%clas_aerosol_abs(row_length, rows, model_levels,         &
       diag%n_clas_aerosol_abs))
  diag%clas_aerosol_abs(:,:,:,:) = 0.0
END IF

IF (diag%n_clas_aerosol_sca > 0) THEN   
  ALLOCATE(diag%clas_aerosol_sca(row_length, rows, model_levels,         &
       diag%n_clas_aerosol_sca))
  diag%clas_aerosol_sca(:,:,:,:) = 0.0
END IF

IF (diag%n_ukca_aerosol_ext > 0) THEN   
  ALLOCATE(diag%ukca_aerosol_ext(row_length, rows, model_levels,         &
       diag%n_ukca_aerosol_ext))
  diag%ukca_aerosol_ext(:,:,:,:) = 0.0
END IF

IF (diag%n_ukca_aerosol_abs > 0) THEN   
  ALLOCATE(diag%ukca_aerosol_abs(row_length, rows, model_levels,         &
       diag%n_ukca_aerosol_abs))
  diag%ukca_aerosol_abs(:,:,:,:) = 0.0
END IF

IF (diag%n_ukca_aerosol_sca > 0) THEN   
  ALLOCATE(diag%ukca_aerosol_sca(row_length, rows, model_levels,         &
       diag%n_ukca_aerosol_sca))
  diag%ukca_aerosol_sca(:,:,:,:) = 0.0
END IF

IF (diag%n_ukca_aerosol_gsca > 0) THEN   
  ALLOCATE(diag%ukca_aerosol_gsca(row_length, rows, model_levels,         &
       diag%n_ukca_aerosol_gsca))
  diag%ukca_aerosol_gsca(:,:,:,:) = 0.0
END IF

IF (diag%l_ls_qcl_rad) THEN
  ALLOCATE(diag%ls_qcl_rad(row_length, rows, model_levels))
  diag%ls_qcl_rad(:,:,:)=0.0
END IF

IF (diag%l_ls_qcf_rad) THEN
  ALLOCATE(diag%ls_qcf_rad(row_length, rows, model_levels))
  diag%ls_qcf_rad(:,:,:)=0.0
END IF

IF (diag%l_cc_qcl_rad) THEN
  ALLOCATE(diag%cc_qcl_rad(row_length, rows, model_levels))
  diag%cc_qcl_rad(:,:,:)=0.0
END IF

IF (diag%l_cc_qcf_rad) THEN
  ALLOCATE(diag%cc_qcf_rad(row_length, rows, model_levels))
  diag%cc_qcf_rad(:,:,:)=0.0
END IF

IF (diag%l_ls_cl_rad) THEN
  ALLOCATE(diag%ls_cl_rad(row_length, rows, model_levels))
  diag%ls_cl_rad(:,:,:)=0.0
END IF

IF (diag%l_ls_cf_rad) THEN
  ALLOCATE(diag%ls_cf_rad(row_length, rows, model_levels))
  diag%ls_cf_rad(:,:,:)=0.0
END IF

IF (diag%l_cc_cl_rad) THEN
  ALLOCATE(diag%cc_cl_rad(row_length, rows, model_levels))
  diag%cc_cl_rad(:,:,:)=0.0
END IF

IF (diag%l_cc_cf_rad) THEN
  ALLOCATE(diag%cc_cf_rad(row_length, rows, model_levels))
  diag%cc_cf_rad(:,:,:)=0.0
END IF

IF (diag%l_ccore_clt_rad) THEN
  ALLOCATE(diag%ccore_clt_rad(row_length, rows, model_levels))
  diag%ccore_clt_rad(:,:,:)=0.0
END IF

IF (diag%l_ccore_qcl_rad) THEN
  ALLOCATE(diag%ccore_qcl_rad(row_length, rows, model_levels))
  diag%ccore_qcl_rad(:,:,:)=0.0
END IF

IF (diag%l_ccore_qcf_rad) THEN
  ALLOCATE(diag%ccore_qcf_rad(row_length, rows, model_levels))
  diag%ccore_qcf_rad(:,:,:)=0.0
END IF

IF (diag%l_ls_qcl_rad_path) THEN
  ALLOCATE(diag%ls_qcl_rad_path(row_length, rows))
  diag%ls_qcl_rad_path(:,:)=0.0
END IF

IF (diag%l_ls_qcf_rad_path) THEN
  ALLOCATE(diag%ls_qcf_rad_path(row_length, rows))
  diag%ls_qcf_rad_path(:,:)=0.0
END IF

IF (diag%l_cc_qcl_rad_path) THEN
  ALLOCATE(diag%cc_qcl_rad_path(row_length, rows))
  diag%cc_qcl_rad_path(:,:)=0.0
END IF

IF (diag%l_cc_qcf_rad_path) THEN
  ALLOCATE(diag%cc_qcf_rad_path(row_length, rows))
  diag%cc_qcf_rad_path(:,:)=0.0
END IF

IF (diag%l_ccore_qcl_rad_path) THEN
  ALLOCATE(diag%ccore_qcl_rad_path(row_length, rows))
  diag%ccore_qcl_rad_path(:,:)=0.0
END IF

IF (diag%l_ccore_qcf_rad_path) THEN
  ALLOCATE(diag%ccore_qcf_rad_path(row_length, rows))
  diag%ccore_qcf_rad_path(:,:)=0.0
END IF

IF (diag%l_ls_del_rad) THEN
  ALLOCATE(diag%ls_del_rad(row_length, rows, model_levels))
  diag%ls_del_rad(:,:,:)=0.0
END IF

IF (diag%l_ls_def_rad) THEN
  ALLOCATE(diag%ls_def_rad(row_length, rows, model_levels))
  diag%ls_def_rad(:,:,:)=0.0
END IF

IF (diag%l_cc_del_rad) THEN
  ALLOCATE(diag%cc_del_rad(row_length, rows, model_levels))
  diag%cc_del_rad(:,:,:)=0.0
END IF

IF (diag%l_cc_def_rad) THEN
  ALLOCATE(diag%cc_def_rad(row_length, rows, model_levels))
  diag%cc_def_rad(:,:,:)=0.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE allocate_diag
!------------------------------------------------------------------------------
SUBROUTINE deallocate_diag(diag)

IMPLICIT NONE

TYPE (StrDiag), INTENT(INOUT) :: diag

! SW and LW
IF (ALLOCATED(diag%flux_up)) &
  DEALLOCATE(diag%flux_up)
IF (ALLOCATED(diag%flux_down)) &
  DEALLOCATE(diag%flux_down)
IF (ALLOCATED(diag%flux_up_clear)) &
  DEALLOCATE(diag%flux_up_clear)
IF (ALLOCATED(diag%flux_down_clear)) &
  DEALLOCATE(diag%flux_down_clear)
IF (ALLOCATED(diag%flux_up_clean)) &
  DEALLOCATE(diag%flux_up_clean)
IF (ALLOCATED(diag%flux_down_clean)) &
  DEALLOCATE(diag%flux_down_clean)
IF (ALLOCATED(diag%flux_up_clear_clean)) &
  DEALLOCATE(diag%flux_up_clear_clean)
IF (ALLOCATED(diag%flux_down_clear_clean)) &
  DEALLOCATE(diag%flux_down_clear_clean)
IF (ALLOCATED(diag%flux_up_band)) &
  DEALLOCATE(diag%flux_up_band)
IF (ALLOCATED(diag%flux_down_band)) &
  DEALLOCATE(diag%flux_down_band)
IF (ALLOCATED(diag%flux_up_clear_band)) &
  DEALLOCATE(diag%flux_up_clear_band)
IF (ALLOCATED(diag%flux_down_clear_band)) &
  DEALLOCATE(diag%flux_down_clear_band)
IF (ALLOCATED(diag%flux_up_clean_band)) &
  DEALLOCATE(diag%flux_up_clean_band)
IF (ALLOCATED(diag%flux_down_clean_band)) &
  DEALLOCATE(diag%flux_down_clean_band)
IF (ALLOCATED(diag%flux_up_clear_clean_band)) &
  DEALLOCATE(diag%flux_up_clear_clean_band)
IF (ALLOCATED(diag%flux_down_clear_clean_band)) &
  DEALLOCATE(diag%flux_down_clear_clean_band)
IF (ALLOCATED(diag%flux_up_forc)) &
  DEALLOCATE(diag%flux_up_forc)
IF (ALLOCATED(diag%flux_down_forc)) &
  DEALLOCATE(diag%flux_down_forc)
IF (ALLOCATED(diag%flux_up_clear_forc)) &
  DEALLOCATE(diag%flux_up_clear_forc)
IF (ALLOCATED(diag%flux_down_clear_forc)) &
  DEALLOCATE(diag%flux_down_clear_forc)
IF (ALLOCATED(diag%flux_up_forc_band)) &
  DEALLOCATE(diag%flux_up_forc_band)
IF (ALLOCATED(diag%flux_down_forc_band)) &
  DEALLOCATE(diag%flux_down_forc_band)
IF (ALLOCATED(diag%flux_up_clear_forc_band)) &
  DEALLOCATE(diag%flux_up_clear_forc_band)
IF (ALLOCATED(diag%flux_down_clear_forc_band)) &
  DEALLOCATE(diag%flux_down_clear_forc_band)
IF (ALLOCATED(diag%surf_down_clr)) &
  DEALLOCATE(diag%surf_down_clr)
IF (ALLOCATED(diag%clear_hr)) &
  DEALLOCATE(diag%clear_hr)
IF (ALLOCATED(diag%net_flux_trop)) &
  DEALLOCATE(diag%net_flux_trop)
IF (ALLOCATED(diag%toa_radiance)) &
  DEALLOCATE(diag%toa_radiance)
IF (ALLOCATED(diag%aerosol_optical_depth)) &
  DEALLOCATE(diag%aerosol_optical_depth)
IF (ALLOCATED(diag%aerosol_scat_optical_depth)) &
  DEALLOCATE(diag%aerosol_scat_optical_depth)
IF (ALLOCATED(diag%aerosol_asymmetry_scat)) &
  DEALLOCATE(diag%aerosol_asymmetry_scat)
IF (ALLOCATED(diag%rad_mask)) &
  DEALLOCATE(diag%rad_mask)
IF (ALLOCATED(diag%emission_spectrum)) &
  DEALLOCATE(diag%emission_spectrum)
IF (ALLOCATED(diag%emission_spectrum_clear)) &
  DEALLOCATE(diag%emission_spectrum_clear)
IF (ALLOCATED(diag%emission_spectrum_clean)) &
  DEALLOCATE(diag%emission_spectrum_clean)
IF (ALLOCATED(diag%emission_spectrum_clear_clean)) &
  DEALLOCATE(diag%emission_spectrum_clear_clean)
IF (ALLOCATED(diag%toa_clear_weighted)) &
  DEALLOCATE(diag%toa_clear_weighted)
IF (ALLOCATED(diag%total_clear_area)) &
  DEALLOCATE(diag%total_clear_area)
IF (ALLOCATED(diag%easyaerosol_extinction)) &
  DEALLOCATE(diag%easyaerosol_extinction)
IF (ALLOCATED(diag%easyaerosol_absorption)) &
  DEALLOCATE(diag%easyaerosol_absorption)
IF (ALLOCATED(diag%easyaerosol_scattering)) &
  DEALLOCATE(diag%easyaerosol_scattering)
IF (ALLOCATED(diag%easyaerosol_asytimscat)) &
  DEALLOCATE(diag%easyaerosol_asytimscat)

! SW
IF (ALLOCATED(diag%flux_direct_sph)) &
  DEALLOCATE(diag%flux_direct_sph)
IF (ALLOCATED(diag%flux_direct_div)) &
  DEALLOCATE(diag%flux_direct_div)
IF (ALLOCATED(diag%flux_direct_clear_sph)) &
  DEALLOCATE(diag%flux_direct_clear_sph)
IF (ALLOCATED(diag%flux_direct_clear_div)) &
  DEALLOCATE(diag%flux_direct_clear_div)
IF (ALLOCATED(diag%flux_direct_clean_sph)) &
  DEALLOCATE(diag%flux_direct_clean_sph)
IF (ALLOCATED(diag%flux_direct_clean_div)) &
  DEALLOCATE(diag%flux_direct_clean_div)
IF (ALLOCATED(diag%flux_direct_clear_clean_sph)) &
  DEALLOCATE(diag%flux_direct_clear_clean_sph)
IF (ALLOCATED(diag%flux_direct_clear_clean_div)) &
  DEALLOCATE(diag%flux_direct_clear_clean_div)
IF (ALLOCATED(diag%flux_direct_sph_band)) &
  DEALLOCATE(diag%flux_direct_sph_band)
IF (ALLOCATED(diag%flux_direct_div_band)) &
  DEALLOCATE(diag%flux_direct_div_band)
IF (ALLOCATED(diag%flux_direct_clear_sph_band)) &
  DEALLOCATE(diag%flux_direct_clear_sph_band)
IF (ALLOCATED(diag%flux_direct_clear_div_band)) &
  DEALLOCATE(diag%flux_direct_clear_div_band)
IF (ALLOCATED(diag%flux_direct_clean_sph_band)) &
  DEALLOCATE(diag%flux_direct_clean_sph_band)
IF (ALLOCATED(diag%flux_direct_clean_div_band)) &
  DEALLOCATE(diag%flux_direct_clean_div_band)
IF (ALLOCATED(diag%flux_direct_clear_clean_sph_band)) &
  DEALLOCATE(diag%flux_direct_clear_clean_sph_band)
IF (ALLOCATED(diag%flux_direct_clear_clean_div_band)) &
  DEALLOCATE(diag%flux_direct_clear_clean_div_band)
IF (ALLOCATED(diag%transmission_spectrum)) &
  DEALLOCATE(diag%transmission_spectrum)
IF (ALLOCATED(diag%transmission_spectrum_clear)) &
  DEALLOCATE(diag%transmission_spectrum_clear)
IF (ALLOCATED(diag%transmission_spectrum_clean)) &
  DEALLOCATE(diag%transmission_spectrum_clean)
IF (ALLOCATED(diag%transmission_spectrum_clear_clean)) &
  DEALLOCATE(diag%transmission_spectrum_clear_clean)
IF (ALLOCATED(diag%flux_direct_sph_forc)) &
  DEALLOCATE(diag%flux_direct_sph_forc)
IF (ALLOCATED(diag%flux_direct_div_forc)) &
  DEALLOCATE(diag%flux_direct_div_forc)
IF (ALLOCATED(diag%flux_direct_clear_sph_forc)) &
  DEALLOCATE(diag%flux_direct_clear_sph_forc)
IF (ALLOCATED(diag%flux_direct_clear_div_forc)) &
  DEALLOCATE(diag%flux_direct_clear_div_forc)
IF (ALLOCATED(diag%flux_direct_sph_forc_band)) &
  DEALLOCATE(diag%flux_direct_sph_forc_band)
IF (ALLOCATED(diag%flux_direct_div_forc_band)) &
  DEALLOCATE(diag%flux_direct_div_forc_band)
IF (ALLOCATED(diag%flux_direct_clear_sph_forc_band)) &
  DEALLOCATE(diag%flux_direct_clear_sph_forc_band)
IF (ALLOCATED(diag%flux_direct_clear_div_forc_band)) &
  DEALLOCATE(diag%flux_direct_clear_div_forc_band)
IF (ALLOCATED(diag%spherical_path)) &
  DEALLOCATE(diag%spherical_path)
IF (ALLOCATED(diag%solar_out_toa)) &
  DEALLOCATE(diag%solar_out_toa)
IF (ALLOCATED(diag%solar_out_clear)) &
  DEALLOCATE(diag%solar_out_clear)
IF (ALLOCATED(diag%surface_down_flux)) &
  DEALLOCATE(diag%surface_down_flux)
IF (ALLOCATED(diag%surf_up_clr)) &
  DEALLOCATE(diag%surf_up_clr)
IF (ALLOCATED(diag%up_flux_trop)) &
  DEALLOCATE(diag%up_flux_trop)
! UV-Fluxes
IF (ALLOCATED(diag%uvflux_direct)) &
  DEALLOCATE(diag%uvflux_direct)
IF (ALLOCATED(diag%uvflux_up)) &
  DEALLOCATE(diag%uvflux_up)
IF (ALLOCATED(diag%uvflux_down)) &
  DEALLOCATE(diag%uvflux_down)
IF (ALLOCATED(diag%surf_uv)) &
  DEALLOCATE(diag%surf_uv)
IF (ALLOCATED(diag%surf_uv_clr)) &
  DEALLOCATE(diag%surf_uv_clr)
! Surface Albedos
IF (ALLOCATED(diag%direct_albedo)) &
  DEALLOCATE(diag%direct_albedo)
IF (ALLOCATED(diag%diffuse_albedo)) &
  DEALLOCATE(diag%diffuse_albedo)
IF (ALLOCATED(diag%vis_albedo_sc)) &
  DEALLOCATE(diag%vis_albedo_sc)
IF (ALLOCATED(diag%nir_albedo_sc)) &
  DEALLOCATE(diag%nir_albedo_sc)
! Direct and diffuse downward SW Fluxes
IF (ALLOCATED(diag%flux_direct)) &
  DEALLOCATE(diag%flux_direct)
IF (ALLOCATED(diag%flux_diffuse)) &
  DEALLOCATE(diag%flux_diffuse)
! Microphysical diagnostics
IF (ALLOCATED(diag%re_strat)) &
  DEALLOCATE(diag%re_strat)
IF (ALLOCATED(diag%wgt_strat)) &
  DEALLOCATE(diag%wgt_strat)
IF (ALLOCATED(diag%lwp_strat)) &
  DEALLOCATE(diag%lwp_strat)
IF (ALLOCATED(diag%re_conv)) &
  DEALLOCATE(diag%re_conv)
IF (ALLOCATED(diag%wgt_conv)) &
  DEALLOCATE(diag%wgt_conv)
IF (ALLOCATED(diag%ntot_diag)) &
  DEALLOCATE(diag%ntot_diag)
IF (ALLOCATED(diag%strat_lwc_diag)) &
  DEALLOCATE(diag%strat_lwc_diag)
IF (ALLOCATED(diag%so4_ccn_diag)) &
  DEALLOCATE(diag%so4_ccn_diag)
IF (ALLOCATED(diag%cond_samp_wgt)) &
  DEALLOCATE(diag%cond_samp_wgt)
IF (ALLOCATED(diag%weighted_re)) &
  DEALLOCATE(diag%weighted_re)
IF (ALLOCATED(diag%sum_weight_re)) &
  DEALLOCATE(diag%sum_weight_re)
IF (ALLOCATED(diag%weighted_warm_re)) &
  DEALLOCATE(diag%weighted_warm_re)
IF (ALLOCATED(diag%sum_weight_warm_re)) &
  DEALLOCATE(diag%sum_weight_warm_re)
IF (ALLOCATED(diag%cdnc_ct_diag)) &
  DEALLOCATE(diag%cdnc_ct_diag)
IF (ALLOCATED(diag%cdnc_ct_weight)) &
  DEALLOCATE(diag%cdnc_ct_weight)
IF (ALLOCATED(diag%nc_diag)) &
  DEALLOCATE(diag%nc_diag)
IF (ALLOCATED(diag%nc_weight)) &
  DEALLOCATE(diag%nc_weight)
IF (ALLOCATED(diag%FlxSolBelow690nmSurf)) &
  DEALLOCATE(diag%FlxSolBelow690nmSurf)
IF (ALLOCATED(diag%FlxSeaBelow690nmSurf)) &
  DEALLOCATE(diag%FlxSeaBelow690nmSurf)
IF (ALLOCATED(diag%orog_corr)) &
  DEALLOCATE(diag%orog_corr)

IF (ALLOCATED(diag%cloud_extinction)) &
  DEALLOCATE(diag%cloud_extinction)
IF (ALLOCATED(diag%cloud_weight_extinction)) &
  DEALLOCATE(diag%cloud_weight_extinction)
IF (ALLOCATED(diag%ls_cloud_extinction)) &
  DEALLOCATE(diag%ls_cloud_extinction)
IF (ALLOCATED(diag%ls_cloud_weight_extinction)) &
  DEALLOCATE(diag%ls_cloud_weight_extinction)
IF (ALLOCATED(diag%cnv_cloud_extinction)) &
  DEALLOCATE(diag%cnv_cloud_extinction)
IF (ALLOCATED(diag%cnv_cloud_weight_extinction)) &
  DEALLOCATE(diag%cnv_cloud_weight_extinction)


! LW
IF (ALLOCATED(diag%total_cloud_cover)) &
  DEALLOCATE(diag%total_cloud_cover)
IF (ALLOCATED(diag%clear_olr)) &
  DEALLOCATE(diag%clear_olr)
IF (ALLOCATED(diag%down_flux_trop)) &
  DEALLOCATE(diag%down_flux_trop)
IF (ALLOCATED(diag%total_cloud_on_levels)) &
  DEALLOCATE(diag%total_cloud_on_levels)

IF (ALLOCATED(diag%cloud_absorptivity)) &
  DEALLOCATE(diag%cloud_absorptivity)
IF (ALLOCATED(diag%cloud_weight_absorptivity)) &
  DEALLOCATE(diag%cloud_weight_absorptivity)
IF (ALLOCATED(diag%ls_cloud_absorptivity)) &
  DEALLOCATE(diag%ls_cloud_absorptivity)
IF (ALLOCATED(diag%ls_cloud_weight_absorptivity)) &
  DEALLOCATE(diag%ls_cloud_weight_absorptivity)
IF (ALLOCATED(diag%cnv_cloud_absorptivity)) &
  DEALLOCATE(diag%cnv_cloud_absorptivity)
IF (ALLOCATED(diag%cnv_cloud_weight_absorptivity)) &
  DEALLOCATE(diag%cnv_cloud_weight_absorptivity)

! Deallocate the UKCA 3D aerosol-radiation diagnostics

IF (ALLOCATED(diag%ukca_aerosol_gsca)) &
  DEALLOCATE(diag%ukca_aerosol_gsca)
IF (ALLOCATED(diag%ukca_aerosol_sca)) &
  DEALLOCATE(diag%ukca_aerosol_sca)
IF (ALLOCATED(diag%ukca_aerosol_abs)) &
  DEALLOCATE(diag%ukca_aerosol_abs)
IF (ALLOCATED(diag%ukca_aerosol_ext)) &
  DEALLOCATE(diag%ukca_aerosol_ext)

! Deallocate the CLASSIC 3D aerosol-radiation diagnostics

IF (ALLOCATED(diag%clas_aerosol_sca)) &
  DEALLOCATE(diag%clas_aerosol_sca)
IF (ALLOCATED(diag%clas_aerosol_abs)) &
  DEALLOCATE(diag%clas_aerosol_abs)
IF (ALLOCATED(diag%clas_aerosol_ext)) &
  DEALLOCATE(diag%clas_aerosol_ext)

! Deallocate the UKCA aerosol optical depth diagnostics

IF (ALLOCATED(diag%aaod_ukca_cor_ins)) &
  DEALLOCATE(diag%aaod_ukca_cor_ins)
IF (ALLOCATED(diag%aaod_ukca_acc_ins)) &
  DEALLOCATE(diag%aaod_ukca_acc_ins)
IF (ALLOCATED(diag%aaod_ukca_ait_ins)) &
  DEALLOCATE(diag%aaod_ukca_ait_ins)
IF (ALLOCATED(diag%aaod_ukca_cor_sol)) &
  DEALLOCATE(diag%aaod_ukca_cor_sol)
IF (ALLOCATED(diag%aaod_ukca_acc_sol)) &
  DEALLOCATE(diag%aaod_ukca_acc_sol)
IF (ALLOCATED(diag%aaod_ukca_ait_sol)) &
  DEALLOCATE(diag%aaod_ukca_ait_sol)

! Deallocate the UKCA stratospheric aerosol optical depth diagnostics

IF (ALLOCATED(diag%sod_ukca_cor_ins)) &
  DEALLOCATE(diag%sod_ukca_cor_ins)
IF (ALLOCATED(diag%sod_ukca_acc_ins)) &
  DEALLOCATE(diag%sod_ukca_acc_ins)
IF (ALLOCATED(diag%sod_ukca_ait_ins)) &
  DEALLOCATE(diag%sod_ukca_ait_ins)
IF (ALLOCATED(diag%sod_ukca_cor_sol)) &
  DEALLOCATE(diag%sod_ukca_cor_sol)
IF (ALLOCATED(diag%sod_ukca_acc_sol)) &
  DEALLOCATE(diag%sod_ukca_acc_sol)
IF (ALLOCATED(diag%sod_ukca_ait_sol)) &
  DEALLOCATE(diag%sod_ukca_ait_sol)

! Deallocate the UKCA aerosol optical depth diagnostics

IF (ALLOCATED(diag%aod_ukca_cor_ins)) &
  DEALLOCATE(diag%aod_ukca_cor_ins)
IF (ALLOCATED(diag%aod_ukca_acc_ins)) &
  DEALLOCATE(diag%aod_ukca_acc_ins)
IF (ALLOCATED(diag%aod_ukca_ait_ins)) &
  DEALLOCATE(diag%aod_ukca_ait_ins)
IF (ALLOCATED(diag%aod_ukca_cor_sol)) &
  DEALLOCATE(diag%aod_ukca_cor_sol)
IF (ALLOCATED(diag%aod_ukca_acc_sol)) &
  DEALLOCATE(diag%aod_ukca_acc_sol)
IF (ALLOCATED(diag%aod_ukca_ait_sol)) &
  DEALLOCATE(diag%aod_ukca_ait_sol)

! Deallocate the CLASSIC absorption optical depth diagnostics

IF (ALLOCATED(diag%aaod_nitrate)) &
  DEALLOCATE(diag%aaod_nitrate)
IF (ALLOCATED(diag%aaod_ocff)) &
  DEALLOCATE(diag%aaod_ocff)
IF (ALLOCATED(diag%aaod_biogenic)) &
  DEALLOCATE(diag%aaod_biogenic)
IF (ALLOCATED(diag%aaod_biomass)) &
  DEALLOCATE(diag%aaod_biomass)
IF (ALLOCATED(diag%aaod_soot)) &
  DEALLOCATE(diag%aaod_soot)
IF (ALLOCATED(diag%aaod_seasalt)) &
  DEALLOCATE(diag%aaod_seasalt)
IF (ALLOCATED(diag%aaod_dust)) &
  DEALLOCATE(diag%aaod_dust)
IF (ALLOCATED(diag%aaod_sulphate)) &
  DEALLOCATE(diag%aaod_sulphate)

! Deallocate the aerosol optical depth diagnostics

IF (ALLOCATED(diag%aod_prog_nitrate)) &
  DEALLOCATE(diag%aod_prog_nitrate)
IF (ALLOCATED(diag%aod_prog_ocff)) &
  DEALLOCATE(diag%aod_prog_ocff)
IF (ALLOCATED(diag%aod_prog_biomass)) &
  DEALLOCATE(diag%aod_prog_biomass)
IF (ALLOCATED(diag%aod_prog_soot)) &
  DEALLOCATE(diag%aod_prog_soot)
IF (ALLOCATED(diag%aod_prog_seasalt)) &
  DEALLOCATE(diag%aod_prog_seasalt)
IF (ALLOCATED(diag%aod_prog_dust)) &
  DEALLOCATE(diag%aod_prog_dust)
IF (ALLOCATED(diag%aod_prog_sulphate)) &
  DEALLOCATE(diag%aod_prog_sulphate)
IF (ALLOCATED(diag%angst_total_radn)) &
  DEALLOCATE(diag%angst_total_radn)
IF (ALLOCATED(diag%aod_total_radn)) &
  DEALLOCATE(diag%aod_total_radn)
IF (ALLOCATED(diag%aod_nitrate)) &
  DEALLOCATE(diag%aod_nitrate)
IF (ALLOCATED(diag%aod_delta)) &
  DEALLOCATE(diag%aod_delta)
IF (ALLOCATED(diag%aod_ocff)) &
  DEALLOCATE(diag%aod_ocff)
IF (ALLOCATED(diag%aod_biogenic)) &
  DEALLOCATE(diag%aod_biogenic)
IF (ALLOCATED(diag%aod_biomass)) &
  DEALLOCATE(diag%aod_biomass)
IF (ALLOCATED(diag%aod_soot)) &
  DEALLOCATE(diag%aod_soot)
IF (ALLOCATED(diag%aod_seasalt)) &
  DEALLOCATE(diag%aod_seasalt)
IF (ALLOCATED(diag%aod_dust)) &
  DEALLOCATE(diag%aod_dust)
IF (ALLOCATED(diag%aod_sulphate)) &
  DEALLOCATE(diag%aod_sulphate)

IF (ALLOCATED(diag%l_ukca_aerosol_gsca)) &
  DEALLOCATE(diag%l_ukca_aerosol_gsca)
IF (ALLOCATED(diag%l_ukca_aerosol_sca)) &
  DEALLOCATE(diag%l_ukca_aerosol_sca)
IF (ALLOCATED(diag%l_ukca_aerosol_abs)) &
  DEALLOCATE(diag%l_ukca_aerosol_abs)
IF (ALLOCATED(diag%l_ukca_aerosol_ext)) &
  DEALLOCATE(diag%l_ukca_aerosol_ext)
IF (ALLOCATED(diag%l_clas_aerosol_sca)) &
  DEALLOCATE(diag%l_clas_aerosol_sca)
IF (ALLOCATED(diag%l_clas_aerosol_abs)) &
  DEALLOCATE(diag%l_clas_aerosol_abs)
IF (ALLOCATED(diag%l_clas_aerosol_ext)) &
  DEALLOCATE(diag%l_clas_aerosol_ext)

! Deallocate the grid-box mean cloud diagnostics

IF (ALLOCATED(diag%ls_qcl_rad)) &
  DEALLOCATE(diag%ls_qcl_rad)
IF (ALLOCATED(diag%ls_qcf_rad)) &
  DEALLOCATE(diag%ls_qcf_rad)
IF (ALLOCATED(diag%cc_qcl_rad)) &
  DEALLOCATE(diag%cc_qcl_rad)
IF (ALLOCATED(diag%cc_qcf_rad)) &
  DEALLOCATE(diag%cc_qcf_rad)
IF (ALLOCATED(diag%ls_cl_rad)) &
  DEALLOCATE(diag%ls_cl_rad)
IF (ALLOCATED(diag%ls_cf_rad)) &
  DEALLOCATE(diag%ls_cf_rad)
IF (ALLOCATED(diag%cc_cl_rad)) &
  DEALLOCATE(diag%cc_cl_rad)
IF (ALLOCATED(diag%cc_cf_rad)) &
  DEALLOCATE(diag%cc_cf_rad)

! Deallocate cloud water paths in radiation
IF (ALLOCATED(diag%ls_qcl_rad_path)) &
  DEALLOCATE(diag%ls_qcl_rad_path)
IF (ALLOCATED(diag%ls_qcf_rad_path)) &
  DEALLOCATE(diag%ls_qcf_rad_path)
IF (ALLOCATED(diag%cc_qcl_rad_path)) &
  DEALLOCATE(diag%cc_qcl_rad_path)
IF (ALLOCATED(diag%cc_qcf_rad_path)) &
  DEALLOCATE(diag%cc_qcf_rad_path)
IF (ALLOCATED(diag%ccore_qcl_rad_path)) &
  DEALLOCATE(diag%ccore_qcl_rad_path)
IF (ALLOCATED(diag%ccore_qcf_rad_path)) &
  DEALLOCATE(diag%ccore_qcf_rad_path)
! =========================================================

! Deallocate the grid-box mean convective core cloud diagnostics

IF (ALLOCATED(diag%ccore_clt_rad)) &
  DEALLOCATE(diag%ccore_clt_rad)
IF (ALLOCATED(diag%ccore_qcl_rad)) &
  DEALLOCATE(diag%ccore_qcl_rad)
IF (ALLOCATED(diag%ccore_qcf_rad)) &
  DEALLOCATE(diag%ccore_qcf_rad)

! Deallocate the weighted effective dimension diagnostics

IF (ALLOCATED(diag%ls_del_rad)) &
  DEALLOCATE(diag%ls_del_rad)
IF (ALLOCATED(diag%ls_def_rad)) &
  DEALLOCATE(diag%ls_def_rad)
IF (ALLOCATED(diag%cc_del_rad)) &
  DEALLOCATE(diag%cc_del_rad)
IF (ALLOCATED(diag%cc_def_rad)) &
  DEALLOCATE(diag%cc_def_rad)

END SUBROUTINE deallocate_diag
!------------------------------------------------------------------------------

END MODULE def_diag
