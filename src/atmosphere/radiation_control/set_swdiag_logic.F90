! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Set the SW diagnostic flags to true if necessary

! Method:

! If a radiation diagnostic has been chosen in STASH then the
! flag of the corresponding radiation diagnostic in the structure
! SW_diag is set to true.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

!-----------------------------------------------------------------------
MODULE set_swdiag_logic_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_SWDIAG_LOGIC_MOD'
CONTAINS

SUBROUTINE set_swdiag_logic(sf, nitems, nsects, j_sw, i_off)

USE sw_diag_mod, ONLY: sw_diag
USE rad_input_mod, ONLY: l_rad_perturb, l_radiance, l_forcing
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

USE model_domain_mod, ONLY: model_type, mt_single_column

IMPLICIT NONE

! arguments with intent in

INTEGER, INTENT(IN) :: nitems    ! item number
INTEGER, INTENT(IN) :: nsects    ! section number
INTEGER, INTENT(IN) :: j_sw      ! call to SW radiation
INTEGER, INTENT(IN) :: i_off     ! offset for diagnostics

LOGICAL, INTENT(IN) :: sf(0:,0:)
!        STASH Flags

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SET_SWDIAG_LOGIC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


SELECT CASE (model_type)

CASE DEFAULT

  IF (l_radiance .AND. (j_sw > 1)) THEN

    ! TOA radiances
    sw_diag(j_sw)%l_toa_radiance      = sf(297+i_off,1)

  ELSE IF (l_forcing .AND. (j_sw > 1)) THEN

    ! Fluxes and Heating Rates with an applied forcing
    sw_diag(j_sw)%l_flux_up                           = sf(217+i_off,1)
    sw_diag(j_sw)%l_flux_down                         = sf(218+i_off,1)
    sw_diag(j_sw)%l_flux_up_clear                     = sf(219+i_off,1)
    sw_diag(j_sw)%l_flux_down_clear                   = sf(220+i_off,1)
    sw_diag(j_sw)%l_clear_hr                          = sf(233+i_off,1)
    sw_diag(j_sw)%l_net_flux_trop                     = sf(237+i_off,1)
    sw_diag(j_sw)%l_up_flux_trop                      = sf(238+i_off,1)
    sw_diag(j_sw)%l_flux_direct_sph                   = sf(272+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_sph             = sf(273+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clean_sph             = sf(274+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_clean_sph       = sf(275+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_div                   = sf(276+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_div             = sf(277+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clean_div             = sf(278+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_clean_div       = sf(279+i_off,1)      
    sw_diag(j_sw)%l_flux_up_band                      = sf(509+i_off,1)
    sw_diag(j_sw)%l_flux_down_band                    = sf(510+i_off,1)
    sw_diag(j_sw)%l_flux_up_clear_band                = sf(511+i_off,1)
    sw_diag(j_sw)%l_flux_down_clear_band              = sf(512+i_off,1)      
    sw_diag(j_sw)%l_emission_spectrum                 = sf(513+i_off,1)      
    sw_diag(j_sw)%l_emission_spectrum_clear           = sf(514+i_off,1)      
    sw_diag(j_sw)%l_emission_spectrum_clean           = sf(515+i_off,1)      
    sw_diag(j_sw)%l_emission_spectrum_clear_clean     = sf(516+i_off,1)      
    sw_diag(j_sw)%l_flux_up_clean                     = sf(517+i_off,1)
    sw_diag(j_sw)%l_flux_down_clean                   = sf(518+i_off,1)
    sw_diag(j_sw)%l_flux_up_clear_clean               = sf(519+i_off,1)
    sw_diag(j_sw)%l_flux_down_clear_clean             = sf(520+i_off,1)
    sw_diag(j_sw)%l_flux_direct_sph_band              = sf(537+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_sph_band        = sf(538+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clean_sph_band        = sf(539+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_clean_sph_band  = sf(540+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_div_band              = sf(541+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_div_band        = sf(542+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clean_div_band        = sf(543+i_off,1)      
    sw_diag(j_sw)%l_flux_direct_clear_clean_div_band  = sf(544+i_off,1)      
    sw_diag(j_sw)%l_flux_up_clean_band                = sf(545+i_off,1)
    sw_diag(j_sw)%l_flux_down_clean_band              = sf(546+i_off,1)
    sw_diag(j_sw)%l_flux_up_clear_clean_band          = sf(547+i_off,1)
    sw_diag(j_sw)%l_flux_down_clear_clean_band        = sf(548+i_off,1)
    sw_diag(j_sw)%l_transmission_spectrum             = sf(555+i_off,1)      
    sw_diag(j_sw)%l_transmission_spectrum_clear       = sf(556+i_off,1)      
    sw_diag(j_sw)%l_transmission_spectrum_clean       = sf(557+i_off,1)      
    sw_diag(j_sw)%l_transmission_spectrum_clear_clean = sf(558+i_off,1)      

    sw_diag(j_sw)%l_diag_call = &
      sw_diag(j_sw)%l_flux_up                           .OR. & 
      sw_diag(j_sw)%l_flux_down                         .OR. & 
      sw_diag(j_sw)%l_flux_up_clear                     .OR. & 
      sw_diag(j_sw)%l_flux_down_clear                   .OR. & 
      sw_diag(j_sw)%l_clear_hr                          .OR. & 
      sw_diag(j_sw)%l_net_flux_trop                     .OR. & 
      sw_diag(j_sw)%l_up_flux_trop                      .OR. & 
      sw_diag(j_sw)%l_flux_direct_sph                   .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_sph             .OR. & 
      sw_diag(j_sw)%l_flux_direct_clean_sph             .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_clean_sph       .OR. & 
      sw_diag(j_sw)%l_flux_direct_div                   .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_div             .OR. & 
      sw_diag(j_sw)%l_flux_direct_clean_div             .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_clean_div       .OR. & 
      sw_diag(j_sw)%l_flux_up_band                      .OR. & 
      sw_diag(j_sw)%l_flux_down_band                    .OR. & 
      sw_diag(j_sw)%l_flux_up_clear_band                .OR. & 
      sw_diag(j_sw)%l_flux_down_clear_band              .OR. & 
      sw_diag(j_sw)%l_emission_spectrum                 .OR. & 
      sw_diag(j_sw)%l_emission_spectrum_clear           .OR. & 
      sw_diag(j_sw)%l_emission_spectrum_clean           .OR. & 
      sw_diag(j_sw)%l_emission_spectrum_clear_clean     .OR. & 
      sw_diag(j_sw)%l_flux_up_clean                     .OR. & 
      sw_diag(j_sw)%l_flux_down_clean                   .OR. & 
      sw_diag(j_sw)%l_flux_up_clear_clean               .OR. & 
      sw_diag(j_sw)%l_flux_down_clear_clean             .OR. & 
      sw_diag(j_sw)%l_flux_direct_sph_band              .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_sph_band        .OR. & 
      sw_diag(j_sw)%l_flux_direct_clean_sph_band        .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_clean_sph_band  .OR. & 
      sw_diag(j_sw)%l_flux_direct_div_band              .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_div_band        .OR. & 
      sw_diag(j_sw)%l_flux_direct_clean_div_band        .OR. & 
      sw_diag(j_sw)%l_flux_direct_clear_clean_div_band  .OR. & 
      sw_diag(j_sw)%l_flux_up_clean_band                .OR. & 
      sw_diag(j_sw)%l_flux_down_clean_band              .OR. & 
      sw_diag(j_sw)%l_flux_up_clear_clean_band          .OR. & 
      sw_diag(j_sw)%l_flux_down_clear_clean_band        .OR. & 
      sw_diag(j_sw)%l_transmission_spectrum             .OR. & 
      sw_diag(j_sw)%l_transmission_spectrum_clear       .OR. & 
      sw_diag(j_sw)%l_transmission_spectrum_clean       .OR. & 
      sw_diag(j_sw)%l_transmission_spectrum_clear_clean

  ELSE

    ! Fluxes and Heating Rates
    sw_diag(j_sw)%l_flux_up           = sf(217,1)
    sw_diag(j_sw)%l_flux_down         = sf(218,1)
    sw_diag(j_sw)%l_solar_out_toa     = sf(208,1)
    sw_diag(j_sw)%l_surface_down_flux = sf(235,1)
    sw_diag(j_sw)%l_net_flux_trop     = sf(237,1)
    sw_diag(j_sw)%l_up_flux_trop      = sf(238,1)

    ! Diffuse & Direct Flux
    sw_diag(j_sw)%l_flux_direct       = sf(230,1)
    sw_diag(j_sw)%l_flux_diffuse      = sf(231,1)

    ! Direct fluxes for spherical geometry
    sw_diag(j_sw)%l_flux_direct_sph   = sf(272,1)      
    sw_diag(j_sw)%l_flux_direct_div   = sf(276,1)      

    ! Diagnostics for photosynthetically active surface radiation
    sw_diag(j_sw)%l_FlxSolBelow690nmSurf = sf(259,1)
    sw_diag(j_sw)%l_FlxSeaBelow690nmSurf = sf(260,1)

    ! Diagnostics for orography correction
    sw_diag(j_sw)%l_orog_corr         = sf(295,1)

    ! Diagnostics for sea salt
    sw_diag(j_sw)%seasalt_film_flag   = sf(247,1)
    sw_diag(j_sw)%seasalt_jet_flag    = sf(248,1)

    ! Albedo scaling to obs diagnostics
    sw_diag(j_sw)%l_vis_albedo_sc     = sf(270,1)
    sw_diag(j_sw)%l_nir_albedo_sc     = sf(271,1)

    IF (.NOT. (l_rad_perturb .AND. (j_sw == 2))) THEN
      ! For the incremental time-stepping scheme (l_rad_perturb) many
      ! of the diagnostics are not calculated on the "cloud only"
      ! radiation calls (j_sw==2).

      ! Mask for radiation calculations (REAL 1.0 or 0.0)
      sw_diag(j_sw)%l_rad_mask          = sf(200,1)

      ! Clear-sky Fluxes and Heating Rates
      sw_diag(j_sw)%l_flux_up_clear     = sf(219,1)
      sw_diag(j_sw)%l_flux_down_clear   = sf(220,1)
      sw_diag(j_sw)%l_toa_clear_weighted= sf(228,1)
      sw_diag(j_sw)%l_total_clear_area  = sf(229,1)
      sw_diag(j_sw)%l_solar_out_clear   = sf(209,1)
      sw_diag(j_sw)%l_surf_down_clr     = sf(210,1)
      sw_diag(j_sw)%l_surf_up_clr       = sf(211,1)
      sw_diag(j_sw)%l_clear_hr          = sf(233,1)

      ! Clear-sky direct fluxes for spherical geometry
      sw_diag(j_sw)%l_flux_direct_clear_sph = sf(273,1)      
      sw_diag(j_sw)%l_flux_direct_clear_div = sf(277,1)      

      ! UV Fluxes
      sw_diag(j_sw)%l_uvflux_direct     = sf(212,1)
      sw_diag(j_sw)%l_uvflux_up         = sf(213,1)
      sw_diag(j_sw)%l_uvflux_down       = sf(214,1)
      sw_diag(j_sw)%l_surf_uv           = sf(288,1)
      sw_diag(j_sw)%l_surf_uv_clr       = sf(289,1)
              
      ! Microphysical diagnostics
      sw_diag(j_sw)%re_strat_flag       = sf(221,1)
      sw_diag(j_sw)%wgt_strat_flag      = sf(223,1)
      sw_diag(j_sw)%lwp_strat_flag      = sf(224,1)
      sw_diag(j_sw)%re_conv_flag        = sf(225,1)
      sw_diag(j_sw)%wgt_conv_flag       = sf(226,1)
      sw_diag(j_sw)%ntot_diag_flag      = sf(241,1)
      sw_diag(j_sw)%strat_lwc_diag_flag = sf(242,1)
      sw_diag(j_sw)%so4_ccn_diag_flag   = sf(243,1)
      sw_diag(j_sw)%cond_samp_wgt_flag  = sf(244,1)
      sw_diag(j_sw)%weighted_re_flag    = sf(245,1)
      sw_diag(j_sw)%sum_weight_re_flag  = sf(246,1)
      sw_diag(j_sw)%wgtd_warm_re_flag   = sf(254,1)
      sw_diag(j_sw)%sum_wgt_warm_re_flag= sf(255,1)
      sw_diag(j_sw)%nc_diag_flag        = sf(280,1)
      sw_diag(j_sw)%nc_weight_flag      = sf(281,1)
      sw_diag(j_sw)%cdnc_ct_diag_flag   = sf(298,1)
      sw_diag(j_sw)%cdnc_ct_weight_flag = sf(299,1)

      ! Extinction and absorptivity diagnostics
      sw_diag(j_sw)%l_cloud_extinction           = sf(262,1)
      sw_diag(j_sw)%l_cloud_weight_extinction    = sf(263,1)
      sw_diag(j_sw)%l_ls_cloud_extinction        = sf(264,1)
      sw_diag(j_sw)%l_ls_cloud_weight_extinction = sf(265,1)
      sw_diag(j_sw)%l_cnv_cloud_extinction       = sf(266,1)
      sw_diag(j_sw)%l_cnv_cloud_weight_extinction= sf(267,1)

      ! Surface Albedo Diagnostics
      sw_diag(j_sw)%l_direct_albedo  = sf(268,1)
      sw_diag(j_sw)%l_diffuse_albedo = sf(269,1)

      ! Aerosol Optical Properties
      sw_diag(j_sw)%l_aerosol_optical_depth      = sf(506,1)
      sw_diag(j_sw)%l_aerosol_scat_optical_depth = sf(507,1)
      sw_diag(j_sw)%l_aerosol_asymmetry_scat     = sf(508,1)

      ! Band-by-band Fluxes
      sw_diag(j_sw)%l_flux_up_band               = sf(509,1)
      sw_diag(j_sw)%l_flux_down_band             = sf(510,1)
      sw_diag(j_sw)%l_flux_up_clear_band         = sf(511,1)
      sw_diag(j_sw)%l_flux_down_clear_band       = sf(512,1)      
      sw_diag(j_sw)%l_emission_spectrum          = sf(513,1)      
      sw_diag(j_sw)%l_emission_spectrum_clear    = sf(514,1)      
      sw_diag(j_sw)%l_flux_direct_sph_band       = sf(537,1)      
      sw_diag(j_sw)%l_flux_direct_clear_sph_band = sf(538,1)      
      sw_diag(j_sw)%l_flux_direct_div_band       = sf(541,1)      
      sw_diag(j_sw)%l_flux_direct_clear_div_band = sf(542,1)      
      sw_diag(j_sw)%l_transmission_spectrum      = sf(555,1)      
      sw_diag(j_sw)%l_transmission_spectrum_clear= sf(556,1)      

      ! Clean-air Fluxes
      sw_diag(j_sw)%l_flux_up_clean                     = sf(517,1)
      sw_diag(j_sw)%l_flux_down_clean                   = sf(518,1)
      sw_diag(j_sw)%l_flux_up_clear_clean               = sf(519,1)
      sw_diag(j_sw)%l_flux_down_clear_clean             = sf(520,1)
      sw_diag(j_sw)%l_flux_up_clean_band                = sf(545,1)
      sw_diag(j_sw)%l_flux_down_clean_band              = sf(546,1)
      sw_diag(j_sw)%l_flux_up_clear_clean_band          = sf(547,1)
      sw_diag(j_sw)%l_flux_down_clear_clean_band        = sf(548,1)
      sw_diag(j_sw)%l_emission_spectrum_clean           = sf(515,1)      
      sw_diag(j_sw)%l_emission_spectrum_clear_clean     = sf(516,1)      
      sw_diag(j_sw)%l_flux_direct_clean_sph             = sf(274,1)      
      sw_diag(j_sw)%l_flux_direct_clear_clean_sph       = sf(275,1)      
      sw_diag(j_sw)%l_flux_direct_clean_sph_band        = sf(539,1)      
      sw_diag(j_sw)%l_flux_direct_clear_clean_sph_band  = sf(540,1)      
      sw_diag(j_sw)%l_flux_direct_clean_div             = sf(278,1)      
      sw_diag(j_sw)%l_flux_direct_clear_clean_div       = sf(279,1)      
      sw_diag(j_sw)%l_flux_direct_clean_div_band        = sf(543,1)      
      sw_diag(j_sw)%l_flux_direct_clear_clean_div_band  = sf(544,1)      
      sw_diag(j_sw)%l_transmission_spectrum_clean       = sf(557,1)      
      sw_diag(j_sw)%l_transmission_spectrum_clear_clean = sf(558,1)      

      ! Fluxes after green-house gas forcing
      sw_diag(j_sw)%l_flux_up_forc                    = sf(521,1)
      sw_diag(j_sw)%l_flux_down_forc                  = sf(522,1)
      sw_diag(j_sw)%l_flux_up_clear_forc              = sf(523,1)
      sw_diag(j_sw)%l_flux_down_clear_forc            = sf(524,1)
      sw_diag(j_sw)%l_flux_up_forc_band               = sf(525,1)
      sw_diag(j_sw)%l_flux_down_forc_band             = sf(526,1)
      sw_diag(j_sw)%l_flux_up_clear_forc_band         = sf(527,1)
      sw_diag(j_sw)%l_flux_down_clear_forc_band       = sf(528,1)
      sw_diag(j_sw)%l_flux_direct_sph_forc            = sf(529,1)
      sw_diag(j_sw)%l_flux_direct_clear_sph_forc      = sf(530,1)
      sw_diag(j_sw)%l_flux_direct_div_forc            = sf(531,1)
      sw_diag(j_sw)%l_flux_direct_clear_div_forc      = sf(532,1)
      sw_diag(j_sw)%l_flux_direct_sph_forc_band       = sf(533,1)
      sw_diag(j_sw)%l_flux_direct_clear_sph_forc_band = sf(534,1)
      sw_diag(j_sw)%l_flux_direct_div_forc_band       = sf(535,1)
      sw_diag(j_sw)%l_flux_direct_clear_div_forc_band = sf(536,1)

      ! EasyAerosol
      sw_diag(j_sw)%l_easyaerosol_extinction     = sf(550,1)
      sw_diag(j_sw)%l_easyaerosol_absorption     = sf(551,1)
      sw_diag(j_sw)%l_easyaerosol_scattering     = sf(552,1)
      sw_diag(j_sw)%l_easyaerosol_asytimscat     = sf(553,1)

    END IF ! .not.(l_rad_perturb.and.(j_sw == 2))

  END IF

CASE (mt_single_column)

  ! Only the most common diagnostics are obtained by default.
  sw_diag(j_sw)%l_solar_out_toa      = .TRUE.
  sw_diag(j_sw)%l_surface_down_flux  = .TRUE.

  IF (.NOT. (l_rad_perturb .AND. (j_sw == 2))) THEN
    ! For the incremental time-stepping scheme (l_rad_perturb) many
    ! of the diagnostics are not calculated on the "cloud only"
    ! radiation calls (j_sw==2).
    sw_diag(j_sw)%l_solar_out_clear  = .TRUE.
    sw_diag(j_sw)%l_surf_down_clr    = .TRUE.
    sw_diag(j_sw)%l_surf_up_clr      = .TRUE.
    sw_diag(j_sw)%l_clear_hr         = .TRUE.
    sw_diag(j_sw)%re_strat_flag      = .TRUE.
    sw_diag(j_sw)%wgt_strat_flag     = .TRUE.
    sw_diag(j_sw)%lwp_strat_flag     = .TRUE.
    sw_diag(j_sw)%re_conv_flag       = .TRUE.
    sw_diag(j_sw)%wgt_conv_flag      = .TRUE.
  END IF

END SELECT

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_swdiag_logic
END MODULE set_swdiag_logic_mod
