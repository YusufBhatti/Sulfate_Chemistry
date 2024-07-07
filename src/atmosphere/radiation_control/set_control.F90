! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_control_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_CONTROL_MOD'
CONTAINS

SUBROUTINE set_control(                                                        &

! Structures for the core radiation code interface
  control, spectrum,                                                           &

! Control flags
  l_climat_aerosol, l_use_sulpc_direct, l_use_soot_direct, l_use_biogenic,     &
  l_use_dust, l_use_bmass_direct, l_use_ocff_direct, l_use_nitrate_direct,     &
  l_use_seasalt_direct, l_murk_rad, l_use_sulpc_indirect, n_arcl_species,      &
  l_use_ukca_radaer, l_use_glomap_clim_radaer,                                 &
  l_easyaerosol,                                                               &

! Diagnostic options
  diag, l_cosp)

! Subroutine to set algorithmic options for the core radiation code
!
! Purpose:
!   Algorithmic options and array sizes to be set interactively
!   are determined.
!

USE rad_pcf
USE def_control,  ONLY: StrCtrl, allocate_control
USE def_spectrum, ONLY: StrSpecData
USE def_diag,     ONLY: StrDiag
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE ereport_mod,  ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE planet_constants_mod, ONLY: l_planet_grey_surface, l_planet_aerosol
USE solinc_data,  ONLY: l_orog
USE gen_phys_inputs_mod, ONLY: l_mr_physics

USE r2_set_uv_weight_mod, ONLY: r2_set_uv_weight
IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(INOUT) :: control

! Spectral data:
TYPE (StrSpecData), INTENT(IN)    :: spectrum

! Aerosol Flags:
LOGICAL, INTENT(IN) ::                                                         &
  l_climat_aerosol,                                                            &
!   Flag for climatological aerosols
  l_use_sulpc_direct,                                                          &
!   Flag to include the direct effect of sulphate aerosols
  l_use_soot_direct,                                                           &
!   Flag to include the direct effect of soot aerosols
  l_use_sulpc_indirect,                                                        &
!   Flag to include the indirect effect of sulphate aerosols
  l_use_seasalt_direct,                                                        &
!   Flag to include the direct effect of seasalt aerosols
  l_use_biogenic,                                                              &
!   Flag to include the direct effect of biogenic aerosols
  l_use_dust,                                                                  &
!   Flag to use direct effect of mineral dust
  l_use_bmass_direct,                                                          &
!   Flag to use direct radiative effect of biomass smoke
  l_use_ocff_direct,                                                           &
!   Flag to use direct radiative effect of fossil-fuel oc
  l_use_nitrate_direct,                                                        &
!   Flag to use direct radiative effect of nitrate aerosols
  l_murk_rad,                                                                  &
!   Flag to include urban aerosols
  l_use_ukca_radaer, l_use_glomap_clim_radaer,                                 &
!   Include radiative effect of GLOMAP MODE aerosols
  l_easyaerosol
!   Include radiative effect of EasyAerosol climatologies

INTEGER, INTENT(IN) :: n_arcl_species
!   Number of species from the NWP aerosol climatology.
!   Zero if the NWP climatology is not used.

! Diagnostic options
TYPE (StrDiag), INTENT(IN) :: diag

LOGICAL, INTENT(IN) :: l_cosp
!   Some diagnostics are required for the COSP simulator


! Local variables.
INTEGER :: i
!   Loop variable

INTEGER                      :: ierr = i_normal
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SET_CONTROL'
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set the last band to use as the last band in the spectral file
control%last_band = spectrum%basic%n_band

! Control flag to indicate mixing ratios are with respect to dry mass
control%l_mixing_ratio = l_mr_physics

! Allocate band-by-band control options
CALL allocate_control(control, spectrum)


! Set diagnostic flags
control%l_flux_up_band               = diag%l_flux_up_band .OR.                &
                                       diag%l_emission_spectrum .OR.           &
                                       diag%l_transmission_spectrum
control%l_flux_down_band             = diag%l_flux_down_band
control%l_flux_up_clear_band         = diag%l_flux_up_clear_band .OR.          &
                                       diag%l_emission_spectrum_clear .OR.     &
                                       diag%l_transmission_spectrum_clear
control%l_flux_down_clear_band       = diag%l_flux_down_clear_band
control%l_flux_direct_sph_band       = diag%l_flux_direct_sph_band .OR.        &
                                       diag%l_transmission_spectrum
control%l_flux_direct_clear_sph_band = diag%l_flux_direct_clear_sph_band .OR.  &
                                       diag%l_transmission_spectrum_clear
control%l_flux_direct_div_band       = diag%l_flux_direct_div_band
control%l_flux_direct_clear_div_band = diag%l_flux_direct_clear_div_band
control%l_aerosol_absorption_band    = diag%l_aerosol_optical_depth
control%l_aerosol_scattering_band    = diag%l_aerosol_optical_depth .OR.       &
                                       diag%l_aerosol_scat_optical_depth
control%l_aerosol_asymmetry_band     = diag%l_aerosol_asymmetry_scat

control%l_clear = &
  diag%l_flux_up_clear         .OR. control%l_flux_up_clear_band         .OR.  &
  diag%l_flux_down_clear       .OR. control%l_flux_down_clear_band       .OR.  &
  diag%l_flux_direct_clear_sph .OR. control%l_flux_direct_clear_sph_band .OR.  &
  diag%l_flux_direct_clear_div .OR. control%l_flux_direct_clear_div_band

SELECT CASE (control%isolir)

CASE (ip_solar)
  control%l_clear = control%l_clear                                            &
    .OR. diag%l_solar_out_clear                                                &
    .OR. diag%l_toa_clear_weighted                                             &
    .OR. diag%l_surf_down_clr                                                  &
    .OR. diag%l_surf_up_clr                                                    &
    .OR. diag%l_clear_hr                                                       &
    .OR. (diag%l_cloud_extinction .AND.                                        &
         diag%l_cloud_weight_extinction)                                       &
    .OR. (diag%l_ls_cloud_extinction .AND.                                     &
         diag%l_ls_cloud_weight_extinction)                                    &
    .OR. (diag%l_cnv_cloud_extinction .AND.                                    &
         diag%l_cnv_cloud_weight_extinction) .OR. l_cosp
  control%l_cloud_extinction     = diag%l_cloud_extinction
  control%l_ls_cloud_extinction  = diag%l_ls_cloud_extinction .OR. l_cosp
  control%l_cnv_cloud_extinction = diag%l_cnv_cloud_extinction .OR. l_cosp

  IF (diag%l_uvflux_direct .OR.                                                &
      diag%l_uvflux_down   .OR.                                                &
      diag%l_uvflux_up     .OR.                                                &
      diag%l_surf_uv       .OR.                                                &
      diag%l_surf_uv_clr) THEN
    CALL r2_set_uv_weight(spectrum%basic%n_band,                               &
      spectrum%basic%l_present,                                                &
      spectrum%basic%n_band_exclude,                                           &
      spectrum%basic%index_exclude,                                            &
      spectrum%basic%wavelength_short,                                         &
      spectrum%basic%wavelength_long,                                          &
      control%weight_diag,                                                     &
      spectrum%dim%nd_band, spectrum%dim%nd_exclude,                           &
      spectrum%dim%nd_type)
  END IF
  control%l_flux_direct_band = control%l_flux_direct_band                      &
    .OR. diag%l_uvflux_direct
  control%l_flux_up_band = control%l_flux_up_band                              &
    .OR. diag%l_uvflux_up
  control%l_flux_down_band = control%l_flux_down_band                          &
    .OR. diag%l_uvflux_down                                                    &
    .OR. diag%l_surf_uv
  control%l_flux_down_clear_band = control%l_flux_down_clear_band              &
    .OR. diag%l_surf_uv_clr
  IF (control%l_spherical_solar) THEN
    control%l_flux_direct_sph_band = control%l_flux_direct_sph_band            &
      .OR. diag%l_surf_uv
    control%l_flux_direct_clear_sph_band = control%l_flux_direct_clear_sph_band&
      .OR. diag%l_surf_uv_clr
  END IF

CASE (ip_infra_red)
  control%l_clear = control%l_clear                                            &
    .OR. diag%l_clear_olr                                                      &
    .OR. diag%l_toa_clear_weighted                                             &
    .OR. diag%l_surf_down_clr                                                  &
    .OR. diag%l_clear_hr                                                       &
    .OR. (diag%l_cloud_absorptivity .AND.                                      &
         diag%l_cloud_weight_absorptivity)                                     &
    .OR. (diag%l_ls_cloud_absorptivity .AND.                                   &
         diag%l_ls_cloud_weight_absorptivity)                                  &
    .OR. (diag%l_cnv_cloud_absorptivity .AND.                                  &
         diag%l_cnv_cloud_weight_absorptivity) .OR. l_cosp
  control%l_cloud_absorptivity     = diag%l_cloud_absorptivity
  control%l_ls_cloud_absorptivity  = diag%l_ls_cloud_absorptivity.OR.l_cosp
  control%l_cnv_cloud_absorptivity = diag%l_cnv_cloud_absorptivity.OR.l_cosp

END SELECT

! Control flag for corrections to the direct solar flux at the surface
! for sloping terrain
control%l_orog = l_orog

! Decide on the final options for aerosols:
IF (l_use_glomap_clim_radaer) THEN
  control%l_aerosol_mode = l_use_glomap_clim_radaer
ELSE
  control%l_aerosol_mode = l_use_ukca_radaer .OR. l_easyaerosol
END IF
control%l_aerosol = l_use_sulpc_direct .OR.                                    &
                    l_use_dust .OR.                                            &
                    l_use_soot_direct .OR.                                     &
                    l_use_bmass_direct .OR.                                    &
                    l_use_ocff_direct .OR.                                     &
                    l_use_nitrate_direct .OR.                                  &
                    l_use_seasalt_direct .OR.                                  &
                    l_use_biogenic .OR.                                        &
                    l_murk_rad .OR.                                            &
                    l_climat_aerosol .OR.                                      &
                    l_planet_aerosol .OR.                                      &
                    (n_arcl_species > 0)

! Whilst l_aerosol_ccn is a generic flag for determining CCN from aerosol,
! the view is currently taken that sulphate aerosols must be included with
! all indirect effects, other aerosols being additional, so l_aerosol_ccn
! is assigned solely from l_use_sulpc_indirect.
control%l_aerosol_ccn=l_use_sulpc_indirect


! Set properties for individual bands.
DO i = 1, spectrum%basic%n_band
  control%map_channel(i)           = 1
  control%weight_band(i)           = 1.0
  control%i_scatter_method_band(i) = control%i_scatter_method
  control%i_gas_overlap_band(i)    = control%i_gas_overlap
  IF (ANY(spectrum%gas%i_scale_fnc(i,:) == ip_scale_ses2)) THEN
    ! If SES2 scaling is used in this band then the overlap must also use SES2:
    control%i_gas_overlap_band(i)  = ip_overlap_mix_ses2
  END IF
  IF (diag%l_surf_uv_clr) THEN
    control%l_clear_band(i) = control%weight_diag(i) > 0.0
  END IF
END DO


IF (control%i_angular_integration == ip_two_stream) THEN

  IF (control%l_rescale) control%n_order_forward=2

  ! We permit tiling of sea-ice points only with the two-stream
  ! option at present. Tiling is only of use in separating
  ! different components of the fluxes at the surface and in
  ! particular is not relevent to the calculation of TOA radiances.
  ! Tiling is not needed if the surface is to be treated as grey.
  IF (l_planet_grey_surface) THEN
    control%l_tile=.FALSE.
  ELSE
    control%l_tile=.TRUE.
  END IF

ELSE IF (control%i_angular_integration == ip_spherical_harmonic) THEN

  IF (control%i_sph_mode == ip_sph_mode_flux) THEN
    ! Map all bands to a single output channel
    control%n_channel = 1
    control%map_channel(1:spectrum%basic%n_band) = 1
  ELSE
    ! Map each band to a separate output channel
    control%n_channel = spectrum%basic%n_band
    DO i = 1, spectrum%basic%n_band
      control%map_channel(i) = i
    END DO
  END IF

  IF (control%l_rescale) control%n_order_forward=control%ls_global_trunc+1

  ! As currently implemented, Euler's transformation is applied
  ! only in its most basic form, adding just half of the last
  ! term in an alternating series.
  IF (control%l_euler_trnf) THEN
    control%euler_factor=0.5
  ELSE
    control%euler_factor=1.0
  END IF

  ! Clear-sky fluxes are not available from the spherical harmonic
  ! code in the same call as cloudy fluxes yet. If required, they
  ! should be diagnosed by using a separate call to the code with
  ! clouds switched off.
  IF ( control%l_clear ) THEN
    cmessage = 'Clear-sky fluxes not directly available in harmonics'
    ierr=i_err_fatal
    GO TO 9999
  END IF

  ! We permit tiling of sea-ice points only with the two-stream
  ! option at present.
  control%l_tile=.FALSE.

END IF


9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_control
END MODULE set_control_mod
