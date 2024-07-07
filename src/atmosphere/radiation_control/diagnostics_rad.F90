! *****************************COPYRIGHT******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT******************************
!
! Output radiation diagnostics
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------

MODULE diagnostics_rad_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'DIAGNOSTICS_RAD_MOD'
CONTAINS

SUBROUTINE diagnostics_rad(                                                    &
  sect, j_rad, diag, spectrum,                                                 &
  row_length, rows, ozone_levels, cloud_levels, ntiles, at_extremity, i_off,   &
  ! Fields common to SW and LW radiation
  T_n, T_inc, q_n, qcl_n, cf_n, cfl_n,                                         &
  T_latest, q_latest, qcl_latest, cf_latest, cfl_latest,                       &
  T_incr_diagnostic, surfnet, netsea, obs_solid_angle,                         &
  down_sice_weighted_cat, down_sice_weighted,                                  &
  ! SW diagnostic fields
  itoasw, surfsw_cor, toasw_cor, surfdir_cor, surfdif_cor,                     &
  flux_below_690nm_surf, photosynth_act_rad, flxdirparsurf,                    &
  f_orog, sw_net_land,sw_net_sice, sea_salt_film, sea_salt_jet,                &
  up_sice_weighted_cat, up_sice_weighted,                                      &
  albedo_sice_weighted_cat, albedo_sice_weighted,                              &
  salt_dim1, salt_dim2, salt_dim3,                                             &
  cos_zenith_angle, day_fraction, cos_zen_rts, day_frac_rts, sol_azm_rts,      &
  ! LW diagnostic fields
  d_mass, density, layer_heat_capacity,                                        &
  p_layer_boundaries, p_layer_centres, p_extra_layer,                          &
  t_layer_boundaries, t_extra_layer,                                           &
  olr, lw_down, ozone, O3_trop_level, O3_trop_height,                          &
  T_trop_level, T_trop_height,                                                 &
  STASHwork)

USE def_diag, ONLY: StrDiag
USE def_spectrum, ONLY: StrSpecData
USE timestep_mod, ONLY: timestep
USE rad_input_mod, ONLY: l_rad_perturb
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE solinc_data, ONLY:                                                         &
    horiz_ang, n_horiz_ang, sky, l_skyview, slope_aspect, slope_angle, l_orog
USE ereport_mod, ONLY: ereport
USE submodel_mod, ONLY: atmos_im
USE stash_array_mod, ONLY:                                                     &
    len_stlist, sf_calc, stindex, stlist, num_stash_levels, stash_levels, si,  &
    stash_pseudo_levels, num_stash_pseudo
USE errormessagelength_mod, ONLY: errormessagelength
USE nlsizes_namelist_mod, ONLY: model_levels
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE level_heights_mod, ONLY: d_layer

IMPLICIT NONE


INTEGER, INTENT(IN) :: sect
!   STASH section (1=SW radiation, 2=LW radiation)
INTEGER, INTENT(IN) :: j_rad
!   Number of call to radiation (1=main call, 2=increment or forcing call)

TYPE (StrDiag),     INTENT(IN) :: diag
TYPE (StrSpecData), INTENT(IN) :: spectrum

LOGICAL, INTENT(IN) :: at_extremity(4)
!   Indicates if this processor is at north,
!   south, east or west of the processor grid

INTEGER, INTENT(IN) :: row_length       ! number of points on a row
INTEGER, INTENT(IN) :: rows             ! number of rows in a theta field
INTEGER, INTENT(IN) :: cloud_levels     ! number of cloudy levels
INTEGER, INTENT(IN) :: ozone_levels     ! number of levels where ozone is held
INTEGER, INTENT(IN) :: ntiles           ! number of land surface tiles
INTEGER, INTENT(IN) :: salt_dim1        !
INTEGER, INTENT(IN) :: salt_dim2        ! Dimensions for sea-salt aerosol
INTEGER, INTENT(IN) :: salt_dim3        !

INTEGER, INTENT(IN) :: i_off
!   Offset to diagnostics in multiple calls to radiation

REAL, INTENT(IN) :: T_n(row_length, rows, model_levels)
REAL, INTENT(IN) :: T_inc(row_length, rows, model_levels)
REAL, INTENT(IN) :: q_n(row_length, rows, model_levels)
REAL, INTENT(IN) :: qcl_n(row_length, rows, model_levels)
REAL, INTENT(IN) :: cf_n(row_length, rows, model_levels)
REAL, INTENT(IN) :: cfl_n(row_length, rows, model_levels)
REAL, INTENT(IN) :: T_latest(row_length, rows, model_levels)
REAL, INTENT(IN) :: q_latest(row_length, rows, model_levels)
REAL, INTENT(IN) :: qcl_latest(row_length, rows, model_levels)
REAL, INTENT(IN) :: cf_latest(row_length, rows, model_levels)
REAL, INTENT(IN) :: cfl_latest(row_length, rows, model_levels)
REAL, INTENT(IN) :: T_incr_diagnostic(row_length,rows,model_levels)
REAL, INTENT(IN) :: surfnet (row_length, rows)
REAL, INTENT(IN) :: netsea(row_length, rows)
REAL, INTENT(IN) :: obs_solid_angle(row_length, rows)

REAL, INTENT(IN) :: itoasw (row_length, rows)
REAL, INTENT(IN) :: surfsw_cor(row_length,rows)
REAL, INTENT(IN) :: toasw_cor(row_length,rows)
REAL, INTENT(IN) :: surfdir_cor(row_length,rows)
REAL, INTENT(IN) :: surfdif_cor(row_length,rows)
REAL, INTENT(IN) :: flux_below_690nm_surf(row_length, rows)
REAL, INTENT(IN) :: sw_net_land(row_length, rows)
REAL, INTENT(IN) :: sw_net_sice(row_length, rows)
REAL, INTENT(IN) :: photosynth_act_rad(row_length, rows)
REAL, INTENT(IN) :: flxdirparsurf(row_length, rows)
REAL, INTENT(IN) :: sea_salt_film(salt_dim1, salt_dim2, salt_dim3)
REAL, INTENT(IN) :: sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)
REAL, INTENT(IN) :: f_orog(row_length,rows)

REAL, INTENT(IN) :: olr (row_length, rows)
REAL, INTENT(IN) :: lw_down(row_length, rows)
REAL, INTENT(IN) :: ozone(row_length,rows,ozone_levels)
REAL, INTENT(IN) :: O3_trop_level(row_length,rows)
REAL, INTENT(IN) :: O3_trop_height(row_length,rows)
REAL, INTENT(IN) :: T_trop_level(row_length,rows)
REAL, INTENT(IN) :: T_trop_height(row_length,rows)

! Diagnostics weighted by sea ice concentration:
REAL, INTENT(IN) ::                                                            &
   up_sice_weighted_cat(row_length,rows,nice_use),                             &
   down_sice_weighted_cat(row_length,rows,nice_use),                           &
   up_sice_weighted(row_length,rows),                                          &
   down_sice_weighted(row_length,rows),                                        &
   albedo_sice_weighted_cat(row_length,rows,nice_use),                         &
   albedo_sice_weighted(row_length,rows)

! Cosine solar zenith angle for model timestep
REAL, INTENT(IN) :: cos_zenith_angle(row_length,rows)
! Fraction of model timestep for which sun is above the horizon
REAL, INTENT(IN) :: day_fraction(row_length,rows)
! Cosine solar zenith angle for radiation timestep
REAL, INTENT(IN) :: cos_zen_rts(row_length,rows)
! Fraction of radiation timestep for which sun is above the horizon
REAL, INTENT(IN) :: day_frac_rts(row_length,rows)
! Solar azimuth angle for radiation timestep
REAL, INTENT(IN) :: sol_azm_rts(row_length,rows)

! Fields calculated by set_thermodynamic:
REAL, INTENT(IN) :: d_mass(row_length, rows, model_levels+1)
!   Mass of layer (kg m-2)
REAL, INTENT(IN) :: density(row_length, rows, model_levels+1)
!   Density of layer (kg m-3)
REAL, INTENT(IN) :: layer_heat_capacity(row_length, rows, model_levels)
!   Heat capacity of layer
REAL, INTENT(IN) :: p_layer_boundaries(row_length, rows, 0:model_levels)
!   pressure at layer boundaries
REAL, INTENT(IN) :: p_layer_centres(row_length, rows, 0:model_levels)
!   pressure at layer centres
REAL, INTENT(IN) :: p_extra_layer(row_length, rows)
!   Pressure at centre of extra top layer
REAL, INTENT(IN) :: t_layer_boundaries(row_length, rows, 0:model_levels)
!   Temperature at layer boundaries
REAL, INTENT(IN) :: t_extra_layer(row_length, rows)
!   Temperature at centre of extra top layer


! Diagnostic variables
REAL ::  STASHwork(*)    ! STASH workspace


! Local variables
REAL :: work_3d(row_length, rows, model_levels)

INTEGER :: i, j, k, pslevel
INTEGER :: icode           ! Return code  =0 Normal exit  >1 Error

CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'DIAGNOSTICS_RAD'

INTEGER :: item                    ! STASH item
INTEGER, PARAMETER :: im_index=1   ! internal model index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

! Availability codes
INTEGER, PARAMETER :: ip_main_call = 1
  ! Diagnostics that are not contained in the diag structure can only be
  ! output once and are only valid for the main radiation call.
INTEGER, PARAMETER :: ip_main_and_forcing_call = 2
  ! Diagnostics that are contained in the diag structure and are available
  ! for the main and forcing calls but not the cloud increment calls.
INTEGER, PARAMETER :: ip_all_calls = 3
  ! Diagnostics that are contained in the diag structure and are available
  ! for the main, forcing and cloud increment calls.

LOGICAL :: l_main = .TRUE.

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
icode = 0  ! Initialise error status

CALL diagnostics_common()
IF (sect == 1) CALL diagnostics_sw()
IF (sect == 2) CALL diagnostics_lw()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

CONTAINS


SUBROUTINE diagnostics_common()
! Output diagnostics common to SW and LW
IMPLICIT NONE

! ----------------------------------------------------------------------
! Availability : ip_main_call
! First treat diagnostics that are not contained in the diag
! structure. These can only be output once and are valid after the
! last call to radiation where j_rad==1
! ----------------------------------------------------------------------

! Temperature
item = 4
IF (sf_calc(item,sect) .AND. j_rad == 1) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work_3d,T_n,T_inc)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work_3d(i,j,k) = T_n(i,j,k) + T_inc(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  CALL output_diag_3d(item, ip_main_call, l_main, work_3d, model_levels)
END IF

! Solid angle of grid-box as seen from an observer at 1AU.
CALL output_diag(130, ip_main_call, l_main, obs_solid_angle)

! Temperature Increment: increment diagnostics = modified - previous
CALL output_diag_3d(161, ip_main_call, l_main, T_incr_diagnostic, model_levels)

! Temperature increment including condensation
CALL output_diag_3d_inc(181, ip_main_call, l_main, T_latest, T_n, model_levels)

! Vapour increment
CALL output_diag_3d_inc(182, ip_main_call, l_main, q_latest, q_n, model_levels)

! Liquid water content increment
CALL output_diag_3d_inc(183, ip_main_call, l_main, qcl_latest, qcl_n, &
  model_levels)

! Total cloud fraction increment
CALL output_diag_3d_inc(192, ip_main_call, l_main, cf_latest, cf_n, &
  model_levels)

! Liquid cloud fraction increment
CALL output_diag_3d_inc(193, ip_main_call, l_main, cfl_latest, cfl_n, &
  model_levels)

! Liquid water content increment: positive
CALL output_diag_3d_pinc(194, ip_main_call, l_main, qcl_latest, qcl_n, &
  model_levels)

! Liquid water content increment: negative
CALL output_diag_3d_ninc(195, ip_main_call, l_main, qcl_latest, qcl_n, &
  model_levels)

! Liquid cloud fraction increment: positive
CALL output_diag_3d_pinc(198, ip_main_call, l_main, cfl_latest, cfl_n, &
  model_levels)

! Liquid cloud fraction increment: negative
CALL output_diag_3d_ninc(199, ip_main_call, l_main, cfl_latest, cfl_n, &
  model_levels)

! Surface net downwards flux
CALL output_diag(201, ip_main_call, l_main, surfnet)

! Net flux over sea
CALL output_diag(203, ip_main_call, l_main, netsea)

! Heating: radiation temperature increment per timestep / timestep
item = 232
IF (sf_calc(item,sect) .AND. j_rad == 1) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work_3d,T_incr_diagnostic,timestep)
  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        work_3d(i,j,k) = T_incr_diagnostic(i,j,k) / timestep
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END PARALLEL DO
  CALL output_diag_3d(item, ip_main_call, l_main, work_3d, model_levels)
END IF

! Downward flux on ice categories weighted by ice fraction
CALL output_diag_pl(500, ip_main_call, l_main, down_sice_weighted_cat, nice)

! Downward flux weighted by ice fraction
CALL output_diag(501, ip_main_call, l_main, down_sice_weighted)


! ----------------------------------------------------------------------
! Availability : ip_main_and_forcing_call
! Diagnostics contained within the diag structure are available for
! the main and forcing calls to radiation.
! For the incremental time-stepping scheme many of the diagnostics are
! not calculated on the "cloud only" radiation calls. For these the
! diagnostics from the last full radiation call are used. 
! ----------------------------------------------------------------------

! Clear-sky upward flux on levels
CALL output_diag_3d(219+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clear, diag%flux_up_clear, &
  model_levels+1)

! Clear-sky downward flux on levels
CALL output_diag_3d(220+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clear, diag%flux_down_clear, &
  model_levels+1)

! Clear-sky at top-of-atmosphere weighted by clear area
CALL output_diag(228+i_off, ip_main_and_forcing_call, &
  diag%l_toa_clear_weighted, diag%toa_clear_weighted)

! Total clear-sky area
CALL output_diag(229+i_off, ip_main_and_forcing_call, &
  diag%l_total_clear_area, diag%total_clear_area)

! Clear-sky Heating Rates
CALL output_diag_3d(233+i_off, ip_main_and_forcing_call, &
  diag%l_clear_hr, diag%clear_hr, &
  model_levels)

! Total aerosol optical depth in layers and bands
CALL output_diag_4d(506+i_off, ip_main_and_forcing_call, &
  diag%l_aerosol_optical_depth, diag%aerosol_optical_depth, &
  model_levels, spectrum%basic%n_band)

! Total aerosol scattering optical depth in layers and bands
CALL output_diag_4d(507+i_off, ip_main_and_forcing_call, &
  diag%l_aerosol_scat_optical_depth, diag%aerosol_scat_optical_depth, &
  model_levels, spectrum%basic%n_band)

! Total aerosol asymmetry in layers and bands
! weighted by scattering optical depth
CALL output_diag_4d(508+i_off, ip_main_and_forcing_call, &
  diag%l_aerosol_asymmetry_scat, diag%aerosol_asymmetry_scat, &
  model_levels, spectrum%basic%n_band)

! Upward flux on levels and bands
CALL output_diag_4d(509+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_band, diag%flux_up_band, &
  model_levels+1, spectrum%basic%n_band)

! Downward flux on levels and bands
CALL output_diag_4d(510+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_band, diag%flux_down_band, &
  model_levels+1, spectrum%basic%n_band)

! Clear-sky upward flux on levels and bands
CALL output_diag_4d(511+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clear_band, diag%flux_up_clear_band, &
  model_levels+1, spectrum%basic%n_band)

! Clear-sky downward flux on levels and bands
CALL output_diag_4d(512+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clear_band, diag%flux_down_clear_band, &
  model_levels+1, spectrum%basic%n_band)

! Emission spectrum towards observer in Wm-2 at 1 AU
CALL output_diag_ps(513+i_off, ip_main_and_forcing_call, &
  diag%l_emission_spectrum, diag%emission_spectrum, &
  spectrum%basic%n_band)

! Clear-sky emission spectrum
CALL output_diag_ps(514+i_off, ip_main_and_forcing_call, &
  diag%l_emission_spectrum_clear, diag%emission_spectrum_clear, &
  spectrum%basic%n_band)

! Clean-air emission spectrum
CALL output_diag_ps(515+i_off, ip_main_and_forcing_call, &
  diag%l_emission_spectrum_clean, diag%emission_spectrum_clean, &
  spectrum%basic%n_band)
CALL output_diag_ps(516+i_off, ip_main_and_forcing_call, &
  diag%l_emission_spectrum_clear_clean, diag%emission_spectrum_clear_clean, &
  spectrum%basic%n_band)

! Clean-air upward flux on levels
CALL output_diag_3d(517+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clean, diag%flux_up_clean, &
  model_levels+1)

! Clean-air downward flux on levels
CALL output_diag_3d(518+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clean, diag%flux_down_clean, &
  model_levels+1)

! Clear-sky, clean-air upward flux on levels
CALL output_diag_3d(519+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clear_clean, diag%flux_up_clear_clean, &
  model_levels+1)

! Clear-sky, clean-air downward flux on levels
CALL output_diag_3d(520+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clear_clean, diag%flux_down_clear_clean, &
  model_levels+1)

! GHG forcing upward flux on levels
CALL output_diag_3d(521+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_forc, diag%flux_up_forc, &
  model_levels+1)

! GHG forcing downward flux on levels
CALL output_diag_3d(522+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_forc, diag%flux_down_forc, &
  model_levels+1)

! Clear-sky, GHG forcing upward flux on levels
CALL output_diag_3d(523+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clear_forc, diag%flux_up_clear_forc, &
  model_levels+1)

! Clear-sky, GHG forcing downward flux on levels
CALL output_diag_3d(524+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clear_forc, diag%flux_down_clear_forc, &
  model_levels+1)

! GHG forcing upward flux on levels and bands
CALL output_diag_4d(525+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_forc_band, diag%flux_up_forc_band, &
  model_levels+1, spectrum%basic%n_band)

! GHG forcing downward flux on levels and bands
CALL output_diag_4d(526+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_forc_band, diag%flux_down_forc_band, &
  model_levels+1, spectrum%basic%n_band)

! Clear-sky GHG forcing upward flux on levels and bands
CALL output_diag_4d(527+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clear_forc_band, diag%flux_up_clear_forc_band, &
  model_levels+1, spectrum%basic%n_band)

! Clear-sky GHG forcing downward flux on levels and bands
CALL output_diag_4d(528+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clear_forc_band, diag%flux_down_clear_forc_band, &
  model_levels+1, spectrum%basic%n_band)

! Clean-air band-by-band fluxes
CALL output_diag_4d(545+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clean_band, diag%flux_up_clean_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_4d(546+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clean_band, diag%flux_down_clean_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_4d(547+i_off, ip_main_and_forcing_call, &
  diag%l_flux_up_clear_clean_band, diag%flux_up_clear_clean_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_4d(548+i_off, ip_main_and_forcing_call, &
  diag%l_flux_down_clear_clean_band, diag%flux_down_clear_clean_band, &
  model_levels+1, spectrum%basic%n_band)

! EasyAerosol extinction 3D profiles
CALL output_diag_4d(550+i_off, ip_main_and_forcing_call, &
  diag%l_easyaerosol_extinction, diag%easyaerosol_extinction, &
  model_levels, spectrum%basic%n_band)

! EasyAerosol absorption 3D profiles
CALL output_diag_4d(551+i_off, ip_main_and_forcing_call, &
  diag%l_easyaerosol_absorption, diag%easyaerosol_absorption, &
  model_levels, spectrum%basic%n_band)

! EasyAerosol scattering 3D profiles
CALL output_diag_4d(552+i_off, ip_main_and_forcing_call, &
  diag%l_easyaerosol_scattering, diag%easyaerosol_scattering, &
  model_levels, spectrum%basic%n_band)

! EasyAerosol asymmetry times scattering 3D profiles
CALL output_diag_4d(553+i_off, ip_main_and_forcing_call, &
  diag%l_easyaerosol_asytimscat, diag%easyaerosol_asytimscat, &
  model_levels, spectrum%basic%n_band)

! Radiances at TOA
CALL output_diag_ps(297+i_off, ip_main_and_forcing_call, &
  diag%l_toa_radiance, diag%toa_radiance, &
  spectrum%basic%n_band)


! ----------------------------------------------------------------------
! Availability : ip_all_calls
! The following diagnostics are available on "cloud only" radiation
! calls for the incremental time-stepping scheme, and all calls for
! other schemes.
! ----------------------------------------------------------------------

! Upward flux on levels
CALL output_diag_3d(217+i_off, ip_all_calls, &
  diag%l_flux_up, diag%flux_up, &
  model_levels+1)

! Downward flux on levels
CALL output_diag_3d(218+i_off, ip_all_calls, &
  diag%l_flux_down, diag%flux_down, &
  model_levels+1)

! Net flux at Tropopause
CALL output_diag(237+i_off, ip_all_calls, &
  diag%l_net_flux_trop, diag%net_flux_trop)


END SUBROUTINE diagnostics_common


SUBROUTINE diagnostics_sw()
! Output diagnostics that are unique to the SW
IMPLICIT NONE

! ----------------------------------------------------------------------
! Availability : ip_main_call
! ----------------------------------------------------------------------

! Horizon angles
IF (l_skyview) THEN
  DO i = 1, n_horiz_ang
    item = 100+i
    CALL output_diag(100+i, ip_main_call, l_main, horiz_ang(:,:,1,i))
  END DO
END IF

! Cosine solar zenith angle for model timestep
CALL output_diag(140, ip_main_call, l_main, cos_zenith_angle)

! Fraction of model timestep for which sun is above the horizon
CALL output_diag(141, ip_main_call, l_main, day_fraction)

! Cosine solar zenith angle for radiation timestep
CALL output_diag(142, ip_main_call, l_main, cos_zen_rts)

! Fraction of radiation timestep for which sun is above the horizon
CALL output_diag(143, ip_main_call, l_main, day_frac_rts)

! Incoming SW at TOA
CALL output_diag(207, ip_main_call, l_main, itoasw)

! Surface SW corrected for solar zenith angle
CALL output_diag(202, ip_main_call, l_main, surfsw_cor)

! Outgoing SW corrected for solar zenith angle
CALL output_diag(205, ip_main_call, l_main, toasw_cor)

! Direct surface SW corrected for solar zenith angle
CALL output_diag(215, ip_main_call, l_main, surfdir_cor)

! Diffuse surface SW corrected for solar zenith angle
CALL output_diag(216, ip_main_call, l_main, surfdif_cor)

! Flux Below 690nm at Surface
CALL output_diag(204, ip_main_call, l_main, flux_below_690nm_surf)

! Surface Net Land
CALL output_diag(257, ip_main_call, l_main, sw_net_land)

! SW Net Sice
CALL output_diag(258, ip_main_call, l_main, sw_net_sice)

! Photosynth_act_rad (total PAR flux at surface)
CALL output_diag(290, ip_main_call, l_main, photosynth_act_rad)

! Flux_direct_par (direct component of PAR flux at surface)
CALL output_diag(291, ip_main_call, l_main, flxdirparsurf)

! Solar azimuth angle for radiation timestep (radians clockwise from grid north)
CALL output_diag(292, ip_main_call, l_main, sol_azm_rts)

! Slope Aspect
CALL output_diag(293, ip_main_call, l_orog, slope_aspect)

! Slope Angle
CALL output_diag(294, ip_main_call, l_orog, slope_angle)

! F_orog
CALL output_diag(296, ip_main_call, l_orog, f_orog)

! Sea Salt film
CALL output_diag_3d(247, ip_main_call, l_main, sea_salt_film, model_levels)

! Sea Salt Jet
CALL output_diag_3d(248, ip_main_call, l_main, sea_salt_jet, model_levels)

! Upward SW on ice categories weighted by ice fraction
CALL output_diag_pl(502, ip_main_call, l_main, up_sice_weighted_cat, nice)

! Upward SW weighted by ice fraction
CALL output_diag(503, ip_main_call, l_main, up_sice_weighted)

! Sea ice albdeo on ice categories weighted by ice fraction
CALL output_diag_pl(504, ip_main_call, l_main, albedo_sice_weighted_cat, nice)

! Sea ice albedo weighted by ice fraction
CALL output_diag(505, ip_main_call, l_main, albedo_sice_weighted)


! ----------------------------------------------------------------------
! Availability : ip_main_and_forcing_call
! ----------------------------------------------------------------------

! Mask for radiation calculations
CALL output_diag(200+i_off, ip_main_and_forcing_call, &
  diag%l_rad_mask, diag%rad_mask)

! Outwards Solar Clear Flux at TOA
CALL output_diag(209+i_off, ip_main_and_forcing_call, &
  diag%l_solar_out_clear, diag%solar_out_clear)

! Surface Down Clear
CALL output_diag(210+i_off, ip_main_and_forcing_call, &
  diag%l_surf_down_clr, diag%surf_down_clr)

! Surface up Clear
CALL output_diag(211+i_off, ip_main_and_forcing_call, &
  diag%l_surf_up_clr, diag%surf_up_clr)

! Direct UV-Flux
CALL output_diag_3d(212+i_off, ip_main_and_forcing_call, &
  diag%l_uvflux_direct, diag%uvflux_direct, &
  model_levels+1)

! UV Flux Up
CALL output_diag_3d(213+i_off, ip_main_and_forcing_call, &
  diag%l_uvflux_up, diag%uvflux_up, &
  model_levels+1)

! Down UV Flux
CALL output_diag_3d(214+i_off, ip_main_and_forcing_call, &
  diag%l_uvflux_down, diag%uvflux_down, &
  model_levels+1)

! Surface Down UV Flux
CALL output_diag(288+i_off, ip_main_and_forcing_call, &
  diag%l_surf_uv, diag%surf_uv)

! Surface Down clear-sky UV Flux
CALL output_diag(289+i_off, ip_main_and_forcing_call, &
  diag%l_surf_uv_clr, diag%surf_uv_clr)

! Cloud Extinction
CALL output_diag_3d(262+i_off, ip_main_and_forcing_call, &
  diag%l_cloud_extinction, diag%cloud_extinction, &
  cloud_levels)

! Cloud Weight Extinction
CALL output_diag_3d(263+i_off, ip_main_and_forcing_call, &
  diag%l_cloud_weight_extinction, &
  diag%cloud_weight_extinction, &
  cloud_levels)

! Large-Scale Cloud Extinction
CALL output_diag_3d(264+i_off, ip_main_and_forcing_call, &
  diag%l_ls_cloud_extinction, diag%ls_cloud_extinction, &
  cloud_levels)

! Large-Scale Cloud Weight Extinction
CALL output_diag_3d(265+i_off, ip_main_and_forcing_call, &
  diag%l_ls_cloud_weight_extinction, &
  diag%ls_cloud_weight_extinction, &
  cloud_levels)

! Convective Cloud Extinction
CALL output_diag_3d(266+i_off, ip_main_and_forcing_call, &
  diag%l_cnv_cloud_extinction, diag%cnv_cloud_extinction, &
  cloud_levels)

! Convective Cloud weight Extinction
CALL output_diag_3d(267+i_off, ip_main_and_forcing_call, &
  diag%l_cnv_cloud_weight_extinction, &
  diag%cnv_cloud_weight_extinction, &
  cloud_levels)

! Re. Strat
CALL output_diag_3d(221+i_off, ip_main_and_forcing_call, &
  diag%re_strat_flag, diag%re_strat, &
  cloud_levels)

! Wgt. Strat
CALL output_diag_3d(223+i_off, ip_main_and_forcing_call, &
  diag%wgt_strat_flag, diag%wgt_strat, &
  cloud_levels)

! LWP. Strat
CALL output_diag_3d(224+i_off, ip_main_and_forcing_call, &
  diag%lwp_strat_flag, diag%lwp_strat, &
  cloud_levels)

! Re. Conv
CALL output_diag_3d(225+i_off, ip_main_and_forcing_call, &
  diag%re_conv_flag, diag%re_conv, &
  cloud_levels)

! Wgt. Conv
CALL output_diag_3d(226+i_off, ip_main_and_forcing_call, &
  diag%wgt_conv_flag, diag%wgt_conv, &
  cloud_levels)

! Ntot. Diag
CALL output_diag_3d(241+i_off, ip_main_and_forcing_call, &
  diag%ntot_diag_flag, diag%ntot_diag, &
  cloud_levels)

! Strat. LWC Diag
CALL output_diag_3d(242+i_off, ip_main_and_forcing_call, &
  diag%strat_lwc_diag_flag, diag%strat_lwc_diag, &
  cloud_levels)

! SO4 Cloud Condensation Nuclei
CALL output_diag_3d(243+i_off, ip_main_and_forcing_call, &
  diag%so4_ccn_diag_flag, diag%so4_ccn_diag, &
  cloud_levels)

! Cond. Samp. Wgt
CALL output_diag_3d(244+i_off, ip_main_and_forcing_call, &
  diag%cond_samp_wgt_flag, diag%cond_samp_wgt, &
  cloud_levels)

! Weighted Re
CALL output_diag(245+i_off, ip_main_and_forcing_call, &
  diag%weighted_re_flag, diag%weighted_re)

! Sum Weighted Re
CALL output_diag(246+i_off, ip_main_and_forcing_call, &
  diag%sum_weight_re_flag, diag%sum_weight_re)

! Weighted Warm Re
CALL output_diag(254+i_off, ip_main_and_forcing_call, &
  diag%wgtd_warm_re_flag, diag%weighted_warm_re)

! Sum Weighted Warm Re
CALL output_diag(255+i_off, ip_main_and_forcing_call, &
  diag%sum_wgt_warm_re_flag, diag%sum_weight_warm_re)

! Weighted CDNC @cloud top
CALL output_diag(298+i_off, ip_main_and_forcing_call, &
  diag%cdnc_ct_diag_flag, diag%cdnc_ct_diag)

! Weight for CDNC @cloud top
CALL output_diag(299+i_off, ip_main_and_forcing_call, &
  diag%cdnc_ct_weight_flag, diag%cdnc_ct_weight)

! Nc. Diag
CALL output_diag(280+i_off, ip_main_and_forcing_call, &
  diag%Nc_diag_flag, diag%Nc_diag)

! Nc. Weight
CALL output_diag(281+i_off, ip_main_and_forcing_call, &
  diag%Nc_weight_flag, diag%Nc_weight)

! Direct Surface Albedo on SW bands
CALL output_diag_ps(268+i_off, ip_main_and_forcing_call, &
  diag%l_direct_albedo, diag%direct_albedo, &
  spectrum%basic%n_band)

! Diffuse Surface Albedo on SW bands. Stash 269,1
CALL output_diag_ps(269+i_off, ip_main_and_forcing_call, &
  diag%l_diffuse_albedo, diag%diffuse_albedo, &
  spectrum%basic%n_band)

! Direct fluxes for spherical geometry
CALL output_diag_03d(273+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_sph, diag%flux_direct_clear_sph, &
  model_levels+1)
CALL output_diag_04d(537+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_sph_band, diag%flux_direct_sph_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_04d(538+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_sph_band, diag%flux_direct_clear_sph_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_3d(277+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_div, diag%flux_direct_clear_div, &
  model_levels)
CALL output_diag_4d(541+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_div_band, diag%flux_direct_div_band, &
  model_levels, spectrum%basic%n_band)
CALL output_diag_4d(542+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_div_band, diag%flux_direct_clear_div_band, &
  model_levels, spectrum%basic%n_band)

! Transmission spectrum
CALL output_diag_ps(555+i_off, ip_main_and_forcing_call, &
  diag%l_transmission_spectrum, diag%transmission_spectrum, &
  spectrum%basic%n_band)
CALL output_diag_ps(556+i_off, ip_main_and_forcing_call, &
  diag%l_transmission_spectrum_clear, diag%transmission_spectrum_clear, &
  spectrum%basic%n_band)

! Clean-air direct fluxes for spherical geometry
CALL output_diag_03d(274+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clean_sph, diag%flux_direct_clean_sph, &
  model_levels+1)
CALL output_diag_03d(275+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_clean_sph, diag%flux_direct_clear_clean_sph, &
  model_levels+1)
CALL output_diag_04d(539+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clean_sph_band, diag%flux_direct_clean_sph_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_04d(540+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_clean_sph_band, &
  diag%flux_direct_clear_clean_sph_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_3d(278+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clean_div, diag%flux_direct_clean_div, &
  model_levels)
CALL output_diag_3d(279+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_clean_div, diag%flux_direct_clear_clean_div, &
  model_levels)
CALL output_diag_4d(543+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clean_div_band, diag%flux_direct_clean_div_band, &
  model_levels, spectrum%basic%n_band)
CALL output_diag_4d(544+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_clean_div_band, &
  diag%flux_direct_clear_clean_div_band, &
  model_levels, spectrum%basic%n_band)

! Clean-air transmission spectrum
CALL output_diag_ps(557+i_off, ip_main_and_forcing_call, &
  diag%l_transmission_spectrum_clean, &
  diag%transmission_spectrum_clean, &
  spectrum%basic%n_band)
CALL output_diag_ps(558+i_off, ip_main_and_forcing_call, &
  diag%l_transmission_spectrum_clear_clean, &
  diag%transmission_spectrum_clear_clean, &
  spectrum%basic%n_band)

! GHG forcing direct fluxes for spherical geometry
CALL output_diag_03d(529+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_sph_forc, diag%flux_direct_sph_forc, &
  model_levels+1)
CALL output_diag_03d(530+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_sph_forc, diag%flux_direct_clear_sph_forc, &
  model_levels+1)
CALL output_diag_3d(531+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_div_forc, diag%flux_direct_div_forc, &
  model_levels)
CALL output_diag_3d(532+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_div_forc, diag%flux_direct_clear_div_forc, &
  model_levels)
CALL output_diag_04d(533+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_sph_forc_band, diag%flux_direct_sph_forc_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_04d(534+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_sph_forc_band, &
  diag%flux_direct_clear_sph_forc_band, &
  model_levels+1, spectrum%basic%n_band)
CALL output_diag_4d(535+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_div_forc_band, diag%flux_direct_div_forc_band, &
  model_levels, spectrum%basic%n_band)
CALL output_diag_4d(536+i_off, ip_main_and_forcing_call, &
  diag%l_flux_direct_clear_div_forc_band, &
  diag%flux_direct_clear_div_forc_band, &
  model_levels, spectrum%basic%n_band)


! ----------------------------------------------------------------------
! Availability : ip_all_calls
! ----------------------------------------------------------------------

! Outwards Solar Flux at TOA
CALL output_diag(208+i_off, ip_all_calls, &
  diag%l_solar_out_toa, diag%solar_out_toa)

! Surface Down Flux
CALL output_diag(235+i_off, ip_all_calls, &
  diag%l_surface_down_flux, diag%surface_down_flux)

! Direct flux on levels
CALL output_diag_3d(230+i_off, ip_all_calls, &
  diag%l_flux_direct, diag%flux_direct, &
  model_levels+1)

! Diffuse flux on levels
CALL output_diag_3d(231+i_off, ip_all_calls, &
  diag%l_flux_diffuse, diag%flux_diffuse, &
  model_levels+1)

! SW Up Flux at Tropopause
CALL output_diag(238+i_off, ip_all_calls, &
  diag%l_up_flux_trop, diag%up_flux_trop)

! Fl_solid_below_690nm_surf
CALL output_diag(259+i_off, ip_all_calls, &
  diag%l_FlxSolBelow690nmSurf, diag%FlxSolBelow690nmSurf)

! Fl_sea_below_690nm_surf
CALL output_diag(260+i_off, ip_all_calls, &
  diag%l_FlxSeaBelow690nmSurf, diag%FlxSeaBelow690nmSurf)

! Orographic Correction
CALL output_diag(295+i_off, ip_all_calls, &
  diag%l_orog_corr, diag%orog_corr)

! VIS Albedo scaling to obs, on land surface tiles
CALL output_diag_ps(270+i_off, ip_all_calls, &
  diag%l_vis_albedo_sc, diag%vis_albedo_sc, &
  ntiles)

! NIR Albedo scaling to obs, on land surface tiles
CALL output_diag_ps(271+i_off, ip_all_calls, &
  diag%l_nir_albedo_sc, diag%nir_albedo_sc, &
  ntiles)

! Direct flux for spherical geometry
CALL output_diag_03d(272+i_off, ip_all_calls, &
  diag%l_flux_direct_sph, diag%flux_direct_sph, &
  model_levels+1)

! Direct flux divergence for spherical geometry
CALL output_diag_3d(276+i_off, ip_all_calls, &
  diag%l_flux_direct_div, diag%flux_direct_div, &
  model_levels)

END SUBROUTINE diagnostics_sw


SUBROUTINE diagnostics_lw()
 ! Output diagnostics that are unique to the LW
IMPLICIT NONE

! ----------------------------------------------------------------------
! Availability : ip_main_call
! ----------------------------------------------------------------------

! Skyview factor
CALL output_diag(101, ip_main_call, l_skyview, sky)

! Layer mass
CALL output_diag_3d(110, ip_main_call, l_main, d_mass, model_levels+1)

! Layer density
CALL output_diag_3d(111, ip_main_call, l_main, density, model_levels+1)

! Layer heat capacity
CALL output_diag_3d(112, ip_main_call, l_main, layer_heat_capacity, &
  model_levels)

! Layer pressure
item = 113
IF (sf_calc(item,sect) .AND. j_rad == 1) THEN
!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work_3d,p_layer_centres,p_extra_layer)
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work_3d(i,j,k) = p_layer_centres(i, j, k)
      END DO
    END DO
  END DO
!$OMP END DO
  k = model_levels+1
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      work_3d(i,j,k) = p_extra_layer(i, j)
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
  CALL output_diag_3d(item, ip_main_call, l_main, work_3d, model_levels+1)
END IF

! Layer temperature
item = 114
IF (sf_calc(item,sect) .AND. j_rad == 1) THEN
!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work_3d,t_n,t_extra_layer)
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        work_3d(i,j,k) = t_n(i, j, k)
      END DO
    END DO
  END DO
!$OMP END DO
  k = model_levels+1
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      work_3d(i,j,k) = t_extra_layer(i, j)
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL
  CALL output_diag_3d(item, ip_main_call, l_main, work_3d, model_levels+1)
END IF

! Level pressure
item = 115
IF (sf_calc(item,sect) .AND. j_rad == 1) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,p_layer_boundaries,work_3d)
  DO k = 1, model_levels+1
    DO j = 1, rows
      DO i = 1, row_length
        work_3d(i,j,k) = p_layer_boundaries(i, j, k-1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  CALL output_diag_3d(item, ip_main_call, l_main, work_3d, model_levels+1)
END IF

! Level temperature
item = 116
IF (sf_calc(item,sect) .AND. j_rad == 1) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(model_levels,rows,row_length,work_3d,t_layer_boundaries)
  DO k = 1, model_levels+1
    DO j = 1, rows
      DO i = 1, row_length
        work_3d(i,j,k) = t_layer_boundaries(i, j, k-1)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  CALL output_diag_3d(item, ip_main_call, l_main, work_3d, model_levels+1)
END IF

! OLR (outgoing longwave radiation at top-of-atmosphere)
CALL output_diag(205, ip_main_call, l_main, olr)

! Surface Down LW Flux
CALL output_diag(207, ip_main_call, l_main, lw_down)

! Ozone
CALL output_diag_3d(260, ip_main_call, l_main, ozone, ozone_levels)

! Ozone troposphere level
CALL output_diag(280, ip_main_call, l_main, O3_trop_level)

! Ozone troposphere height
CALL output_diag(281, ip_main_call, l_main, O3_trop_height)

! Temperature troposphere level
CALL output_diag(282, ip_main_call, l_main, T_trop_level)

! Temperature troposphere height
CALL output_diag(283, ip_main_call, l_main, T_trop_height)

! Depth of radiative layers (metres)
CALL output_diag_3d(505, ip_main_call, l_main, d_layer, model_levels+1)


! ----------------------------------------------------------------------
! Availability : ip_main_and_forcing_call
! ----------------------------------------------------------------------

! Clear OLR
CALL output_diag(206+i_off, ip_main_and_forcing_call, &
  diag%l_clear_olr, diag%clear_olr)

! Clear Sky Surface Down Flux
CALL output_diag(208+i_off, ip_main_and_forcing_call, &
  diag%l_surf_down_clr, diag%surf_down_clr)


! UKCA absorption aerosol optical depth diagnostics

! UKCA AAOD Aitken soluble
CALL output_diag_ps(240+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_ukca_ait_sol, diag%aaod_ukca_ait_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AAOD accumulation soluble
CALL output_diag_ps(241+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_ukca_acc_sol, diag%aaod_ukca_acc_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AAOD coarse soluble
CALL output_diag_ps(242+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_ukca_cor_sol, diag%aaod_ukca_cor_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AAOD Aitken insoluble
CALL output_diag_ps(243+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_ukca_ait_ins, diag%aaod_ukca_ait_ins, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AAOD accumulation insoluble
CALL output_diag_ps(244+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_ukca_acc_ins, diag%aaod_ukca_acc_ins, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AAOD coarse insoluble
CALL output_diag_ps(245+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_ukca_cor_ins, diag%aaod_ukca_cor_ins, &
  spectrum%aerosol%n_aod_wavel)


! Cloud absorptivity diagnostics

! Cloud Absorptivity
CALL output_diag_3d(262+i_off, ip_main_and_forcing_call, &
  diag%l_cloud_absorptivity, diag%cloud_absorptivity, &
  cloud_levels)

! Cloud Weight Absorptivity
CALL output_diag_3d(263+i_off, ip_main_and_forcing_call, &
  diag%l_cloud_weight_absorptivity, diag%cloud_weight_absorptivity, &
  cloud_levels)

! Large Scale Cloud Absorptivity
CALL output_diag_3d(264+i_off, ip_main_and_forcing_call, &
  diag%l_ls_cloud_absorptivity, diag%ls_cloud_absorptivity, &
  cloud_levels)

! Large Scale Cloud Weight Absorptivity
CALL output_diag_3d(265+i_off, ip_main_and_forcing_call, &
  diag%l_ls_cloud_weight_absorptivity, diag%ls_cloud_weight_absorptivity, &
  cloud_levels)

! Convective Cloud Absorptivity
CALL output_diag_3d(266+i_off, ip_main_and_forcing_call, &
  diag%l_cnv_cloud_absorptivity, diag%cnv_cloud_absorptivity, &
  cloud_levels)

! Convective Cloud Weight Absorptivity
CALL output_diag_3d(267+i_off, ip_main_and_forcing_call, &
  diag%l_cnv_cloud_weight_absorptivity, diag%cnv_cloud_weight_absorptivity, &
  cloud_levels)


! CLASSIC aerosol optical depth diagnostics

! AOD Sulphate
CALL output_diag_ps(284+i_off, ip_main_and_forcing_call, &
  diag%l_aod_sulphate, diag%aod_sulphate, &
  spectrum%aerosol%n_aod_wavel)

! AOD Dust
CALL output_diag_ps(285+i_off, ip_main_and_forcing_call, &
  diag%l_aod_dust, diag%aod_dust, &
  spectrum%aerosol%n_aod_wavel)

! AOD Seasalt
CALL output_diag_ps(286+i_off, ip_main_and_forcing_call, &
  diag%l_aod_seasalt, diag%aod_seasalt, &
  spectrum%aerosol%n_aod_wavel)

! AOD Soot
CALL output_diag_ps(287+i_off, ip_main_and_forcing_call, &
  diag%l_aod_soot, diag%aod_soot, &
  spectrum%aerosol%n_aod_wavel)

! AOD Biomass
CALL output_diag_ps(288+i_off, ip_main_and_forcing_call, &
  diag%l_aod_biomass, diag%aod_biomass, &
  spectrum%aerosol%n_aod_wavel)

! AOD Biogenic
CALL output_diag_ps(289+i_off, ip_main_and_forcing_call, &
  diag%l_aod_biogenic, diag%aod_biogenic, &
  spectrum%aerosol%n_aod_wavel)

! AOD Fossil-fuel organic carbon
CALL output_diag_ps(295+i_off, ip_main_and_forcing_call, &
  diag%l_aod_ocff, diag%aod_ocff, &
  spectrum%aerosol%n_aod_wavel)

! AOD Delta aerosol
CALL output_diag_ps(296+i_off, ip_main_and_forcing_call, &
  diag%l_aod_delta, diag%aod_delta, &
  spectrum%aerosol%n_aod_wavel)

! AOD Nitrate
CALL output_diag_ps(297+i_off, ip_main_and_forcing_call, &
  diag%l_aod_nitrate, diag%aod_nitrate, &
  spectrum%aerosol%n_aod_wavel)

! Total AOD in Radiation
CALL output_diag_ps(298+i_off, ip_main_and_forcing_call, &
  diag%l_aod_total_radn, diag%aod_total_radn, &
  spectrum%aerosol%n_aod_wavel)

! Angstrom Exponent of the Total AOD in Radiation
CALL output_diag_ps(299+i_off, ip_main_and_forcing_call, &
  diag%l_angst_total_radn, diag%angst_total_radn, &
  spectrum%aerosol%n_aod_wavel)


! UKCA aerosol optical depth diagnostics

! UKCA AOD Aitken soluble
CALL output_diag_ps(300+i_off, ip_main_and_forcing_call, &
  diag%l_aod_ukca_ait_sol, diag%aod_ukca_ait_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AOD accum. soluble
CALL output_diag_ps(301+i_off, ip_main_and_forcing_call, &
  diag%l_aod_ukca_acc_sol, diag%aod_ukca_acc_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AOD coarse soluble
CALL output_diag_ps(302+i_off, ip_main_and_forcing_call, &
  diag%l_aod_ukca_cor_sol, diag%aod_ukca_cor_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AOD Aitken insoluble
CALL output_diag_ps(303+i_off, ip_main_and_forcing_call, &
  diag%l_aod_ukca_ait_ins, diag%aod_ukca_ait_ins, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AOD accum. insoluble
CALL output_diag_ps(304+i_off, ip_main_and_forcing_call, &
  diag%l_aod_ukca_acc_ins, diag%aod_ukca_acc_ins, &
  spectrum%aerosol%n_aod_wavel)

! UKCA AOD coarse insoluble
CALL output_diag_ps(305+i_off, ip_main_and_forcing_call, &
  diag%l_aod_ukca_cor_ins, diag%aod_ukca_cor_ins, &
  spectrum%aerosol%n_aod_wavel)

! UKCA stratospheric AOD Aitken soluble
CALL output_diag_ps(251+i_off, ip_main_and_forcing_call, &
  diag%l_sod_ukca_ait_sol, diag%sod_ukca_ait_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA stratospheric AOD accum. soluble
CALL output_diag_ps(252+i_off, ip_main_and_forcing_call, &
  diag%l_sod_ukca_acc_sol, diag%sod_ukca_acc_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA stratospheric AOD coarse soluble
CALL output_diag_ps(253+i_off, ip_main_and_forcing_call, &
  diag%l_sod_ukca_cor_sol, diag%sod_ukca_cor_sol, &
  spectrum%aerosol%n_aod_wavel)

! UKCA stratospheric AOD Aitken insoluble
CALL output_diag_ps(254+i_off, ip_main_and_forcing_call, &
  diag%l_sod_ukca_ait_ins, diag%sod_ukca_ait_ins, &
  spectrum%aerosol%n_aod_wavel)

! UKCA stratospheric AOD accum. insoluble
CALL output_diag_ps(255+i_off, ip_main_and_forcing_call, &
  diag%l_sod_ukca_acc_ins, diag%sod_ukca_acc_ins, &
  spectrum%aerosol%n_aod_wavel)

! UKCA stratospheric AOD coarse insoluble
CALL output_diag_ps(256+i_off, ip_main_and_forcing_call, &
  diag%l_sod_ukca_cor_ins, diag%sod_ukca_cor_ins, &
  spectrum%aerosol%n_aod_wavel)


! Convective core diagnostics

! Convective core 3D cloud fraction
CALL output_diag_3d(317+i_off, ip_main_and_forcing_call, &
  diag%l_ccore_clt_rad, diag%ccore_clt_rad, &
  model_levels)

! Convective core gridbox-mean liquid mixing ratio
CALL output_diag_3d(318+i_off, ip_main_and_forcing_call, &
  diag%l_ccore_qcl_rad, diag%ccore_qcl_rad, &
  model_levels)

! Convective core gridbox-mean ice mixing ratio
CALL output_diag_3d(319+i_off, ip_main_and_forcing_call, &
  diag%l_ccore_qcf_rad, diag%ccore_qcf_rad, &
  model_levels)

! Convective core liquid water path
CALL output_diag(395+i_off, ip_main_and_forcing_call, &
  diag%l_ccore_qcl_rad_path, diag%ccore_qcl_rad_path)

! Convective core ice water path
CALL output_diag(396+i_off, ip_main_and_forcing_call, &
  diag%l_ccore_qcf_rad_path, diag%ccore_qcf_rad_path)


! Prognostic aerosol optical depth diagnostics

! Prognositc Sulphate AOD
CALL output_diag_ps(421+i_off, ip_main_and_forcing_call, &
  diag%l_aod_prog_sulphate, diag%aod_prog_sulphate, &
  spectrum%aerosol%n_aod_wavel)

! Prognostic Dust AOD
CALL output_diag_ps(422+i_off, ip_main_and_forcing_call, &
  diag%l_aod_prog_dust, diag%aod_prog_dust, &
  spectrum%aerosol%n_aod_wavel)

! Diagnosed Seasalt AOD
CALL output_diag_ps(423+i_off, ip_main_and_forcing_call, &
  diag%l_aod_prog_seasalt, diag%aod_prog_seasalt, &
  spectrum%aerosol%n_aod_wavel)

! Prognostic Soot AOD
CALL output_diag_ps(424+i_off, ip_main_and_forcing_call, &
  diag%l_aod_prog_soot, diag%aod_prog_soot, &
  spectrum%aerosol%n_aod_wavel)

! Prognostic Biomass AOD
CALL output_diag_ps(425+i_off, ip_main_and_forcing_call, &
  diag%l_aod_prog_biomass, diag%aod_prog_biomass, &
  spectrum%aerosol%n_aod_wavel)

! AOD Fossil-fuel organic carbon
CALL output_diag_ps(426+i_off, ip_main_and_forcing_call, &
  diag%l_aod_prog_ocff, diag%aod_prog_ocff, &
  spectrum%aerosol%n_aod_wavel)

! Prognostic Nitrate AOD
CALL output_diag_ps(427+i_off, ip_main_and_forcing_call, &
  diag%l_aod_prog_nitrate, diag%aod_prog_nitrate, &
  spectrum%aerosol%n_aod_wavel)


! UKCA RADAER 3D aerosol-radiation diagnostics

! UKCA aerosol extinction 3D profiles
CALL output_diag_4d(530+i_off, ip_main_and_forcing_call, &
  diag%n_ukca_aerosol_ext > 0, diag%ukca_aerosol_ext, &
  model_levels, diag%n_ukca_aerosol_ext)

! UKCA aerosol absorption 3D profiles
CALL output_diag_4d(531+i_off, ip_main_and_forcing_call, &
  diag%n_ukca_aerosol_abs > 0, diag%ukca_aerosol_abs, &
  model_levels, diag%n_ukca_aerosol_abs)

! UKCA aerosol scattering 3D profiles
CALL output_diag_4d(532+i_off, ip_main_and_forcing_call, &
  diag%n_ukca_aerosol_sca > 0, diag%ukca_aerosol_sca, &
  model_levels, diag%n_ukca_aerosol_sca)

! UKCA aerosol asymmetry * scattering 3D profiles
CALL output_diag_4d(533+i_off, ip_main_and_forcing_call, &
  diag%n_ukca_aerosol_gsca > 0, diag%ukca_aerosol_gsca, &
  model_levels, diag%n_ukca_aerosol_gsca)


! CLASSIC RADAER 3D aerosol-radiation diagnostics

! CLASSIC aerosol extinction 3D profiles
CALL output_diag_4d(540+i_off, ip_main_and_forcing_call, &
  diag%n_clas_aerosol_ext > 0, diag%clas_aerosol_ext, &
  model_levels, diag%n_clas_aerosol_ext)

! CLASSIC aerosol absorption 3D profiles
CALL output_diag_4d(541+i_off, ip_main_and_forcing_call, &
  diag%n_clas_aerosol_abs > 0, diag%clas_aerosol_abs, &
  model_levels, diag%n_clas_aerosol_abs)

! CLASSIC aerosol scattering 3D profiles
CALL output_diag_4d(542+i_off, ip_main_and_forcing_call, &
  diag%n_clas_aerosol_sca > 0, diag%clas_aerosol_sca, &
  model_levels, diag%n_clas_aerosol_sca)


! CLASSIC absorption aerosol optical depth diagnostics

! AAOD Sulphate
CALL output_diag_ps(584+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_sulphate, diag%aaod_sulphate, &
  spectrum%aerosol%n_aod_wavel)

! AAOD Dust
CALL output_diag_ps(585+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_dust, diag%aaod_dust, &
  spectrum%aerosol%n_aod_wavel)

! AAOD Seasalt
CALL output_diag_ps(586+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_seasalt, diag%aaod_seasalt, &
  spectrum%aerosol%n_aod_wavel)

! AAOD Soot
CALL output_diag_ps(587+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_soot, diag%aaod_soot, &
  spectrum%aerosol%n_aod_wavel)

! AAOD Biomass
CALL output_diag_ps(588+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_biomass, diag%aaod_biomass, &
  spectrum%aerosol%n_aod_wavel)

! AAOD Biogenic
CALL output_diag_ps(589+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_biogenic, diag%aaod_biogenic, &
  spectrum%aerosol%n_aod_wavel)

! AAOD Fossil-fuel organic carbon
CALL output_diag_ps(595+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_ocff, diag%aaod_ocff, &
  spectrum%aerosol%n_aod_wavel)

! AAOD Nitrate
CALL output_diag_ps(597+i_off, ip_main_and_forcing_call, &
  diag%l_aaod_nitrate, diag%aaod_nitrate, &
  spectrum%aerosol%n_aod_wavel)


! ----------------------------------------------------------------------
! Availability : ip_all_calls
! ----------------------------------------------------------------------

! LW Down Flux at Tropopause
CALL output_diag(238+i_off, ip_all_calls, &
  diag%l_down_flux_trop, diag%down_flux_trop)

! Total Cloud Cover
CALL output_diag(204+i_off, ip_all_calls, &
  diag%l_total_cloud_cover, diag%total_cloud_cover)

! Total Cloud on Levels
CALL output_diag_3d(261+i_off, ip_all_calls, &
  diag%l_total_cloud_on_levels, diag%total_cloud_on_levels, &
  cloud_levels)


! Cloud water mixing ratios (3D)

! LS cloud liquid water mixing ratio
CALL output_diag_3d(308+i_off, ip_all_calls, &
  diag%l_ls_qcl_rad, diag%ls_qcl_rad, &
  model_levels)

! LS cloud ice water mixing ratio
CALL output_diag_3d(309+i_off, ip_all_calls, &
  diag%l_ls_qcf_rad, diag%ls_qcf_rad, &
  model_levels)

! Convective cloud liquid water mixing ratio
CALL output_diag_3d(310+i_off, ip_all_calls, &
  diag%l_cc_qcl_rad, diag%cc_qcl_rad, &
  model_levels)

! Convective cloud ice water mixing ratio
CALL output_diag_3d(311+i_off, ip_all_calls, &
  diag%l_cc_qcf_rad, diag%cc_qcf_rad, &
  model_levels)


! Cloud water paths (2D)

! Large scale liquid water path
CALL output_diag(391+i_off, ip_all_calls, &
  diag%l_ls_qcl_rad_path, diag%ls_qcl_rad_path)

! Large scale ice water path
CALL output_diag(392+i_off, ip_all_calls, &
  diag%l_ls_qcf_rad_path, diag%ls_qcf_rad_path)

! Convective liquid water path
CALL output_diag(393+i_off, ip_all_calls, &
  diag%l_cc_qcl_rad_path, diag%cc_qcl_rad_path)

! Convective ice water path
CALL output_diag(394+i_off, ip_all_calls, &
  diag%l_cc_qcf_rad_path, diag%cc_qcf_rad_path)


! Cloud amounts

! LS liquid cloud fraction of gridbox seen by radiation
CALL output_diag_3d(312+i_off, ip_all_calls, &
  diag%l_ls_cl_rad, diag%ls_cl_rad, &
  model_levels)

! LS ice cloud fraction of gridbox seen by radiation
CALL output_diag_3d(313+i_off, ip_all_calls, &
  diag%l_ls_cf_rad, diag%ls_cf_rad, &
  model_levels)

! CONV liquid cloud fraction of gridbox seen by radiation
CALL output_diag_3d(314+i_off, ip_all_calls, &
  diag%l_cc_cl_rad, diag%cc_cl_rad, &
  model_levels)

! CONV ice cloud fraction of gridbox seen by radiation
CALL output_diag_3d(315+i_off, ip_all_calls, &
  diag%l_cc_cf_rad, diag%cc_cf_rad, &
  model_levels)


! Cloud effective dimension diagnostics

! LS liquid effective dimension
CALL output_diag_3d(397+i_off, ip_all_calls, &
  diag%l_ls_del_rad, diag%ls_del_rad, &
  model_levels)

! LS ice effective dimension
CALL output_diag_3d(398+i_off, ip_all_calls, &
  diag%l_ls_def_rad, diag%ls_def_rad, &
  model_levels)

! CNV liquid effective dimension
CALL output_diag_3d(399+i_off, ip_all_calls, &
  diag%l_cc_del_rad, diag%cc_del_rad, &
  model_levels)

! CNV ice effective dimension
CALL output_diag_3d(400+i_off, ip_all_calls, &
  diag%l_cc_def_rad, diag%cc_def_rad, &
  model_levels)

END SUBROUTINE diagnostics_lw


SUBROUTINE output_diag(item,availability,flag,field)
  ! Output single level fields
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field(:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
    ! DEPENDS ON: copydiag
    CALL copydiag (STASHwork(si(item,sect,im_index)), field,                   &
      row_length, rows, 0, 0, 0, 0, at_extremity,                              &
      atmos_im, sect, item, icode, cmessage)
    IF (icode /= 0) THEN
      WRITE(cmessage,'(A,I3,A)') 'Error in copydiag (item ',item,')'
      CALL ereport(RoutineName,icode,cmessage)
    END IF
  END IF
END SUBROUTINE output_diag


SUBROUTINE output_diag_ps(item,availability,flag,field,pseudo_levels)
  ! Output single level fields on all pseudo levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, pseudo_levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field(:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
    DO pslevel = 1, pseudo_levels
      CALL copydiag(STASHwork(si(item,sect,im_index)                           &
        +(row_length*rows*(pslevel-1))), field(:,:,pslevel),                   &
        row_length, rows, 0, 0, 0, 0, at_extremity,                            &
        atmos_im, sect, item, icode, cmessage)
      IF (icode /= 0) THEN
        WRITE(cmessage,'(A,I3,A)') 'Error in copydiag (item ',item,')'
        CALL ereport(RoutineName,icode,cmessage)
      END IF
    END DO
  END IF
END SUBROUTINE output_diag_ps


SUBROUTINE output_diag_pl(item,availability,flag,field,pseudo_levels)
  ! Output single level fields on list of pseudo levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, pseudo_levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field(:,:,:)

  LOGICAL :: l_output
  LOGICAL :: pll(pseudo_levels) ! pseudolevel list
  INTEGER :: pslevel_out        ! index for pseudolevels sent to STASH

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
    ! DEPENDS ON: set_pseudo_list
    CALL set_pseudo_list(pseudo_levels, len_stlist,                            &
         stlist(1,stindex(1,item,sect,im_index)), pll,                         &
         stash_pseudo_levels, num_stash_pseudo, icode, cmessage)
    IF (icode /= 0) THEN
      WRITE(cmessage,'(A,I3,A)') 'Error in set_pseudo_list (item ',item,')'
      CALL ereport(RoutineName,icode,cmessage)
    END IF
    pslevel_out=0
    DO pslevel = 1, pseudo_levels
      IF (pll(pslevel)) THEN
        pslevel_out=pslevel_out+1
        CALL copydiag(stashwork(si(item,sect,im_index)                         &
          +(pslevel_out-1)*row_length*rows), field(:,:,pslevel),               &
          row_length, rows, 0, 0, 0, 0, at_extremity,                          &
          atmos_im, sect, item, icode, cmessage)
        IF (icode /= 0) THEN
          WRITE(cmessage,'(A,I3,A)') 'Error in copydiag (item ',item,')'
          CALL ereport(RoutineName,icode,cmessage)
        END IF
      END IF
    END DO
  END IF
END SUBROUTINE output_diag_pl


SUBROUTINE output_diag_3d(item,availability,flag,field,levels)
  ! Output fields on levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field(:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
    ! DEPENDS ON: copydiag_3d
    CALL copydiag_3d(STASHwork(si(item,sect,im_index)), field,                 &
      row_length, rows, levels, 0, 0, 0, 0, at_extremity,                      &
      stlist(1,stindex(1,item,sect,im_index)), len_stlist,                     &
      stash_levels, num_stash_levels+1,                                        &
      atmos_im, sect, item, icode, cmessage)
    IF (icode /= 0) THEN
      WRITE(cmessage,'(A,I3,A)') 'Error in copydiag_3d (item ',item,')'
      CALL ereport(RoutineName,icode,cmessage)
    END IF
  END IF
END SUBROUTINE output_diag_3d


SUBROUTINE output_diag_03d(item,availability,flag,field,levels)
  ! Output fields on levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field(:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
    ! DEPENDS ON: copydiag_03d
    CALL copydiag_03d(STASHwork(si(item,sect,im_index)), field,                &
      row_length, rows, levels, 0, 0, 0, 0, at_extremity,                      &
      stlist(1,stindex(1,item,sect,im_index)), len_stlist,                     &
      stash_levels, num_stash_levels+1,                                        &
      atmos_im, sect, item, icode, cmessage)
    IF (icode /= 0) THEN
      WRITE(cmessage,'(A,I3,A)') 'Error in copydiag_03d (item ',item,')'
      CALL ereport(RoutineName,icode,cmessage)
    END IF
  END IF
END SUBROUTINE output_diag_03d


SUBROUTINE output_diag_4d(item,availability,flag,field,levels,pseudo_levels)
  ! Output fields on levels and all pseudo levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, levels, pseudo_levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field(:,:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
    DO pslevel = 1, pseudo_levels
      CALL copydiag_3d(STASHwork(si(item,sect,im_index)                        &
        +(row_length*rows*levels*(pslevel-1))),                                &
        field(:,:,:,pslevel),                                                  &
        row_length, rows, levels, 0, 0, 0, 0, at_extremity,                    &
        stlist(1,stindex(1,item,sect,im_index)), len_stlist,                   &
        stash_levels, num_stash_levels+1,                                      &
        atmos_im, sect, item, icode, cmessage)
      IF (icode /= 0) THEN
        WRITE(cmessage,'(A,I3,A)') 'Error in copydiag_3d (item ',item,')'
        CALL ereport(RoutineName,icode,cmessage)
      END IF
    END DO
  END IF
END SUBROUTINE output_diag_4d


SUBROUTINE output_diag_04d(item,availability,flag,field,levels,pseudo_levels)
  ! Output fields on levels and all pseudo levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, levels, pseudo_levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field(:,:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
    DO pslevel = 1, pseudo_levels
      CALL copydiag_03d(STASHwork(si(item,sect,im_index)                       &
        +(row_length*rows*(levels+1)*(pslevel-1))),                            &
        field(:,:,:,pslevel),                                                  &
        row_length, rows, levels, 0, 0, 0, 0, at_extremity,                    &
        stlist(1,stindex(1,item,sect,im_index)), len_stlist,                   &
        stash_levels, num_stash_levels+1,                                      &
        atmos_im, sect, item, icode, cmessage)
      IF (icode /= 0) THEN
        WRITE(cmessage,'(A,I3,A)') 'Error in copydiag_03d (item ',item,')'
        CALL ereport(RoutineName,icode,cmessage)
      END IF
    END DO
  END IF
END SUBROUTINE output_diag_04d


SUBROUTINE output_diag_3d_inc(item,availability,flag,field1,field2,levels)
  ! Output increment between input fields on levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field1(:,:,:), field2(:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(levels,rows,row_length,work_3d,field1,field2)
    DO k = 1, levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = field1(i,j,k) - field2(i,j,k)
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO
    CALL output_diag_3d(item,availability,flag,work_3d,levels)
  END IF
END SUBROUTINE output_diag_3d_inc


SUBROUTINE output_diag_3d_pinc(item,availability,flag,field1,field2,levels)
  ! Output positive increment between input fields on levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field1(:,:,:), field2(:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(levels,rows,row_length,work_3d,field1,field2)
    DO k = 1, levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = MAX(0.0, field1(i,j,k) - field2(i,j,k))
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO
    CALL output_diag_3d(item,availability,flag,work_3d,levels)
  END IF
END SUBROUTINE output_diag_3d_pinc


SUBROUTINE output_diag_3d_ninc(item,availability,flag,field1,field2,levels)
  ! Output negative increment between input fields on levels
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: item, availability, levels
  LOGICAL, INTENT(IN) :: flag
  REAL,    INTENT(IN) :: field1(:,:,:), field2(:,:,:)
  LOGICAL :: l_output

  CALL check_availability(item,availability,flag,l_output)
  IF (l_output) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(levels,rows,row_length,work_3d,field1,field2)
    DO k = 1, levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = MIN(0.0, field1(i,j,k) - field2(i,j,k))
        END DO ! i
      END DO   ! j
    END DO     ! k
!$OMP END PARALLEL DO
    CALL output_diag_3d(item,availability,flag,work_3d,levels)
  END IF
END SUBROUTINE output_diag_3d_ninc


SUBROUTINE check_availability(item,availability,flag,l_output)
  ! Check if diagnostic should be output for this radiation call
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: item, availability
  LOGICAL, INTENT(IN)  :: flag
  LOGICAL, INTENT(OUT) :: l_output

  l_output = .FALSE.
  SELECT CASE (availability)
  CASE (ip_main_call)
    IF (sf_calc(item,sect) .AND. (j_rad == 1) .AND. flag) THEN
      l_output = .TRUE.
    END IF
  CASE (ip_main_and_forcing_call)
    IF (sf_calc(item,sect) .AND. l_rad_perturb .AND. (j_rad == 1)) THEN
      l_output = .TRUE.
    END IF
    IF (flag .AND. .NOT.l_rad_perturb) THEN
      l_output = .TRUE.
    END IF
  CASE (ip_all_calls)
    IF (flag .AND. .NOT.(l_rad_perturb .AND. (j_rad == 1))) THEN
      l_output = .TRUE.
    END IF
  CASE DEFAULT
    icode=1
    WRITE(cmessage,'(A,I3)') 'The wrong availability code has been set '//     &
      'for item ',item+i_off
    CALL ereport(RoutineName,icode,cmessage)
  END SELECT
END SUBROUTINE check_availability


END SUBROUTINE diagnostics_rad
END MODULE diagnostics_rad_mod
