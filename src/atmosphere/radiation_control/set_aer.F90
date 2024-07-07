! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the aerosol fields and calculate aerosol optical depths
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_aer_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_AER_MOD'
CONTAINS

SUBROUTINE set_aer(control, dimen, spectrum, atm, cld, aer,                    &
! Model grid
  list,                                                                        &
! Properties of CLASSIC aerosols
  l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,                      &
  bl_depth, n_levels_bl, l_murk_rad, aero_meso,                                &
  l_dust, l_use_dust, dust_dim1, dust_dim2,                                    &
  dust_1, dust_2, dust_3, dust_4, dust_5, dust_6,                              &
  l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic,                      &
  l_sulpc_so2, l_use_sulpc_direct, sulp_dim1, sulp_dim2,                       &
  accum_sulphate, aitken_sulphate,                                             &
  l_use_seasalt_direct, salt_dim1, salt_dim2, salt_dim3,                       &
  sea_salt_film, sea_salt_jet,                                                 &
  l_soot, l_use_soot_direct, soot_dim1, soot_dim2,                             &
  fresh_soot, aged_soot,                                                       &
  l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2,                       &
  fresh_bmass, aged_bmass,                                                     &
  l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2,                             &
  fresh_ocff, aged_ocff,                                                       &
  l_nitrate, l_use_nitrate_direct, nitrate_dim1, nitrate_dim2,                 &
  accum_nitrate,                                                               &
  n_arcl_species, n_arcl_compnts, i_arcl_compnts,                              &
  l_use_arcl, arcl_dim1, arcl_dim2, arcl,                                      &
  land, lying_snow, pstar,                                                     &
  p_layer_boundaries, trindx, alat, previous_time,                             &
! Properties of UKCA aerosols
  l_ukca_radaer, l_glomap_clim_radaer, ukca_radaer, ukca_dim1, ukca_dim2,      &
  ukca_mmr, ukca_cvl, ukca_dry, ukca_wet,                                      &
  ukca_rho, ukca_vol, ukca_wtv, ukca_nbr,                                      &
! Properties of EasyAerosol
  l_easyaerosol_rad, easyaerosol_rad,                                          &
! Diagnostics
  diag, row_list, col_list,                                                    &
! Dimensions of arrays
  nd_field, n_ukca_cpnt, n_ukca_mode)

USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts
USE rad_pcf
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_atm,      ONLY: StrAtm
USE def_cld,      ONLY: StrCld
USE def_aer,      ONLY: StrAer, allocate_aer, allocate_aer_prsc
USE def_diag,     ONLY: StrDiag
USE ukca_radaer_struct_mod, ONLY:                                              &
    ip_ukca_mode_aitken, ip_ukca_mode_accum, ip_ukca_mode_coarse,              &
    ukca_radaer_struct
USE dust_parameters_mod, ONLY: l_twobin_dust
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE ereport_mod,  ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE aodtype_mod,  ONLY: ip_type_allaod, ip_type_sulphate, ip_type_dust,        &
                        ip_type_seasalt, ip_type_soot, ip_type_biomass,        &
                        ip_type_biogenic, ip_type_ocff, ip_type_delta,         &
                        ip_type_nitrate, ip_type_twobdust
USE rad_input_mod, ONLY: l_extra_top
USE nlsizes_namelist_mod, ONLY: model_levels

USE ukca_radaer_band_average_mod, ONLY: ukca_radaer_band_average
USE ukca_radaer_compute_aod_mod, ONLY: ukca_radaer_compute_aod
USE ukca_radaer_3d_diags_mod, ONLY: ukca_radaer_3d_diags
USE ukca_radaer_prepare_mod, ONLY: ukca_radaer_prepare
USE ukca_radaer_set_aerosol_field_mod, ONLY: ukca_radaer_set_aerosol_field
USE def_easyaerosol, ONLY: t_easyaerosol_rad
USE easyaerosol_mod, ONLY: set_easyaerosol_field

USE compute_all_aod_mod, ONLY: compute_all_aod
USE r2_set_aerosol_field_mod, ONLY: r2_set_aerosol_field
USE set_moist_aerosol_properties_mod, ONLY: set_moist_aerosol_properties
IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(IN)  :: control

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Atmospheric properties:
TYPE(StrAtm),       INTENT(IN)  :: atm

! Cloud properties:
TYPE(StrCld),       INTENT(IN)  :: cld

! Aerosol properties:
TYPE(StrAer),       INTENT(OUT) :: aer

INTEGER, INTENT(IN) ::                                                         &
  nd_field,                                                                    &
!   Field size in calling program
  n_levels_bl,                                                                 &
!   Number of layers occupied by boundary-layer aerosol
!   if l_clim_aero_hgt is false.
  n_ukca_mode,                                                                 &
!   Number of aerosol modes in UKCA_RADAER
  n_ukca_cpnt
!   Number of aerosol components in UKCA_RADAER

INTEGER, INTENT(IN) ::                                                         &
  list(nd_field)
!   List of points in full model array on which radiation is calculated

! CLASSIC aerosols:
LOGICAL, INTENT(IN) ::                                                         &
  l_climat_aerosol,                                                            &
!   Flag for climatological aerosol
  l_clim_aero_hgt,                                                             &
!   Flag to use the depth of the boundary layer to set
!   the climatological aerosol
  l_hadgem1_clim_aero,                                                         &
!   Flag to use HadGEM1 setting for climatological aerosols
  l_murk_rad,                                                                  &
!   Flag for mesoscale model aerosol
  l_sulpc_so2,                                                                 &
!   Sulphur cycle available for effecting fluxes or diagnostics
  l_use_sulpc_direct,                                                          &
!   Flag to use sulphur cycle for direct effect
  l_dust,                                                                      &
!   Dust is available for effecting fluxes or diagnostics
  l_use_dust,                                                                  &
!   Flag to use direct rad effect of mineral dust
  l_use_biogenic,                                                              &
!   Flag to use biogenic for direct effect
  l_soot,                                                                      &
!   Soot is available for effecting fluxes or diagnostics
  l_use_soot_direct,                                                           &
!   Use direct rad. effect of soot aerosol
  l_biomass,                                                                   &
!   Biomass is available for effecting fluxes or diagnostics
  l_use_bmass_direct,                                                          &
!   Flag to use direct rad. effect of biomass smoke
  l_ocff,                                                                      &
!   OCFF is available for effecting fluxes or diagnostics
  l_use_ocff_direct,                                                           &
!   Flag to use direct rad. effect of ocff
  l_use_seasalt_direct,                                                        &
!   Flag to use sea-salt for direct effect
  l_nitrate,                                                                   &
!   Nitrate is available for effecting fluxes or diagnostics
  l_use_nitrate_direct,                                                        &
!   Flag to use nitrate for direct effect
  l_easyaerosol_rad
!   Flag to use EasyAerosol for direct effect

INTEGER, INTENT(IN) ::                                                         &
  sulp_dim1,sulp_dim2,                                                         &
!   Dimensions for _sulphate arrays, (P_FIELD,P_LEVELS or 1,1)
  dust_dim1, dust_dim2,                                                        &
!   Dimensions for mineral dust arrays (p_field,p_levels or 1,1)
  biogenic_dim1, biogenic_dim2,                                                &
!   dimensions for biogenic array passed down to
!   r2_set_aerosol_field if direct effect required.
  soot_dim1, soot_dim2,                                                        &
!   Dimensions for soot arrays (P_FIELD,P_LEVELS or 1,1)
  bmass_dim1, bmass_dim2,                                                      &
!   Dimensions for biomass arrays (P_FIELD,P_LEVELS or 1,1)
  ocff_dim1, ocff_dim2,                                                        &
!   Dimensions for ocff arrays (P_FIELD,P_LEVELS or 1,1)
  nitrate_dim1, nitrate_dim2,                                                  &
!   Dimensions for nitrate arrays (P_FIELD,P_LEVELS or 1,1)
  salt_dim1, salt_dim2, salt_dim3
!   dimensions for salt arrays on input (salt_dim1*salt_dim2=p_field
!   and salt_dim3=p_levels, or else 1,1,1)

REAL, INTENT(IN) ::                                                            &
  accum_sulphate(sulp_dim1, sulp_dim2),                                        &
!   Mass mixing ratio of accumulation mode aerosol
  aitken_sulphate(sulp_dim1, sulp_dim2),                                       &
!   Mass mixing ratio of aitken mode aerosol
  dust_1(dust_dim1, dust_dim2),                                                &
!   Mass mixing ratio of div1 dust
  dust_2(dust_dim1, dust_dim2),                                                &
!   Mass mixing ratio of div2 dust
  dust_3(dust_dim1, dust_dim2),                                                &
!   Mass mixing ratio of div3 dust
  dust_4(dust_dim1, dust_dim2),                                                &
!   Mass mixing ratio of div4 dust
  dust_5(dust_dim1, dust_dim2),                                                &
!   Mass mixing ratio of div5 dust
  dust_6(dust_dim1, dust_dim2),                                                &
!   Mass mixing ratio of div6 dust
  biogenic(biogenic_dim1, biogenic_dim2),                                      &
!   Mixing ratios of biogenic aerosol
  sea_salt_film(salt_dim1, salt_dim2, salt_dim3),                              &
!   Number concentration of film-mode sea-salt aerosol
  sea_salt_jet(salt_dim1, salt_dim2, salt_dim3),                               &
!   Number concentration of jet-mode sea-salt aerosol
  fresh_soot(soot_dim1, soot_dim2),                                            &
!   Soot mixing ratios
  aged_soot(soot_dim1, soot_dim2),                                             &
!   Soot mixing ratios
  fresh_bmass(bmass_dim1, bmass_dim2),                                         &
!   Mass mixing ratio of fresh biomass smoke
  aged_bmass(bmass_dim1, bmass_dim2),                                          &
!   Mass mixing ratio of aged biomass smoke
  fresh_ocff(ocff_dim1, ocff_dim2),                                            &
!   Mass mixing ratio of fresh fossil-fuel organic carbon aer
  aged_ocff(ocff_dim1, ocff_dim2),                                             &
!   Mass mixing ratio of aged fossil-fuel organic carbon aer
  accum_nitrate(nitrate_dim1, nitrate_dim2),                                   &
!   Mass mixing ratio of accumulation nitrate aerosol
  bl_depth(nd_field),                                                          &
!   Depth of the boundary layer
  aero_meso(nd_field, model_levels)
!   Mixing ratio of 'urban' aerosol of mesoscale model


! Aerosol climatology for NWP
INTEGER, INTENT(IN) ::                                                         &
  n_arcl_species,                                                              &
!   Number of requested species within the climatology
  n_arcl_compnts,                                                              &
!   Corresponding number of requested components
  i_arcl_compnts(npd_arcl_compnts),                                            &
!   Index of each component
  arcl_dim1, arcl_dim2
!   Dimensions

LOGICAL, INTENT(IN) ::                                                         &
  l_use_arcl(npd_arcl_species)
!   Model switch for each species

REAL, INTENT(IN) ::                                                            &
  arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)
!   Mass-mixing ratios

! EasyAerosol climatology
TYPE (t_easyaerosol_rad), INTENT(IN) :: easyaerosol_rad

! Land mask
LOGICAL, INTENT(IN) :: land(nd_field)

! Mass loading of lying snow
REAL, INTENT(IN) :: lying_snow(nd_field)

! Surface pressures
REAL, INTENT(IN) :: pstar(nd_field)

! Pressure at boundaries of layers
REAL, INTENT(IN) :: p_layer_boundaries(nd_field, 0:model_levels)

! The layer boundary of the tropopause
INTEGER, INTENT(IN) :: trindx(nd_field)

! Latitude in degrees
REAL, INTENT(IN) :: alat(dimen%nd_profile)

! Model time at beginning of timestep
INTEGER, INTENT(IN) :: previous_time(7)


! Fields for UKCA aerosols

! Model switch
LOGICAL, INTENT(IN) :: l_ukca_radaer
LOGICAL, INTENT(IN) :: l_glomap_clim_radaer

! UKCA_RADAER structure
TYPE (ukca_radaer_struct), INTENT(IN) :: ukca_radaer

! Dimensions
INTEGER, INTENT(IN) :: ukca_dim1
INTEGER, INTENT(IN) :: ukca_dim2

! Component mass mixing ratios and volumes
REAL, INTENT(IN) :: ukca_mmr(ukca_dim1, ukca_dim2, n_ukca_cpnt)
REAL, INTENT(IN) :: ukca_cvl(ukca_dim1, ukca_dim2, n_ukca_cpnt)

! Dry and wet modal diameters
REAL, INTENT(IN) :: ukca_dry(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL, INTENT(IN) :: ukca_wet(ukca_dim1, ukca_dim2, n_ukca_mode)

! Modal densities, volumes, volume of water, and number conc
REAL, INTENT(IN) :: ukca_rho(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL, INTENT(IN) :: ukca_vol(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL, INTENT(IN) :: ukca_wtv(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL, INTENT(IN) :: ukca_nbr(ukca_dim1, ukca_dim2, n_ukca_mode)

! Diagnostics:
TYPE (StrDiag), INTENT(INOUT) :: diag

INTEGER, INTENT(IN) :: row_list(dimen%nd_profile)
!                          list of row indices of lit points
INTEGER, INTENT(IN) :: col_list(dimen%nd_profile)
!                          list of column indices of lit points


! Local arguments:

! UKCA modal optical depth diagnostics: full column
REAL ::                                                                        &
  aod_ukca_ait_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  aod_ukca_acc_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  aod_ukca_cor_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  aod_ukca_ait_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  aod_ukca_acc_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  aod_ukca_cor_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel)

! UKCA modal optical depth diagnostics: stratosphere
REAL ::                                                                        &
  sod_ukca_ait_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  sod_ukca_acc_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  sod_ukca_cor_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  sod_ukca_ait_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  sod_ukca_acc_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel),               &
  sod_ukca_cor_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel)

! UKCA modal absorption optical depth diagnostics: full column
REAL ::                                                                        &
  aaod_ukca_ait_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),              &
  aaod_ukca_acc_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),              &
  aaod_ukca_cor_sol(dimen%nd_profile, spectrum%dim%nd_aod_wavel),              &
  aaod_ukca_ait_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel),              &
  aaod_ukca_acc_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel),              &
  aaod_ukca_cor_ins(dimen%nd_profile, spectrum%dim%nd_aod_wavel)

! UKCA 3D extinction diagnostic
REAL ::                                                                        &
  ukca_aerosol_ext(dimen%nd_profile, dimen%nd_layer)

! UKCA 3D absorption diagnostic
REAL ::                                                                        &
  ukca_aerosol_abs(dimen%nd_profile, dimen%nd_layer)

! UKCA 3D scattering diagnostic
REAL ::                                                                        &
  ukca_aerosol_sca(dimen%nd_profile, dimen%nd_layer)

! UKCA 3D asymmetry * scattering paraemter diagnostic
REAL ::                                                                        &
  ukca_aerosol_gsca(dimen%nd_profile, dimen%nd_layer)

INTEGER :: i_pointer_water
!   Pointer to water vapour

INTEGER :: i_humidity_pointer(dimen%nd_profile, dimen%nd_layer)
!   Pointer to look-up table of humidities for aerosols

INTEGER :: i, j, l, k
!   Loop variables 

INTEGER :: iwv_ext, iwv_abs, iwv_sca, iwv_gsca
!   Indicies for 3D aerosol-radiation arrays

LOGICAL :: l_moist_aerosol
!   Flag for moist aerosol

REAL :: delta_humidity
!   Increment in look-up table for hum.

REAL :: ukca_modal_number(dimen%nd_profile, dimen%nd_layer, n_ukca_mode)
!   UKCA aerosol modal number concentrations

REAL :: ukca_mix_ratio(dimen%nd_profile, dimen%nd_layer, n_ukca_cpnt)
REAL :: ukca_comp_vol(dimen%nd_profile, dimen%nd_layer, n_ukca_cpnt)
REAL :: ukca_dry_diam(dimen%nd_profile, dimen%nd_layer, n_ukca_mode)
REAL :: ukca_wet_diam(dimen%nd_profile, dimen%nd_layer, n_ukca_mode)
REAL :: ukca_modal_rho(dimen%nd_profile, dimen%nd_layer, n_ukca_mode)
REAL :: ukca_modal_vol(dimen%nd_profile, dimen%nd_layer, n_ukca_mode)
REAL :: ukca_modal_wtv(dimen%nd_profile, dimen%nd_layer, n_ukca_mode)
REAL :: ukca_modal_nbr(dimen%nd_profile, dimen%nd_layer, n_ukca_mode)
!   UKCA aerosol component mass-mixing ratios, modal dry and
!   wet diameters, modal densities, volumes, and volumes of water.

REAL :: rm
INTEGER :: trindxrad(dimen%nd_profile)
!   Level of tropopause on radiation scheme grid for UKCA-MODE aerosols

LOGICAL, PARAMETER :: soluble_wanted = .TRUE.
LOGICAL, PARAMETER :: soluble_unwanted = .FALSE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                       :: ierr = i_normal
CHARACTER (LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'SET_AER'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL allocate_aer(aer, dimen, spectrum)
CALL allocate_aer_prsc(aer, dimen, spectrum)

! Set the mixing ratios of aerosols.
IF (control%l_aerosol .OR. control%l_aerosol_ccn) THEN
  CALL r2_set_aerosol_field(ierr,                                              &
    atm%n_profile, model_levels, atm%n_layer, spectrum%aerosol%n_aerosol,      &
    spectrum%aerosol%n_aerosol_mr, spectrum%aerosol%type_aerosol,              &
    list, row_list, col_list, l_extra_top,                                     &
    l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,                    &
    bl_depth, atm%t, n_levels_bl, l_murk_rad, aero_meso,                       &
    l_dust, l_use_dust, dust_dim1, dust_dim2,                                  &
    dust_1, dust_2, dust_3, dust_4, dust_5, dust_6,                            &
    l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic,                    &
    l_sulpc_so2, l_use_sulpc_direct, sulp_dim1, sulp_dim2,                     &
    accum_sulphate, aitken_sulphate,                                           &
    l_use_seasalt_direct, salt_dim1, salt_dim2, salt_dim3,                     &
    sea_salt_film, sea_salt_jet, atm%p,                                        &
    l_soot, l_use_soot_direct, soot_dim1, soot_dim2,                           &
    fresh_soot, aged_soot,                                                     &
    l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2,                     &
    fresh_bmass, aged_bmass,                                                   &
    l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2,                           &
    fresh_ocff, aged_ocff,                                                     &
    l_nitrate,l_use_nitrate_direct, nitrate_dim1, nitrate_dim2,                &
    accum_nitrate,                                                             &
    n_arcl_species, n_arcl_compnts, i_arcl_compnts,                            &
    l_use_arcl, arcl_dim1, arcl_dim2, arcl,                                    &
    land, lying_snow, pstar,                                                   &
    p_layer_boundaries, trindx, alat, previous_time,                           &
    aer%mix_ratio, aer%mr_source, aer%mr_type_index,                           &
    nd_field, dimen%nd_profile, dimen%nd_layer,                                &
    spectrum%dim%nd_aerosol_species,                                           &
    spectrum%dim%nd_aerosol_mr, 1)
  IF (ierr /= i_normal) THEN
    cmessage = 'Error following call to r2_set_aerosol_field'
    GO TO 9999
  END IF
END IF

! Set the mixing ratios of UKCA aerosols
IF (l_ukca_radaer .OR. l_glomap_clim_radaer) THEN
  CALL ukca_radaer_set_aerosol_field(                                          &
    list, l_extra_top, atm%n_layer, atm%n_profile,                             &
    ukca_dim1, ukca_dim2,                                                      &
    ukca_mmr, ukca_cvl,                                                        &
    ukca_dry, ukca_wet, ukca_rho, ukca_vol, ukca_wtv,                          &
    ukca_nbr, trindx,                                                          &
    ukca_mix_ratio, ukca_comp_vol,                                             &
    ukca_dry_diam, ukca_wet_diam, ukca_modal_rho,                              &
    ukca_modal_vol, ukca_modal_wtv,                                            &
    ukca_modal_nbr, trindxrad,                                                 &
    nd_field, dimen%nd_profile, dimen%nd_layer,                                &
    n_ukca_cpnt, n_ukca_mode)
  aer%n_mode = n_ukca_mode
ELSE
  aer%n_mode = 0
END IF

! Calculate the clear-sky mean relative humidity for moist aerosols.
l_moist_aerosol=.FALSE.
DO j=1, spectrum%aerosol%n_aerosol
  l_moist_aerosol=l_moist_aerosol.OR.                                          &
    (spectrum%aerosol%i_aerosol_parm(j) == ip_aerosol_param_moist).OR.         &
    (spectrum%aerosol%i_aerosol_parm(j) == ip_aerosol_param_phf_moist)
END DO

IF (l_moist_aerosol) THEN
  i_pointer_water=MAX(spectrum%cont%index_water, 1)
  CALL set_moist_aerosol_properties(ierr,                                      &
    atm%n_profile, atm%n_layer,                                                &
    spectrum%aerosol%n_aerosol, spectrum%aerosol%i_aerosol_parm,               &
    spectrum%aerosol%nhumidity, control%l_mixing_ratio,                        &
    atm%gas_mix_ratio(1, 1, i_pointer_water), atm%t, atm%p, cld%w_cloud,       &
    delta_humidity, aer%mean_rel_humidity, i_humidity_pointer,                 &
    dimen%nd_profile, dimen%nd_layer, dimen%id_cloud_top,                      &
    spectrum%dim%nd_aerosol_species)
END IF

! For the CLASSIC aerosols, all AOD calculations are dealt with via a
! wrapper subroutine
IF (control%l_aerosol .OR. control%l_aerosol_ccn) THEN
  CALL compute_all_aod(diag, row_list, col_list,                               &
    spectrum%aerosol%n_aerosol,      spectrum%aerosol%n_aerosol_mr,            &
    spectrum%dim%nd_aerosol_species, spectrum%dim%nd_aerosol_mr,               &
    spectrum%aerosol%n_aod_wavel,    spectrum%dim%nd_aod_wavel,                &
    spectrum%aerosol%aod_wavel,                                                &
    atm%n_profile,  dimen%nd_profile,                                          &
    atm%n_layer, 1, dimen%nd_layer,                                            &
    spectrum%dim%nd_humidity,        spectrum%dim%nd_humidity,                 &
    spectrum%aerosol%type_aerosol,   spectrum%aerosol%i_aod_type,              &
    aer%mix_ratio, atm%mass, atm%density,                                      &
    aer%mr_source, aer%mr_type_index,                                          &
    spectrum%aerosol%i_aerosol_parm, spectrum%aerosol%aod_abs,                 &
    spectrum%aerosol%aod_scat,                                                 &
    i_humidity_pointer, aer%mean_rel_humidity,                                 &
    spectrum%aerosol%humidities,     delta_humidity,    l_use_arcl)
END IF

! For UKCA aerosols, call RADAER routines to calculate optical properties.
IF (l_ukca_radaer .OR. l_glomap_clim_radaer) THEN
  ! Calculations shared by ukca_radaer_band_average() and
  ! ukca_radaer_compute_aod().
  CALL ukca_radaer_prepare(                                                    &
    ! Actual array dimensions
    atm%n_profile, atm%n_layer, n_ukca_mode, n_ukca_cpnt,                      &
    ! UKCA_RADAER structure
    ukca_radaer,                                                               &
    ! Component mass-mixing ratios
    ukca_mix_ratio,                                                            &
    ! Modal mass-mixing ratios
    aer%mode_mix_ratio,                                                        &
    ! Input modal number concentrations
    ukca_modal_nbr,                                                            &
    ! Output modal number concentrations
    ukca_modal_number,                                                         &
    ! Pressure and temperature
    atm%p, atm%t,                                                              &
    ! Fixed array dimensions
    dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode)

  ! Compute the band-averaged optical properties for UKCA aerosols
  CALL ukca_radaer_band_average(                                               &
    ! Spectral information
    spectrum%basic%n_band, control%isolir, spectrum%basic%l_present(14),       &
    spectrum%basic%n_band_exclude, spectrum%basic%index_exclude,               &
    ! Actual array dimensions
    atm%n_profile, atm%n_layer, n_ukca_mode, n_ukca_cpnt,                      &
    ! UKCA_RADAER structure
    ukca_radaer,                                                               &
    ! Modal mass-mixing ratios
    aer%mode_mix_ratio,                                                        &
    ! Modal number concentrations
    ukca_modal_number,                                                         &
    ! Modal diameters from UKCA module
    ukca_dry_diam, ukca_wet_diam,                                              &
    ! Other inputs from UKCA module
    ukca_comp_vol, ukca_modal_vol, ukca_modal_rho,                             &
    ukca_modal_wtv,                                                            &
    ! Model level of the tropopause
    trindxrad,                                                                 &
    ! Band-averaged optical properties (outputs)
    aer%mode_absorption, aer%mode_scattering, aer%mode_asymmetry,              &
    ! Fixed array dimensions
    dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode,                   &
    spectrum%dim%nd_band, spectrum%dim%nd_exclude) 

  ! Calculate AOD diagnostics for Aitken soluble mode if requested
  IF (diag%l_aod_ukca_ait_sol .OR. diag%l_sod_ukca_ait_sol .OR.                &
       diag%l_aaod_ukca_ait_sol ) THEN

    CALL ukca_radaer_compute_aod(                                              &
      ! Actual array dimension
      atm%n_profile, atm%n_layer, n_ukca_mode,                                 &
      n_ukca_cpnt, spectrum%aerosol%n_aod_wavel,                               &
      ! UKCA_RADAER structure
      ukca_radaer,                                                             &
      ! Modal diameters from UKCA module
      ukca_dry_diam, ukca_wet_diam,                                            &
      ! Mass thickness of layers
      atm%mass,                                                                &
      ! Component volumes
      ukca_comp_vol,                                                           &
      ! Modal volumes, densities, and water content
      ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,                          &
      ! Modal mass-mixing ratios
      aer%mode_mix_ratio,                                                      &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Type selection
      ip_ukca_mode_aitken, soluble_wanted,                                     &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Modal extinction aerosol opt depth (output)
      aod_ukca_ait_sol, sod_ukca_ait_sol,                                      &
      ! Modal absorption aerosol opt depth (output)
      aaod_ukca_ait_sol,                                                       &
      ! Fixed array dimensions
      dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode,                 &
      spectrum%dim%nd_aod_wavel)

    ! Put requested AOD diagnostics in diag structure

    IF (diag%l_aod_ukca_ait_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aod_ukca_ait_sol(col_list(l), row_list(l), i)                   &
            = aod_ukca_ait_sol(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_sod_ukca_ait_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%sod_ukca_ait_sol(col_list(l), row_list(l), i)                   &
            = sod_ukca_ait_sol(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_aaod_ukca_ait_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aaod_ukca_ait_sol(col_list(l), row_list(l), i)                  &
            = aaod_ukca_ait_sol(l, i)
        END DO
      END DO
    END IF
  END IF 

  ! Calculate AOD diagnostics for accumulation soluble mode if requested
  IF (diag%l_aod_ukca_acc_sol .OR. diag%l_sod_ukca_acc_sol .OR.                &
       diag%l_aaod_ukca_acc_sol ) THEN

    CALL ukca_radaer_compute_aod(                                              &
      ! Actual array dimension
      atm%n_profile, atm%n_layer, n_ukca_mode,                                 &
      n_ukca_cpnt, spectrum%aerosol%n_aod_wavel,                               &
      ! UKCA_RADAER structure
      ukca_radaer,                                                             &
      ! Modal diameters from UKCA module
      ukca_dry_diam, ukca_wet_diam,                                            &
      ! Mass thickness of layers
      atm%mass,                                                                &
      ! Component volumes
      ukca_comp_vol,                                                           &
      ! Modal volumes, densities, and water content
      ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,                          &
      ! Modal mass-mixing ratios
      aer%mode_mix_ratio,                                                      &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Type selection
      ip_ukca_mode_accum, soluble_wanted,                                      &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Modal extinction aerosol opt depth (output)
      aod_ukca_acc_sol, sod_ukca_acc_sol,                                      &
      ! Modal absorption aerosol opt depth (output)
      aaod_ukca_acc_sol,                                                       &
      ! Fixed array dimensions
      dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode,                 &
      spectrum%dim%nd_aod_wavel)

    ! Put requested AOD diagnostics in diag structure

    IF (diag%l_aod_ukca_acc_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aod_ukca_acc_sol(col_list(l), row_list(l), i)                   &
            = aod_ukca_acc_sol(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_sod_ukca_acc_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%sod_ukca_acc_sol(col_list(l), row_list(l), i)                   &
            = sod_ukca_acc_sol(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_aaod_ukca_acc_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aaod_ukca_acc_sol(col_list(l), row_list(l), i)                  &
            = aaod_ukca_acc_sol(l, i)
        END DO
      END DO
    END IF
  END IF 

  ! Calculate AOD diagnostics for coarse soluble mode if requested
  IF (diag%l_aod_ukca_cor_sol .OR. diag%l_sod_ukca_cor_sol .OR.                &
       diag%l_aaod_ukca_cor_sol ) THEN

    CALL ukca_radaer_compute_aod(                                              &
      ! Actual array dimension
      atm%n_profile, atm%n_layer, n_ukca_mode,                                 &
      n_ukca_cpnt, spectrum%aerosol%n_aod_wavel,                               &
      ! UKCA_RADAER structure
      ukca_radaer,                                                             &
      ! Modal diameters from UKCA module
      ukca_dry_diam, ukca_wet_diam,                                            &
      ! Mass thickness of layers
      atm%mass,                                                                &
      ! Component volumes
      ukca_comp_vol,                                                           &
      ! Modal volumes, densities, and water content
      ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,                          &
      ! Modal mass-mixing ratios
      aer%mode_mix_ratio,                                                      &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Type selection
      ip_ukca_mode_coarse, soluble_wanted,                                     &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Modal extinction aerosol opt depth (output)
      aod_ukca_cor_sol, sod_ukca_cor_sol,                                      &
      ! Modal absorption aerosol opt depth (output)
      aaod_ukca_cor_sol,                                                       &
      ! Fixed array dimensions
      dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode,                 &
      spectrum%dim%nd_aod_wavel)

    ! Put requested AOD diagnostics in diag structure

    IF (diag%l_aod_ukca_cor_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aod_ukca_cor_sol(col_list(l), row_list(l), i)                   &
            = aod_ukca_cor_sol(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_sod_ukca_cor_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%sod_ukca_cor_sol(col_list(l), row_list(l), i)                   &
            = sod_ukca_cor_sol(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_aaod_ukca_cor_sol) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aaod_ukca_cor_sol(col_list(l), row_list(l), i)                  &
            = aaod_ukca_cor_sol(l, i)
        END DO
      END DO
    END IF
  END IF 

  ! Calculate AOD diagnostics for Aitken insoluble mode if requested
  IF (diag%l_aod_ukca_ait_ins .OR. diag%l_sod_ukca_ait_ins .OR.                &
       diag%l_aaod_ukca_ait_ins ) THEN

    CALL ukca_radaer_compute_aod(                                              &
      ! Actual array dimension
      atm%n_profile, atm%n_layer, n_ukca_mode,                                 &
      n_ukca_cpnt, spectrum%aerosol%n_aod_wavel,                               &
      ! UKCA_RADAER structure
      ukca_radaer,                                                             &
      ! Modal diameters from UKCA module
      ukca_dry_diam, ukca_wet_diam,                                            &
      ! Mass thickness of layers
      atm%mass,                                                                &
      ! Component volumes
      ukca_comp_vol,                                                           &
      ! Modal volumes, densities, and water content
      ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,                          &
      ! Modal mass-mixing ratios
      aer%mode_mix_ratio,                                                      &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Type selection
      ip_ukca_mode_aitken, soluble_unwanted,                                   &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Modal extinction aerosol opt depth (output)
      aod_ukca_ait_ins, sod_ukca_ait_ins,                                      &
      ! Modal absorption aerosol opt depth (output)
      aaod_ukca_ait_ins,                                                       &
      ! Fixed array dimensions
      dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode,                 &
      spectrum%dim%nd_aod_wavel)

    ! Put requested AOD diagnostics in diag structure

    IF (diag%l_aod_ukca_ait_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aod_ukca_ait_ins(col_list(l), row_list(l), i)                   &
            = aod_ukca_ait_ins(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_sod_ukca_ait_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%sod_ukca_ait_ins(col_list(l), row_list(l), i)                   &
            = sod_ukca_ait_ins(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_aaod_ukca_ait_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aaod_ukca_ait_ins(col_list(l), row_list(l), i)                  &
            = aaod_ukca_ait_ins(l, i)
        END DO
      END DO
    END IF
  END IF 

  ! Calculate AOD diagnostics for accumulation insoluble mode if requested
  IF (diag%l_aod_ukca_acc_ins .OR. diag%l_sod_ukca_acc_ins .OR.                &
       diag%l_aaod_ukca_acc_ins ) THEN

    CALL ukca_radaer_compute_aod(                                              &
      ! Actual array dimension
      atm%n_profile, atm%n_layer, n_ukca_mode,                                 &
      n_ukca_cpnt, spectrum%aerosol%n_aod_wavel,                               &
      ! UKCA_RADAER structure
      ukca_radaer,                                                             &
      ! Modal diameters from UKCA module
      ukca_dry_diam, ukca_wet_diam,                                            &
      ! Mass thickness of layers
      atm%mass,                                                                &
      ! Component volumes
      ukca_comp_vol,                                                           &
      ! Modal volumes, densities, and water content
      ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,                          &
      ! Modal mass-mixing ratios
      aer%mode_mix_ratio,                                                      &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Type selection
      ip_ukca_mode_accum, soluble_unwanted,                                    &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Modal extinction aerosol opt depth (output)
      aod_ukca_acc_ins, sod_ukca_acc_ins,                                      &
      ! Modal absorption aerosol opt depth (output)
      aaod_ukca_acc_ins,                                                       &
      ! Fixed array dimensions
      dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode,                 &
      spectrum%dim%nd_aod_wavel)

    ! Put requested AOD diagnostics in diag structure

    IF (diag%l_aod_ukca_acc_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aod_ukca_acc_ins(col_list(l), row_list(l), i)                   &
            = aod_ukca_acc_ins(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_sod_ukca_acc_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%sod_ukca_acc_ins(col_list(l), row_list(l), i)                   &
            = sod_ukca_acc_ins(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_aaod_ukca_acc_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aaod_ukca_acc_ins(col_list(l), row_list(l), i)                  &
            = aaod_ukca_acc_ins(l, i)
        END DO
      END DO
    END IF
  END IF 

  ! Calculate AOD diagnostics for coarse insoluble mode if requested
  IF (diag%l_aod_ukca_cor_ins .OR. diag%l_sod_ukca_cor_ins .OR.                &
       diag%l_aaod_ukca_cor_ins ) THEN

    CALL ukca_radaer_compute_aod(                                              &
      ! Actual array dimension
      atm%n_profile, atm%n_layer, n_ukca_mode,                                 &
      n_ukca_cpnt, spectrum%aerosol%n_aod_wavel,                               &
      ! UKCA_RADAER structure
      ukca_radaer,                                                             &
      ! Modal diameters from UKCA module
      ukca_dry_diam, ukca_wet_diam,                                            &
      ! Mass thickness of layers
      atm%mass,                                                                &
      ! Component volumes
      ukca_comp_vol,                                                           &
      ! Modal volumes, densities, and water content
      ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,                          &
      ! Modal mass-mixing ratios
      aer%mode_mix_ratio,                                                      &
      ! Modal number concentrations
      ukca_modal_number,                                                       &
      ! Type selection
      ip_ukca_mode_coarse, soluble_unwanted,                                   &
      ! Model level of the tropopause
      trindxrad,                                                               &
      ! Modal extinction aerosol opt depth (output)
      aod_ukca_cor_ins, sod_ukca_cor_ins,                                      &
      ! Modal absorption aerosol opt depth (output)
      aaod_ukca_cor_ins,                                                       &
      ! Fixed array dimensions
      dimen%nd_profile, dimen%nd_layer, dimen%nd_aerosol_mode,                 &
      spectrum%dim%nd_aod_wavel)

    ! Put requested AOD diagnostics in diag structure

    IF (diag%l_aod_ukca_cor_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aod_ukca_cor_ins(col_list(l), row_list(l), i)                   &
            = aod_ukca_cor_ins(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_sod_ukca_cor_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%sod_ukca_cor_ins(col_list(l), row_list(l), i)                   &
            = sod_ukca_cor_ins(l, i)
        END DO
      END DO
    END IF
    IF (diag%l_aaod_ukca_cor_ins) THEN
      DO i=1, spectrum%aerosol%n_aod_wavel
        DO l=1, atm%n_profile
          diag%aaod_ukca_cor_ins(col_list(l), row_list(l), i)                  &
            = aaod_ukca_cor_ins(l, i)
        END DO
      END DO
    END IF
  END IF 

  ! Counters for array indicies
  iwv_ext = 0
  iwv_abs = 0
  iwv_sca = 0
  iwv_gsca = 0

  DO i = 1, spectrum%aerosol%n_aod_wavel

    ! If requested then calculate 3D aerosol radiation diagnostics
    ! separately for each wavelength from ukca_radaer_3D_diags

    IF (diag%l_ukca_aerosol_ext(i)  .OR.                                       &
        diag%l_ukca_aerosol_abs(i)  .OR.                                       &
        diag%l_ukca_aerosol_sca(i)  .OR.                                       &
        diag%l_ukca_aerosol_gsca(i)) THEN 

      CALL ukca_radaer_3D_diags(                                               &
        ! Fixed array dimensions
        dimen%nd_profile,  dimen%nd_layer, dimen%nd_aerosol_mode,              &
        ! Actual array dimensions
        atm%n_profile, atm%n_layer, n_ukca_mode, n_ukca_cpnt,                  &
        ! UKCA_RADAER structure
        ukca_radaer,                                                           &
        ! Modal diameters from UKCA module
        ukca_dry_diam, ukca_wet_diam,                                          &
        ! Air density
        atm%density,                                                           &
        ! Component volumes
        ukca_comp_vol,                                                         &
        ! Modal volumes, densities, and water content
        ukca_modal_vol, ukca_modal_rho, ukca_modal_wtv,                        &
        ! Modal mass-mixing ratios
        aer%mode_mix_ratio,                                                    &
        ! Modal number concentrations
        ukca_modal_number,                                                     &
        ! Model level of the tropopause
        trindxrad,                                                             &
        ! Index of wavelength to consider
        i,                                                                     &
        ! Aerosol extinction and absorption profiles
        ukca_aerosol_ext, ukca_aerosol_abs,                                    &
        ! Aerosol scattering and asymmetry parameter
        ukca_aerosol_sca, ukca_aerosol_gsca                                    &
      )

      ! Put requested 3D diagnostics in diag structure

      IF (diag%l_ukca_aerosol_ext(i)) THEN
        iwv_ext = iwv_ext + 1
        DO k=1, atm%n_layer 
          DO l=1, atm%n_profile
            diag%ukca_aerosol_ext(col_list(l), row_list(l), k, iwv_ext)        &
               = ukca_aerosol_ext(l, atm%n_layer+1-k)
          END DO
        END DO
      END IF
      IF (diag%l_ukca_aerosol_sca(i)) THEN
        iwv_abs = iwv_abs + 1
        DO k=1, atm%n_layer 
          DO l=1, atm%n_profile
            diag%ukca_aerosol_sca(col_list(l), row_list(l), k, iwv_abs)        &
               = ukca_aerosol_sca(l, atm%n_layer+1-k)
          END DO
        END DO
      END IF
      IF (diag%l_ukca_aerosol_abs(i)) THEN
        iwv_sca = iwv_sca + 1
        DO k=1, atm%n_layer 
          DO l=1, atm%n_profile
            diag%ukca_aerosol_abs(col_list(l), row_list(l), k, iwv_sca)        &
               = ukca_aerosol_abs(l, atm%n_layer+1-k)
          END DO
        END DO
      END IF
      IF (diag%l_ukca_aerosol_gsca(i)) THEN
        iwv_gsca = iwv_gsca + 1
        DO k=1, atm%n_layer 
          DO l=1, atm%n_profile
            diag%ukca_aerosol_gsca(col_list(l), row_list(l), k, iwv_gsca)      &
               = ukca_aerosol_gsca(l, atm%n_layer+1-k)
          END DO
        END DO
      END IF
    END IF ! END calculations for UKCA 3D aerosol-radiation diagnostics
  END DO
END IF ! l_ukca_radaer .OR. l_glomap_clim_radaer

IF (l_easyaerosol_rad) THEN

  CALL set_easyaerosol_field(dimen, spectrum%basic%n_band, atm, &
                             aer, col_list, row_list, easyaerosol_rad)

  ! Update diagnostics, in m-1
  IF (diag%l_easyaerosol_extinction) THEN
    DO i = 1, spectrum%basic%n_band
      DO k = 1, model_levels 
        DO l = 1, atm%n_profile
          diag%easyaerosol_extinction(col_list(l), row_list(l), k, i) = &
              (aer%mode_absorption(l, atm%n_layer+1-k, aer%n_mode, i) + &
               aer%mode_scattering(l, atm%n_layer+1-k, aer%n_mode, i)) * &
               atm%density(l, atm%n_layer+1-k)
        END DO ! l
      END DO ! k
    END DO ! i
  END IF
  IF (diag%l_easyaerosol_absorption) THEN
    DO i = 1, spectrum%basic%n_band
      DO k = 1, model_levels 
        DO l = 1, atm%n_profile
          diag%easyaerosol_absorption(col_list(l), row_list(l), k, i) = &
               aer%mode_absorption(l, atm%n_layer+1-k, aer%n_mode, i) * &
               atm%density(l, atm%n_layer+1-k)
        END DO ! l
      END DO ! k
    END DO ! i
  END IF
  IF (diag%l_easyaerosol_scattering) THEN
    DO i = 1, spectrum%basic%n_band
      DO k = 1, model_levels 
        DO l = 1, atm%n_profile
          diag%easyaerosol_scattering(col_list(l), row_list(l), k, i) = &
               aer%mode_scattering(l, atm%n_layer+1-k, aer%n_mode, i) * &
               atm%density(l, atm%n_layer+1-k)
        END DO ! l
      END DO ! k
    END DO ! i
  END IF
  IF (diag%l_easyaerosol_asytimscat) THEN
    DO i = 1, spectrum%basic%n_band
      DO k = 1, model_levels 
        DO l = 1, atm%n_profile
          diag%easyaerosol_asytimscat(col_list(l), row_list(l), k, i) = &
               aer%mode_asymmetry(l, atm%n_layer+1-k, aer%n_mode, i)  * &
               aer%mode_scattering(l, atm%n_layer+1-k, aer%n_mode, i) * &
               atm%density(l, atm%n_layer+1-k)
        END DO ! l
      END DO ! k
    END DO ! i
  END IF

END IF ! l_easyaerosol_rad

9999 CONTINUE
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_aer
END MODULE set_aer_mod
