! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Shortwave interface to the socrates core radiation code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE sw_rad_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SW_RAD_MOD'
CONTAINS

SUBROUTINE sw_rad(                                                             &
!                   Mixing ratios
      h2o, co2, o3, o2_mix_ratio                                               &
    , co2_dim1, co2_dim2, co2_3d, l_co2_3d                                     &
    , n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio                              &
!                   Chemical greenhouse gas fields
    , ngrgas, grgas_field                                                      &
!                   Thermodynamic variables
    , pstar                                                                    &
    , p_layer_boundaries                                                       &
    , p_layer_centres                                                          &
    , t_layer_centres                                                          &
    , t_layer_boundaries, p_extra_layer, t_extra_layer                         &
    , d_mass, density, layer_heat_capacity                                     &
    , r_layer_centres, r_layer_boundaries                                      &
!                   Options for COSP
    , l_cosp                                                                   &
!                   Options for treating clouds
    , l_inhom_cloud, inhom_cloud, dp_corr_strat, dp_corr_conv                  &
!                   Stratiform cloud fields
    , lca_area, lca_bulk, lccwc1, lccwc2, n_drop_pot                           &
!                   Convective cloud fields
    , cca, cccwp, ccw, lcbase, ccb, cct                                        &
!                   Surface fields
    , land_albedo, flandg, sea_ice_albedo                                      &
    , open_sea_albedo, ice_fraction, land, land0p5, lying_snow                 &
!                   Solar fields
    , coszin, lit, solar_constant, list, scs, sindec                           &
    , cos_zen_sph, day_frac_sph                                                &
!                   Aerosol fields
    , l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero                   &
    , bl_depth, n_levels_bl                                                    &
    , l_dust, l_use_dust, dust_dim1, dust_dim2                                 &
    , dust_1, dust_2, dust_3, dust_4, dust_5, dust_6                           &
    , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic                   &
    , l_sulpc_so2, l_use_sulpc_direct, l_use_sulpc_indirect                    &
    , sulp_dim1, sulp_dim2                                                     &
    , accum_sulphate, aitken_sulphate, diss_sulphate                           &
    , sea_salt_film, sea_salt_jet, l_use_seasalt_indirect                      &
    , l_use_seasalt_direct, salt_dim1, salt_dim2, salt_dim3                    &
    , l_soot, l_use_soot_direct, soot_dim1, soot_dim2                          &
    , fresh_soot, aged_soot                                                    &
    , l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2                    &
    , fresh_bmass, aged_bmass, cloud_bmass, l_use_bmass_indirect               &
    , l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2                          &
    , fresh_ocff, aged_ocff, cloud_ocff, l_use_ocff_indirect                   &
    , l_nitrate, l_use_nitrate_direct, nitrate_dim1, nitrate_dim2              &
    , accum_nitrate, diss_nitrate, l_use_nitrate_indirect                      &
    , l_use_arcl, arcl_dim1, arcl_dim2, n_arcl_species                         &
    , n_arcl_compnts, i_arcl_compnts, arcl                                     &
    , aero_meso, l_murk_rad                                                    &
    , l_ukca_radaer, l_glomap_clim_radaer, ukca_radaer, ukca_dim1, ukca_dim2   &
    , ukca_mmr, ukca_cvl, ukca_dry, ukca_wet                                   &
    , ukca_rho, ukca_vol, ukca_wtv, ukca_nbr                                   &
    , l_easyaerosol_rad, easyaerosol_rad                                       &
    , l_easyaerosol_cdnc, easyaerosol_cdnc                                     &
!                   Time
    , previous_time, seconds_since_midnight                                    &
!                   Grid-dependent arrays
    , true_latitude, true_longitude, obs_solid_angle, trans_solid_angle        &
    , dir_flux_to_trans                                                        &
!                   Level of tropopause
    , trindx                                                                   &
!                   Spectrum
    , spectrum                                                                 &
!                   Algorithmic options
    , control, pts                                                             &
!                   diagnostics
    , sw_diag, row_list, col_list                                              &
!                   Physical dimensions
    , dimen, nlit, n_points, n_layer, nclds                                    &
    , nozone, row_length, rows                                                 &
    , nd_field, nd_max_bands, nd_field_flux_diag, nd_field_rad_diag            &
    , n_cca_lev, n_ukca_mode, n_ukca_cpnt                                      &
!                   Output
    , surf_down_sw, flux_below_690nm_surf                                      &
    , netsw, top_absorption, swsea, swout                                      &
!                   COSP input arguments
    , cosp_gbx, cosp_sgx                                                       &
    )

USE rad_pcf
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_atm,      ONLY: StrAtm
USE def_cld,      ONLY: StrCld
USE def_aer,      ONLY: StrAer
USE def_bound,    ONLY: StrBound
USE def_out,      ONLY: StrOut
USE def_diag,     ONLY: StrDiag
USE ukca_radaer_struct_mod, ONLY: ukca_radaer_struct
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE cosp_types_mod, ONLY: cosp_gridbox, cosp_subgrid
USE arcl_mod,    ONLY: npd_arcl_species, npd_arcl_compnts
USE nlsizes_namelist_mod, ONLY: model_levels
USE def_easyaerosol, ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

USE socrates_calc_mod, ONLY: socrates_calc
USE socrates_init_mod, ONLY: socrates_init
USE socrates_postproc_mod, ONLY: socrates_postproc
IMPLICIT NONE


! Dummy arguments

! Dimensions of arrays:
INTEGER, INTENT(IN) :: row_length
!                          length of rows on each domain
INTEGER, INTENT(IN) :: rows
!                          number of rows in the domain
INTEGER ::                                                                     &
          !, intent(in)
    nd_field                                                                   &
!       Field size in calling program
    , nd_max_bands                                                             &
!       Maximum number of SW bands in all SW calls
    , nd_field_flux_diag                                                       &
!       Field size for flux diagnostics
    , nd_field_rad_diag
!       Field size for radiance diagnostics

! Actual sizes used:
INTEGER ::                                                                     &
          !, intent(in)
    n_points                                                                   &
!       Number of points to be diagnosed including unlit points
    , nozone                                                                   &
!       Number of levels with ozone
    , n_layer                                                                  &
!       number of layers seen in the radiation scheme
    , nclds                                                                    &
!       Number of cloudy levels
    , n_levels_bl                                                              &
!       Number of layers occupied by boundary-layer aerosol
!       if l_clim_aero_hgt is false.
    , n_cca_lev                                                                &
!       Number of convective cloud levels
    , n_ukca_mode                                                              &
!       Number of aerosol modes in UKCA_RADAER
    , n_ukca_cpnt
!       Number of aerosol components in UKCA_RADAER

! Spectral data:
TYPE (StrSpecData), INTENT(IN) :: spectrum

! Controlling options:
TYPE (StrCtrl), INTENT(IN) :: control

! Dimensions:
TYPE (StrDim), INTENT(IN) :: dimen

! Gaseous mixing ratios
REAL ::                                                                        &
          !, intent(in)
    h2o(nd_field, model_levels)                                                &
!       Mass mixing ratio of water
    , co2                                                                      &
!       Mass mixing ratio of co2
    , o3(nd_field, nozone)                                                     &
!       Mass mixing ratios of ozone
    , o2_mix_ratio                                                             &
!       Mass mixing ratio of oxygen
    , n2o_mix_ratio                                                            &
!       Mass mixing ratio of nitrous oxide
    , ch4_mix_ratio                                                            &
!       Mass mixing ratio of methane
    , so2_mix_ratio
!       Mass mixing ratio of sulphur dioxide

! Chemical greenhouse gas fields
INTEGER, INTENT(IN) :: ngrgas
REAL, INTENT(IN) :: grgas_field(nd_field, model_levels, ngrgas)

! General atmospheric properties:
REAL ::                                                                        &
          !, intent(in)
    pstar(nd_field)                                                            &
!       Surface pressures
    , p_layer_boundaries(nd_field,0:model_levels)                              &
!        pressure at boundaries of layers
    , p_layer_centres(nd_field,0:model_levels)                                 &
!        pressure at centres of layers
    , t_layer_centres(nd_field, model_levels)
!       Temperatures at centres of layers

REAL, INTENT(IN) :: t_layer_boundaries(nd_field, 0:model_levels)
!   Temperature at layer boundaries
REAL, INTENT(IN) :: p_extra_layer(nd_field)
!   Pressure at centre of extra top layer
REAL, INTENT(IN) :: t_extra_layer(nd_field)
!   Temperature at centre of extra top layer
REAL, INTENT(IN) :: d_mass(nd_field, model_levels+1)
!   Mass of layer (kg m-2)
REAL, INTENT(IN) :: density(nd_field, model_levels+1)
!   Density of layer (kg m-3)
REAL, INTENT(IN) :: layer_heat_capacity(nd_field, model_levels)
!   Heat capacity of layer
REAL, INTENT(IN) :: r_layer_centres(nd_field, model_levels+1)
!   Height above centre of planet / radius at centres of layers
REAL, INTENT(IN) :: r_layer_boundaries(nd_field, 0:model_levels+1)
!   Height above centre of planet / radius at layer boundaries

! Incident solar radiation:
INTEGER ::                                                                     &
          !, intent(in)
    nlit                                                                       &
!       Number of lit points
    , list(nd_field)
!       List of lit points
REAL, INTENT(IN) ::                                                            &
    coszin(nd_field)                                                           &
!       Cosines of zenith angle
    , solar_constant                                                           &
!       Total solar irradiance at 1 AU
    , scs                                                                      &
!       Scaling of solar incident field
    , lit(nd_field)                                                            &
!       Fraction of time point is lit
    , sindec                                                                   &
!       sin(solar declination)
    , cos_zen_sph(nd_field, 0:n_layer+1)                                       &
!       Cosines of zenith angle for each layer
    , day_frac_sph(nd_field, 0:n_layer+1)
!       Fraction of time point is lit within each layer

! Flag for COSP
LOGICAL,INTENT(IN) :: l_cosp

! Options for treating clouds
LOGICAL ::                                                                     &
          !, intent(in)
    l_inhom_cloud
!       Flag to use scaling factors for inhomogeneous cloud
REAL ::                                                                        &
          !, intent(in)
    inhom_cloud(dimen%nd_cloud_component)                                      &
!       Scaling factors for inhomogeneous cloud
    , dp_corr_strat                                                            &
!       Decorrelation pressure scale for large scale cloud
    , dp_corr_conv
!       Decorrelation pressure scale for convective cloud

! Properties of stratiform clouds:
REAL ::                                                                        &
          !, intent(in)
    lccwc1(nd_field, nclds+1/(nclds+1))                                        &
!       Nominal liquid water contents
    , lccwc2(nd_field, nclds+1/(nclds+1))                                      &
!       Nominal ice water contents
    , lca_area(nd_field, nclds+1/(nclds+1))                                    &
!       Area fractions of layer clouds outside convective towers
    , lca_bulk(nd_field, nclds+1/(nclds+1))
!       Bulk fractions of layer clouds outside convective towers

REAL, INTENT(IN) :: n_drop_pot(nd_field, nclds)

! Properties of convective clouds:
INTEGER ::                                                                     &
          !, intent(in)
    ccb(nd_field)                                                              &
!       Base of convective cloud
    , lcbase(nd_field)                                                         &
!       Base of convective cloud (corrected)
    , cct(nd_field)
!       Top of convective cloud
REAL ::                                                                        &
          !, intent(in)
    cccwp(nd_field)                                                            &
!       Water path of convective cloud
    , ccw(nd_field, model_levels)                                              &
!       Convective cloud water
    , cca(nd_field,n_cca_lev)
!       Fraction of convective cloud

! Aerosols:
LOGICAL ::                                                                     &
          !, intent(in)
    l_climat_aerosol                                                           &
!       Flag for climatological aerosol
    , l_clim_aero_hgt                                                          &
!       flag to use the depth of the boundary layer to set
!       the climatological aerosol
    , l_hadgem1_clim_aero                                                      &
!       Flag to use HadGEM1 setting for climatological aerosols
    , l_murk_rad
!       flag for mesoscale model aerosol
LOGICAL ::                                                                     &
          !, intent(in)
    l_sulpc_so2                                                                &
!       Sulphur cycle available for effecting fluxes or diagnostics
    , l_use_sulpc_direct                                                       &
!       Flag to use sulphur cycle for direct effect
    , l_use_sulpc_indirect                                                     &
!       Flag to use sulphur cycle for indirect effect
    , l_dust                                                                   &
!       Dust is available for effecting fluxes or diagnostics
    , l_use_dust                                                               &
!       Flag to use direct rad effect of mineral dust
    , l_use_biogenic                                                           &
!       Flag to use biogenic for direct effect
    , l_soot                                                                   &
!        Soot is available for effecting fluxes or diagnostics
    , l_use_soot_direct                                                        &
!       Use direct rad. effect of soot aerosol
    , l_biomass                                                                &
!        Biomass is available for effecting fluxes or diagnostics
    , l_use_bmass_direct                                                       &
!       Flag to use direct rad. effect of biomass smoke
    , l_use_bmass_indirect                                                     &
!       Flag to use indirect effect of biomass smoke
    , l_ocff                                                                   &
!        OCFF is available for effecting fluxes or diagnostics
    , l_use_ocff_direct                                                        &
!       Flag to use direct rad. effect of ocff
    , l_use_ocff_indirect                                                      &
!       Flag to use indirect effect of ocff
    , l_use_seasalt_indirect                                                   &
!       Flag to use sea-salt for indirect effect
    , l_use_seasalt_direct                                                     &
!       Flag to use sea-salt for direct effect
    , l_nitrate                                                                &
!        Nitrate is available for effecting fluxes or diagnostics
    , l_use_nitrate_direct                                                     &
!       Flag to use nitrate for direct effect
    , l_use_nitrate_indirect                                                   &
!       Flag to use nitrate for indirect effect
    , l_easyaerosol_rad                                                        &
!       Flag to use EasyAerosol for direct effect
    , l_easyaerosol_cdnc
!       Flag to use EasyAerosol cloud droplet number concentrations

INTEGER ::                                                                     &
          !, intent(in)
    sulp_dim1,sulp_dim2                                                        &
!       Dimensions for _sulphate arrays, (P_FIELD,P_LEVELS or 1,1)
    , dust_dim1, dust_dim2                                                     &
!       Dimensions for mineral dust arrays (p_field,p_levels or 1,1)
    , biogenic_dim1, biogenic_dim2                                             &
!       dimensions for biogenic array passed down to
!       r2_set_aerosol_field if direct effect required.
    , soot_dim1, soot_dim2                                                     &
!       Dimensions for soot arrays (P_FIELD,P_LEVELS or 1,1)
    , bmass_dim1, bmass_dim2                                                   &
!       Dimensions for biomass arrays (P_FIELD,P_LEVELS or 1,1)
    , ocff_dim1, ocff_dim2                                                     &
!       Dimensions for ocff arrays (P_FIELD,P_LEVELS or 1,1)
    , nitrate_dim1, nitrate_dim2                                               &
!       Dimensions for nitrate arrays (P_FIELD,P_LEVELS or 1,1)
    , salt_dim1, salt_dim2, salt_dim3
!       dimensions for salt arrays on input (salt_dim1*salt_dim2=p_field
!       and salt_dim3=p_levels, or else 1,1,1)
REAL ::                                                                        &
          !, intent(in)
    accum_sulphate(sulp_dim1, sulp_dim2)                                       &
!       Mass mixing ratio of accumulation mode aerosol
    , aitken_sulphate(sulp_dim1, sulp_dim2)                                    &
!       Mass mixing ratio of aitken mode aerosol
    , diss_sulphate(sulp_dim1, sulp_dim2)                                      &
!       Mixing ratio of dissolved sulphate
    , dust_1(dust_dim1, dust_dim2)                                             &
!       Mass mixing ratio of div1 dust
    , dust_2(dust_dim1, dust_dim2)                                             &
!       Mass mixing ratio of div2 dust
    , dust_3(dust_dim1, dust_dim2)                                             &
!       Mass mixing ratio of div3 dust
    , dust_4(dust_dim1, dust_dim2)                                             &
!       Mass mixing ratio of div4 dust
    , dust_5(dust_dim1, dust_dim2)                                             &
!       Mass mixing ratio of div5 dust
    , dust_6(dust_dim1, dust_dim2)                                             &
!       Mass mixing ratio of div6 dust
    , biogenic(biogenic_dim1, biogenic_dim2)                                   &
!       Mixing ratios of biogenic aerosol
    , sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                           &
!         number concentration of film-mode sea-salt aerosol
    , sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                            &
!         number concentration of jet-mode sea-salt aerosol
    , fresh_soot(soot_dim1, soot_dim2)                                         &
!       Soot mixing ratios
    , aged_soot(soot_dim1, soot_dim2)                                          &
!       Soot mixing ratios
    , fresh_bmass(bmass_dim1, bmass_dim2)                                      &
!       Mass mixing ratio of fresh biomass smoke
    , aged_bmass(bmass_dim1, bmass_dim2)                                       &
!       Mass mixing ratio of aged biomass smoke
    , cloud_bmass(bmass_dim1, bmass_dim2)                                      &
!       Mass mixing ratio of in-cloud biomass smoke
    , fresh_ocff(ocff_dim1, ocff_dim2)                                         &
!       Mass mixing ratio of fresh fossil-fuel organic carbon aer
    , aged_ocff(ocff_dim1, ocff_dim2)                                          &
!       Mass mixing ratio of aged fossil-fuel organic carbon aer
    , cloud_ocff(ocff_dim1, ocff_dim2)                                         &
!       Mass mixing ratio of in-cloud fossil-fuel org carbon aer
    , accum_nitrate(nitrate_dim1, nitrate_dim2)                                &
!       Mass mixing ratio of accumulation nitrate aerosol
    , diss_nitrate(nitrate_dim1, nitrate_dim2)                                 &
!       Mass mixing ratio of dissolved nitrate aerosol
    , bl_depth(nd_field)                                                       &
!       depth of the boundary layer
    , aero_meso(nd_field, model_levels)
!       mixing ratio of 'urban' aerosol of mesoscale model

! Aerosol climatology for NWP

! Number of requested species within the climatology
INTEGER :: n_arcl_species

! Corresponding number of requested components
INTEGER :: n_arcl_compnts

! Model switch for each species
LOGICAL :: l_use_arcl(npd_arcl_species)

! Index of each component
INTEGER :: i_arcl_compnts(npd_arcl_compnts)

! Array dimensions
INTEGER :: arcl_dim1, arcl_dim2

! Mass-mixing ratios
REAL :: arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)

! UKCA_RADAER:

! Model switch
LOGICAL :: l_ukca_radaer
LOGICAL :: l_glomap_clim_radaer

! UKCA_RADAER structure
TYPE (ukca_radaer_struct) :: ukca_radaer

! Dimensions
INTEGER :: ukca_dim1
INTEGER :: ukca_dim2

! Component mass mixing ratios and volumes
REAL :: ukca_mmr(ukca_dim1, ukca_dim2, n_ukca_cpnt)
REAL :: ukca_cvl(ukca_dim1, ukca_dim2, n_ukca_cpnt)

! Dry and wet modal diameters
REAL :: ukca_dry(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL :: ukca_wet(ukca_dim1, ukca_dim2, n_ukca_mode)

! Modal densities, volumes, volume of water, and number conc
REAL :: ukca_rho(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL :: ukca_vol(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL :: ukca_wtv(ukca_dim1, ukca_dim2, n_ukca_mode)
REAL :: ukca_nbr(ukca_dim1, ukca_dim2, n_ukca_mode)

! EasyAerosol climatology
TYPE (t_easyaerosol_rad)  :: easyaerosol_rad
TYPE (t_easyaerosol_cdnc) :: easyaerosol_cdnc

! Carbon cycle:
LOGICAL ::                                                                     &
    l_co2_3d
!       Controls use of 3D CO2 field
INTEGER ::                                                                     &
          !, intent(in)
    co2_dim1, co2_dim2
!       Dimensions for CO2 array, (P_FIELD,P_LEVELS or 1,1)
REAL ::                                                                        &
          !, intent(in)
    co2_3d(co2_dim1, co2_dim2)
!       Mass mixing ratio of carbon dioxide

! Properties of the surface:
LOGICAL, INTENT(IN) :: land(nd_field)
!                         Land mask
LOGICAL, INTENT(IN) :: land0p5(nd_field)
!                         Land mask (TRUE if land fraction > 0.5)
REAL ::                                                                        &
          !, intent(in)
    ice_fraction(nd_field)                                                     &
!         fraction of sea ice in sea portion of grid box
    , land_albedo(nd_field,4)                                                  &
!         land surface albedo fields
!         (*,1) - direct beam visible
!         (*,2) - diffuse visible
!         (*,3) - direct beam near-ir
!         (*,4) - diffuse near-ir
    , flandg(nd_field)                                                         &
!         land fraction in grid box
    , sea_ice_albedo(nd_field,4)                                               &
!         sea ice albedo fields
!         (*,1) - direct beam visible
!         (*,2) - diffuse visible
!         (*,3) - direct beam near-ir
!         (*,4) - diffuse near-ir
    , open_sea_albedo(nd_field, 2, nd_max_bands)                               &
!       Surface albedo field of open sea
!         (direct and diffuse components)
    , lying_snow(nd_field)
!       Mass loading of lying snow

REAL, INTENT(IN) :: true_latitude(nd_field)
REAL, INTENT(IN) :: true_longitude(nd_field)
REAL, INTENT(IN) :: seconds_since_midnight
INTEGER, INTENT(IN) :: previous_time(7)

! Solid angle subtended by grid-box for an observer at 1 AU
REAL, INTENT(IN) :: obs_solid_angle(nd_field)
! Solid angle subtended by grid-box for a transit observer at 1 AU
REAL, INTENT(IN) :: trans_solid_angle(nd_field)
! Conversion factor from direct flux to flux at transit observer
REAL, INTENT(IN) :: dir_flux_to_trans(nd_field)

!                   level of tropopause
INTEGER :: trindx(nd_field)
!       The layer boundary of the tropopause

! Increment of time:
REAL ::                                                                        &
          !, intent(in)
    pts
!       Time increment

! Calculated fluxes:
REAL ::                                                                        &
          !, intent(out)
    swout(nd_field, model_levels+2)                                            &
!       Net downward fluxes
    , swsea(nd_field)                                                          &
!       Sea-surface components of flux
!       weighted by (open sea)/(total sea) fraction
    , netsw(nd_field)                                                          &
!       Net absorbed shortwave radiation
    , flux_below_690nm_surf(nd_field_flux_diag)                                &
!       Net surface flux below 690nm (at points where there
!       is sea-ice this is weighted by the fraction of open sea.)
    , surf_down_sw(nd_field_flux_diag, 4)                                      &
!         surface downward shortwave radiation components
!         (*,1) - direct beam visible
!         (*,2) - diffuse visible
!         (*,3) - direct beam near-ir
!         (*,4) - diffuse near-ir
    , top_absorption(nd_field)
!         radiative absorption above the top of the atmosphere
!         as seen in the main model


! Diagnostics:
TYPE (StrDiag) :: sw_diag

INTEGER, INTENT(IN) :: row_list(nd_field)
!                          list of row indices of lit points
INTEGER, INTENT(IN) :: col_list(nd_field)
!                          list of column indices of lit points

! Structure with COSP inputs
TYPE(cosp_gridbox),INTENT(INOUT) :: cosp_gbx
TYPE(cosp_subgrid),INTENT(INOUT) :: cosp_sgx


! Local variables.

! Atmospheric properties:
TYPE(StrAtm) :: atm

! Cloud properties:
TYPE(StrCld) :: cld

! Aerosol properties:
TYPE(StrAer) :: aer

! Boundary conditions:
TYPE(StrBound) :: bound

! Output fields from core radiation code:
TYPE(StrOut) :: radout

! Output fields from diagnostic calls:
TYPE(StrOut) :: radout_clean
TYPE(StrOut) :: radout_forc


CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SW_RAD'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise the data for shortwave radiation computation. 
CALL socrates_init(                                                            &
  ! Mixing ratios
  h2o, co2, o3, o2_mix_ratio, co2_dim1, co2_dim2, co2_3d, l_co2_3d,            &
  n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio,                                 &
  ! Chemical greenhouse gas fields
  ngrgas, grgas_field,                                                         &
  ! Thermodynamic variables
  pstar, p_layer_boundaries, p_layer_centres, t_layer_centres,                 &
  t_layer_boundaries, p_extra_layer, t_extra_layer, d_mass, density,           &
  r_layer_centres, r_layer_boundaries,                                         &
  ! Options for treating clouds
  l_inhom_cloud, inhom_cloud, dp_corr_strat, dp_corr_conv,                     &
  ! Stratiform cloud fields
  lca_area, lca_bulk, lccwc1, lccwc2, n_drop_pot,                              &
  ! Convective cloud fields
  cca, cccwp, ccw, lcbase, ccb, cct,                                           &
  ! Surface fields
  land_albedo, flandg, sea_ice_albedo,                                         &
  open_sea_albedo, ice_fraction, land, land0p5, lying_snow,                    &
  ! Solar fields
  coszin, lit, solar_constant, list, scs, sindec,                              &
  cos_zen_sph, day_frac_sph,                                                   &
  ! Aerosol fields
  l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,                      &
  bl_depth, n_levels_bl,                                                       &
  l_dust, l_use_dust, dust_dim1, dust_dim2,                                    &
  dust_1, dust_2, dust_3, dust_4, dust_5, dust_6,                              &
  l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic,                      &
  l_sulpc_so2, l_use_sulpc_direct, sulp_dim1, sulp_dim2,                       &
  accum_sulphate, aitken_sulphate, diss_sulphate,                              &
  sea_salt_film, sea_salt_jet, l_use_seasalt_indirect,                         &
  l_use_seasalt_direct, salt_dim1, salt_dim2, salt_dim3,                       &
  l_soot, l_use_soot_direct, soot_dim1, soot_dim2, fresh_soot, aged_soot,      &
  l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2,                       &
  fresh_bmass, aged_bmass, cloud_bmass, l_use_bmass_indirect,                  &
  l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2,                             &
  fresh_ocff, aged_ocff, cloud_ocff, l_use_ocff_indirect,                      &
  l_nitrate, l_use_nitrate_direct, nitrate_dim1, nitrate_dim2,                 &
  accum_nitrate, diss_nitrate, l_use_nitrate_indirect,                         &
  l_use_arcl, arcl_dim1, arcl_dim2, n_arcl_species,                            &
  n_arcl_compnts, i_arcl_compnts, arcl, aero_meso, l_murk_rad,                 &
  l_ukca_radaer, l_glomap_clim_radaer, ukca_radaer, ukca_dim1, ukca_dim2,      &
  ukca_mmr, ukca_cvl, ukca_dry, ukca_wet,                                      &
  ukca_rho, ukca_vol, ukca_wtv, ukca_nbr,                                      &
  l_easyaerosol_rad, easyaerosol_rad,                                          &
  l_easyaerosol_cdnc, easyaerosol_cdnc,                                        &
  ! Time
  previous_time, seconds_since_midnight,                                       &
  ! Grid-dependent arrays
  true_latitude, true_longitude, trindx,                                       &
  ! Spectrum, controlling options and diagnostics
  spectrum, control, sw_diag, row_list, col_list,                              &
  ! Physical dimensions
  dimen, nlit, n_layer, nclds, nozone, row_length, rows,                       &
  nd_field, nd_max_bands, n_cca_lev, n_ukca_mode, n_ukca_cpnt,                 &
  ! Output radiance core data
  atm, cld, aer, bound)


! Compute the short wave calculations
CALL socrates_calc(                                                            &
  ! Controlling options and dimensions
  spectrum, control, sw_diag, dimen,                                           &
  ! Physical inputs
  atm, cld, aer, bound,                                                        &
  ! Outputs
  radout, radout_clean, radout_forc)


! Post-process the short wave data.
CALL socrates_postproc(                                                        &
  ! Thermodynamic variables
  layer_heat_capacity,                                                         &
  ! Options for COSP
  l_cosp,                                                                      &
  ! Surface fields
  flandg, ice_fraction, land,                                                  &
  ! Solar fields and grid-dependent arrays
  coszin, lit, day_frac_sph, list,                                             &
  obs_solid_angle, trans_solid_angle, dir_flux_to_trans,                       &
  trindx,                                                                      &
  ! Spectrum, controlling options and diagnostics
  spectrum, control, pts, sw_diag, row_list, col_list,                         &
  ! Physical dimensions
  dimen, nlit, n_points, n_layer, nclds, row_length, rows,                     &
  nd_field, nd_field_flux_diag,                                                &
  ! Radiance core data
  atm, cld, aer, bound, radout, radout_clean, radout_forc,                     &
  ! Output
  surf_down_sw, flux_below_690nm_surf,                                         &
  netsw, top_absorption, swsea, swout,                                         &
  cosp_gbx, cosp_sgx)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE sw_rad
END MODULE sw_rad_mod
