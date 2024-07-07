! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Initialise the data for the core radiation calculations.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------

MODULE socrates_init_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOCRATES_INIT_MOD'
CONTAINS

SUBROUTINE socrates_init(                                                      &
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
    , d_mass, density                                                          &
    , r_layer_centres, r_layer_boundaries                                      &
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
    , l_sulpc_so2, l_use_sulpc_direct                                          &
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
    , true_latitude, true_longitude                                            &
!                   Level of tropopause
    , trindx                                                                   &
!                   Spectrum
    , spectrum                                                                 &
!                   Algorithmic options
    , control                                                                  &
!                   diagnostics
    , sw_diag, row_list, col_list                                              &
!                   Physical dimensions
    , dimen, nlit, n_layer, nclds                                              &
    , nozone, row_length, rows                                                 &
    , nd_field, nd_max_bands                                                   &
    , n_cca_lev, n_ukca_mode, n_ukca_cpnt                                      &
!                   Output radiance core data
    , atm, cld, aer, bound                                                     &
    )

USE rad_pcf
USE conversions_mod, ONLY: recip_pi_over_180
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_atm,      ONLY: StrAtm
USE def_cld,      ONLY: StrCld,   allocate_cld,                                &
                                  allocate_cld_prsc,                           &
                                  allocate_cld_mcica 
USE def_aer,      ONLY: StrAer
                                                      
USE def_bound,    ONLY: StrBound
USE def_diag,     ONLY: StrDiag
USE mcica_mod
USE ukca_radaer_struct_mod, ONLY: ukca_radaer_struct
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE cosp_types_mod, ONLY: cosp_gridbox, cosp_subgrid
USE ereport_mod, ONLY: ereport
USE arcl_mod,    ONLY: npd_arcl_species, npd_arcl_compnts
USE nlsizes_namelist_mod, ONLY: model_levels
USE errormessagelength_mod, ONLY: errormessagelength
USE set_atm_mod, ONLY: set_atm
USE def_easyaerosol, ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

USE mcica_order_mod, ONLY: mcica_order
USE r2_cloud_level_diag_mod, ONLY: r2_cloud_level_diag
USE r2_set_cloud_field_mod, ONLY: r2_set_cloud_field
USE r2_set_cloud_parametrization_mod, ONLY: r2_set_cloud_parametrization
USE set_aer_mod, ONLY: set_aer
USE set_bound_mod, ONLY: set_bound
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
    , nd_max_bands                                                             
!       Maximum number of SW bands in all SW calls

! Actual sizes used:
INTEGER ::                                                                     &
          !, intent(in)
      nozone                                                                   &
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

! Atmospheric properties:
TYPE(StrAtm), INTENT(INOUT) :: atm

! Cloud properties:
TYPE(StrCld), INTENT(INOUT) :: cld

! Aerosol properties:
TYPE(StrAer), INTENT(INOUT) :: aer

! Boundary conditions:
TYPE(StrBound), INTENT(INOUT) :: bound

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
REAL ::                                                                        &
          !, intent(in)
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

!                   level of tropopause
INTEGER :: trindx(nd_field)
!       The layer boundary of the tropopause

! Diagnostics:
TYPE (StrDiag) :: sw_diag

INTEGER, INTENT(IN) :: row_list(nd_field)
!                          list of row indices of lit points
INTEGER, INTENT(IN) :: col_list(nd_field)
!                          list of column indices of lit points

! Local variables.

INTEGER :: i, j, k, l, ll, lll
!       Loop variables

! General atmospheric properties:
REAL, PARAMETER :: nullmmr = 0.0
!       Null mass mixing ratio
REAL :: alat(dimen%nd_profile)
!       Latitude in degrees

! Cloudy properties:
INTEGER ::                                                                     &
    i_cloud_tmp
!       Cloud Overlap used by r2_cloud_level_diag
REAL ::                                                                        &
    condensed_min_dim(dimen%nd_cloud_component)                                &
!       Minimum dimensions of condensed components
    , condensed_max_dim(dimen%nd_cloud_component)                              &
!       Maximum dimensions of condensed components
    , CDNC(dimen%nd_profile, dimen%id_cloud_top : dimen%nd_layer,              &
           dimen%nd_cloud_component)
!       CDNC array for call to r2_set_cloud_field from SW

! Surface properties:
LOGICAL :: land0p5_g(dimen%nd_profile)
!       Gathered land mask (true if land fraction >0.5)
REAL :: flandg_g(dimen%nd_profile)
!       Gathered land fraction


! Variables required for compatibility with subroutines:
REAL :: dummy1d(1)

! Logical flags for diagnostics
LOGICAL :: l_all_temps
!   If TRUE, the routine has been called to obtain diagnostics
!   for clouds consisting of liquid water at any temperature
!   (as done in MODIS retrievals). If FALSE, only clouds with
!   temperatures above freezing are to be diagnosed (as done
!   in AVHRR retrievals).
LOGICAL :: l_radius
!   If TRUE, the routine has been called for calculating cloud
!   droplet effective radius at cloud top. If FALSE then it's
!   been called for calculating cloud droplet number concentration
!   at cloud top.

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SOCRATES_INIT'
CHARACTER (LEN=errormessagelength) :: cmessage
INTEGER :: ierr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate structures for the core radiation code interface
CALL allocate_cld(cld, dimen, spectrum)
CALL allocate_cld_prsc(cld, dimen, spectrum)

! Initialize the error flag for the radiation code.
ierr=i_normal


! Set the properties at the boundaries (surface and top-of-atmosphere)
CALL set_bound(control, dimen, spectrum, bound, sw_diag,                       &
  nd_field, nlit, n_layer, list, row_list, col_list,                           &
  land, flandg, ice_fraction,                                                  &
  dummy1d, dummy1d, dummy1d, dummy1d, dummy1d,                                 &
  open_sea_albedo, sea_ice_albedo, land_albedo,                                &
  coszin, lit, solar_constant, scs, cos_zen_sph, day_frac_sph)


! Set the grid used by the core radiation code
CALL set_atm(control, dimen, spectrum, atm,                                    &
! Grid
  nlit, n_layer, nozone, nd_field, list, true_latitude, true_longitude,        &
! Thermodynamic fields
  p_layer_boundaries, p_layer_centres, t_layer_boundaries, t_layer_centres,    &
  p_extra_layer, t_extra_layer, d_mass, density,                               &
  r_layer_centres, r_layer_boundaries,                                         &
! Gas mass mixing ratios
  h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio,                   &
  nullmmr, nullmmr, o2_mix_ratio,                                              &
  nullmmr, nullmmr, nullmmr, nullmmr, nullmmr,                                 &
! 3D CO2
  co2_dim1, co2_dim2, co2_3d, l_co2_3d,                                        &
! Chemical greenhouse gas fields
  ngrgas, grgas_field,                                                         &
! Fields to calculate satellite viewing directions
  sindec, seconds_since_midnight)


DO l=1, nlit
  land0p5_g(l)=land0p5(list(l))
  flandg_g(l)=flandg(list(l))
  alat(l) = recip_pi_over_180 * true_latitude(list(l))
END DO



! Assign the properties of clouds.
IF (control%l_cloud) THEN

  ! check the consistency of cloud diagnostics.
  IF (sw_diag%re_conv_flag) THEN
    IF (.NOT. sw_diag%wgt_conv_flag) THEN
      cmessage =                                                             &
        '*** error: microphysical diagnostics for convective ' //            &
        'cloud must include the cloud weighting.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  END IF

  IF ( (sw_diag%re_strat_flag) .OR. (sw_diag%lwp_strat_flag) ) THEN
    IF (.NOT. sw_diag%wgt_strat_flag) THEN
      cmessage =                                                             &
        '*** error: microphysical diagnostics for stratiform ' //            &
        'cloud must include the cloud weighting.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  END IF

  IF (sw_diag%l_cloud_extinction) THEN
    IF (.NOT. sw_diag%l_cloud_weight_extinction) THEN
      cmessage =                                                             &
        '*** error: the cloud extinction ' //                                &
        'may be diagnosed only in conjunction ' //                           &
        'with the corresponding weights.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  END IF

  IF (sw_diag%l_ls_cloud_extinction) THEN
    IF (.NOT. sw_diag%l_ls_cloud_weight_extinction) THEN
      cmessage =                                                             &
        '*** error: the layer cloud extinction ' //                          &
        'may be diagnosed only in conjunction ' //                           &
        'with the corresponding weights.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  END IF

  IF (sw_diag%l_cnv_cloud_extinction) THEN
    IF (.NOT. sw_diag%l_cnv_cloud_weight_extinction) THEN
      cmessage =                                                             &
        '*** error: the conv. cloud extinction ' //                          &
        'may be diagnosed only in conjunction ' //                           &
        'with the corresponding weights.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  END IF


  CALL r2_set_cloud_parametrization(ierr, spectrum%basic%n_band              &
    , control%i_st_water, control%i_cnv_water                                &
    , control%i_st_ice, control%i_cnv_ice                                    &
    , spectrum%drop%l_drop_type, spectrum%drop%i_drop_parm                   &
    , spectrum%drop%n_phf, spectrum%drop%parm_list                           &
    , spectrum%drop%parm_min_dim, spectrum%drop%parm_max_dim                 &
    , spectrum%ice%l_ice_type, spectrum%ice%i_ice_parm                       &
    , spectrum%ice%n_phf, spectrum%ice%parm_list                             &
    , spectrum%ice%parm_min_dim, spectrum%ice%parm_max_dim                   &
    , cld%i_condensed_param, cld%condensed_n_phf                             &
    , cld%condensed_param_list, condensed_min_dim, condensed_max_dim         &
    , spectrum%dim%nd_band                                                   &
    , spectrum%dim%nd_drop_type, spectrum%dim%nd_ice_type                    &
    , spectrum%dim%nd_cloud_parameter, dimen%nd_cloud_component)
  IF (ierr /= i_normal) THEN
    cmessage = 'Error following call to r2_set_cloud_parametrization'
    GO TO 9999
  END IF

  CALL r2_set_cloud_field(nlit, n_layer, nclds, list                         &
    , atm%p, atm%t, atm%mass, alat                                           &
    , ccb, cct, cca, cccwp, ccw, lcbase                                      &
    , lccwc1, lccwc2, lca_area, lca_bulk, n_drop_pot                         &
    , control%l_microphysics, control%l_aerosol_ccn                          &
    , sea_salt_film, sea_salt_jet                                            &
    , l_use_seasalt_indirect, salt_dim1, salt_dim2, salt_dim3                &
    , l_use_biogenic, biogenic, biogenic_dim1, biogenic_dim2                 &
    , sulp_dim1, sulp_dim2, accum_sulphate, diss_sulphate                    &
    , aitken_sulphate, l_use_bmass_indirect                                  &
    , bmass_dim1, bmass_dim2, aged_bmass, cloud_bmass                        &
    , l_use_ocff_indirect, ocff_dim1, ocff_dim2                              &
    , aged_ocff, cloud_ocff                                                  &
    , l_use_nitrate_indirect, nitrate_dim1, nitrate_dim2                     &
    , accum_nitrate, diss_nitrate                                            &
    , l_easyaerosol_cdnc, easyaerosol_cdnc                                   &
    , lying_snow                                                             &
    , land0p5_g, flandg_g                                                    &
    , control%i_cloud_representation, cld%i_condensed_param                  &
    , condensed_min_dim, condensed_max_dim                                   &
    , cld%n_condensed, cld%type_condensed                                    &
    , cld%w_cloud, cld%n_cloud_type, cld%frac_cloud                          &
    , control%l_local_cnv_partition                                          &
    , cld%condensed_mix_ratio, cld%condensed_dim_char, CDNC                  &
    , sw_diag                                                                &
    , col_list, row_list, row_length, rows                                   &
    , nd_field, dimen%nd_profile, dimen%nd_layer                             &
    , spectrum%dim%nd_aerosol_species, dimen%nd_cloud_component              &
    , dimen%nd_cloud_type, dimen%id_cloud_top, n_cca_lev                     &
    )

  SELECT CASE (control%i_cloud)
  CASE (ip_cloud_mcica)
    i_cloud_tmp=ip_cloud_mix_max
  CASE (ip_cloud_part_corr)
    i_cloud_tmp=ip_cloud_mix_max
  CASE (ip_cloud_part_corr_cnv)
    i_cloud_tmp=ip_cloud_triple
  CASE DEFAULT
    i_cloud_tmp=control%i_cloud
  END SELECT

  l_all_temps = .TRUE.
  l_radius = .TRUE.
  IF (sw_diag%weighted_re_flag .AND. sw_diag%sum_weight_re_flag) THEN
    CALL r2_cloud_level_diag(control, dimen, atm, cld, CDNC                  &
      , nclds, list, i_cloud_tmp, l_all_temps, l_radius                      &
      , sw_diag%weighted_re, sw_diag%sum_weight_re                           &
      , col_list, row_list, row_length, rows, nd_field)
  END IF

  l_all_temps = .FALSE.
  l_radius = .TRUE.
  IF (sw_diag%wgtd_warm_re_flag .AND. sw_diag%sum_wgt_warm_re_flag) THEN
    CALL r2_cloud_level_diag(control, dimen, atm, cld, CDNC                  &
      , nclds, list, i_cloud_tmp, l_all_temps, l_radius                      &
      , sw_diag%weighted_warm_re, sw_diag%sum_weight_warm_re                 &
      , col_list, row_list, row_length, rows, nd_field)
  END IF

  l_all_temps = .TRUE.
  l_radius = .FALSE.
  IF (sw_diag%cdnc_ct_diag_flag .AND. sw_diag%cdnc_ct_weight_flag) THEN
    CALL r2_cloud_level_diag(control, dimen, atm, cld, CDNC                  &
      , nclds, list, i_cloud_tmp, l_all_temps, l_radius                      &
      , sw_diag%cdnc_ct_diag, sw_diag%cdnc_ct_weight                         &
      , col_list, row_list, row_length, rows, nd_field)
  END IF



  IF (control%i_cloud==ip_cloud_mcica) THEN
    CALL allocate_cld_mcica(cld, dimen, spectrum)
    DO ll=1, spectrum%dim%nd_k_term
      DO l=1, spectrum%basic%n_band
        cld%subcol_k(l,ll)=sw_subcol_k(l,ll)
      END DO
    END DO
    DO l=1,subcol_need
      cld%subcol_reorder(l)=sw_subcol_reorder(l)
    END DO
    CALL mcica_order(control, spectrum, cld)

    DO l=1,nlit
      i=((row_list(l)-1)*row_length)+col_list(l)
      cld%frac_cloudy(l)=frac_cloudy_full(i)
    END DO
    IF (ALLOCATED(cic_sub_full)) THEN
      DO lll=1,MIN(subcol_need, tot_subcol_gen)
        DO ll=dimen%id_cloud_top,n_layer
          DO l=1,nlit
            i=((row_list(l)-1)*row_length)+col_list(l)
            cld%c_sub(l,ll,lll,1)=clw_sub_full(i,ll,lll)
            cld%c_sub(l,ll,lll,2)=cic_sub_full(i,ll,lll)
          END DO
        END DO
      END DO
    ELSE
      DO lll=1,MIN(subcol_need, tot_subcol_gen)
        DO ll=dimen%id_cloud_top,n_layer
          DO l=1,nlit
            i=((row_list(l)-1)*row_length)+col_list(l)
            cld%c_sub(l,ll,lll,1)=clw_sub_full(i,ll,lll)
            cld%c_sub(l,ll,lll,2)=clw_sub_full(i,ll,lll)
          END DO
        END DO
      END DO
    END IF
  END IF

  ! Scale the condensed water contents to simulate
  ! inhomogeneities in the clouds.
  IF ( (l_inhom_cloud) .AND. ( control%i_cloud /= ip_cloud_mcica ) ) THEN
    IF ( ALLOCATED(cloud_inhom_param_full) ) THEN
      DO k = 1, dimen%nd_cloud_component
        DO j = dimen%id_cloud_top, n_layer
          DO i = 1, nlit
            l=((row_list(i)-1)*row_length)+col_list(i)
            cld%condensed_mix_ratio(i,j,k) = cloud_inhom_param_full(l,j)     &
              * cld%condensed_mix_ratio(i,j,k)
          END DO
        END DO
      END DO
    ELSE
      DO k = 1, dimen%nd_cloud_component
        DO j = dimen%id_cloud_top, n_layer
          DO i = 1, nlit
            cld%condensed_mix_ratio(i,j,k)                                   &
              = inhom_cloud(k) * cld%condensed_mix_ratio(i,j,k)
          END DO
        END DO
      END DO
    END IF
  END IF

  cld%dp_corr_strat = dp_corr_strat
  cld%dp_corr_conv  = dp_corr_conv

END IF ! control%l_cloud


CALL set_aer(control, dimen, spectrum, atm, cld, aer,                        &
! Model grid
  list,                                                                      &
! Properties of CLASSIC aerosols
  l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero,                    &
  bl_depth, n_levels_bl, l_murk_rad, aero_meso,                              &
  l_dust, l_use_dust, dust_dim1, dust_dim2,                                  &
  dust_1, dust_2, dust_3, dust_4, dust_5, dust_6,                            &
  l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic,                    &
  l_sulpc_so2, l_use_sulpc_direct, sulp_dim1, sulp_dim2,                     &
  accum_sulphate, aitken_sulphate,                                           &
  l_use_seasalt_direct, salt_dim1, salt_dim2, salt_dim3,                     &
  sea_salt_film, sea_salt_jet,                                               &
  l_soot, l_use_soot_direct, soot_dim1, soot_dim2,                           &
  fresh_soot, aged_soot,                                                     &
  l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2,                     &
  fresh_bmass, aged_bmass,                                                   &
  l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2,                           &
  fresh_ocff, aged_ocff,                                                     &
  l_nitrate, l_use_nitrate_direct, nitrate_dim1, nitrate_dim2,               &
  accum_nitrate,                                                             &
  n_arcl_species, n_arcl_compnts, i_arcl_compnts,                            &
  l_use_arcl, arcl_dim1, arcl_dim2, arcl,                                    &
  land0p5, lying_snow, pstar,                                                &
  p_layer_boundaries, trindx, alat, previous_time,                           &
! Properties of UKCA aerosols
  l_ukca_radaer, l_glomap_clim_radaer, ukca_radaer, ukca_dim1, ukca_dim2,    &
  ukca_mmr, ukca_cvl, ukca_dry, ukca_wet,                                    &
  ukca_rho, ukca_vol, ukca_wtv, ukca_nbr,                                    &
! Properties of EasyAerosol
  l_easyaerosol_rad, easyaerosol_rad,                                        &
! Diagnostics
  sw_diag, row_list, col_list,                                               &
! Dimensions of arrays
  nd_field, n_ukca_cpnt, n_ukca_mode )


9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE socrates_init
END MODULE socrates_init_mod
