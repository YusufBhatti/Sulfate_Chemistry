! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Longwave interface to the core radiation code.
!
! Method:
!   Principally, arrays are transferred to the appropriate formats.
!   Separate subroutines are called for each physical process.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE lw_rad_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'LW_RAD_MOD'
CONTAINS

SUBROUTINE lw_rad(ierr                                                         &
!                   Gaseous mixing ratios
    , h2o, co2, o3                                                             &
    , co2_dim1, co2_dim2, co2_3d, l_co2_3d                                     &
!                   Chemical greenhouse gas fields
    , ngrgas, grgas_field                                                      &
    , n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio                              &
    , cfc11_mix_ratio, cfc12_mix_ratio, cfc113_mix_ratio                       &
    , cfc114_mix_ratio                                                         &
    , hcfc22_mix_ratio, hfc125_mix_ratio, hfc134a_mix_ratio                    &
!                   Thermodynamic variables
    , t_layer_centres, t_rad_surf, t_rad_land, t_rad_sice, t_rad_sea           &
    , pstar, p_layer_boundaries, p_layer_centres                               &
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
    , land, flandg, ice_fraction                                               &
    , lying_snow, emis_land                                                    &
!                   Solar fields
    , coszin, lit, solar_constant, scs, sindec                                 &
!                   Aerosol fields
    , l_climat_aerosol, l_clim_aero_hgt, l_hadgem1_clim_aero                   &
    , bl_depth, n_levels_bl                                                    &
    , l_dust, l_use_dust, dust_dim1, dust_dim2                                 &
    , dust_1, dust_2, dust_3, dust_4, dust_5, dust_6                           &
    , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic                   &
    , l_sulpc_so2, l_use_sulpc_direct, l_use_sulpc_indirect                    &
    , sulp_dim1,sulp_dim2                                                      &
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
    , true_latitude, true_longitude, obs_solid_angle                           &
!                   Level of tropopause
    , trindx                                                                   &
!                   Spectrum
    , spectrum                                                                 &
!                   Algorithmic options
    , control, pts, list                                                       &
!                   Diagnostics
    , lw_diag, row_list, col_list                                              &
!                   Physical dimensions
    , dimen, n_points, n_layer, nclds                                          &
    , nozone, row_length, rows, nd_field                                       &
    , nd_field_flux_diag, nd_field_rad_diag                                    &
    , n_cca_lev, n_ukca_mode, n_ukca_cpnt                                      &
!                   Output fields
    , olr, lw_down, top_absorption, lwsea, lwout                               &
!                   COSP input arguments
    , cosp_gbx, cosp_sgx, cosp_sgh                                             &
    )

USE rad_pcf
USE conversions_mod, ONLY: recip_pi_over_180, pi
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_atm,      ONLY: StrAtm,                       deallocate_atm
USE def_cld,      ONLY: StrCld,   allocate_cld,       deallocate_cld,          &
                                  allocate_cld_prsc,  deallocate_cld_prsc,     &
                                  allocate_cld_mcica, deallocate_cld_mcica
USE def_aer,      ONLY: StrAer,                       deallocate_aer,          &
                                                      deallocate_aer_prsc
USE def_bound,    ONLY: StrBound,                     deallocate_bound
USE def_out,      ONLY: StrOut,                       deallocate_out
USE def_diag,     ONLY: StrDiag
USE mcica_mod
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ukca_radaer_struct_mod
USE coradoca,     ONLY: c2c_aerosol,                                           &
                        c2c_co2, co2_mmr_scl, co2_mmr_add,                     &
                        c2c_n2o, n2o_mmr_scl, n2o_mmr_add,                     &
                        c2c_ch4, ch4_mmr_scl, ch4_mmr_add,                     &
                        c2c_cfc11, cfc11_mmr_scl, cfc11_mmr_add,               &
                        c2c_cfc12, cfc12_mmr_scl, cfc12_mmr_add,               &
                        c2c_c113, cfc113_mmr_scl, cfc113_mmr_add,              &
                        c2c_hcfc22, hcfc22_mmr_scl, hcfc22_mmr_add,            &
                        c2c_hfc125, hfc125_mmr_scl, hfc125_mmr_add,            &
                        c2c_hfc134, hfc134a_mmr_scl, hfc134a_mmr_add
USE gas_list_pcf, ONLY: ip_co2, ip_n2o, ip_ch4, ip_cfc11, ip_cfc12,            &
                        ip_cfc113, ip_hcfc22, ip_hfc125, ip_hfc134a
USE ereport_mod, ONLY: ereport
USE arcl_mod, ONLY: npd_arcl_species, npd_arcl_compnts
USE cosp_constants_mod, ONLY: i_cvcliq, i_cvcice, i_lscliq, i_lscice
USE cosp_types_mod, ONLY: cosp_gridbox, cosp_subgrid, cosp_sghydro
USE rad_input_mod, ONLY: l_extra_top
USE nlsizes_namelist_mod, ONLY: model_levels
USE errormessagelength_mod, ONLY: errormessagelength
USE set_atm_mod, ONLY: set_atm
USE def_easyaerosol, ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

USE mcica_order_mod, ONLY: mcica_order
USE r2_calc_total_cloud_cover_mod, ONLY: r2_calc_total_cloud_cover
USE r2_set_cloud_field_mod, ONLY: r2_set_cloud_field
USE r2_set_cloud_parametrization_mod, ONLY: r2_set_cloud_parametrization
USE set_aer_mod, ONLY: set_aer
USE set_bound_mod, ONLY: set_bound
USE socrates_calc_mod, ONLY: socrates_calc
USE set_diag_mod, ONLY: set_diag

IMPLICIT NONE


! Dummy arguments

INTEGER ::                                                                     &
          !, intent(out)
    ierr
!       Error flag

! Dimensions of arrays:
INTEGER, INTENT(IN) :: row_length
!                          length of rows on each domain
INTEGER, INTENT(IN) :: rows
!                          number of rows in the domain
INTEGER ::                                                                     &
          !, intent(in)
    nd_field                                                                   &
!       Field size in calling program
    , nd_field_flux_diag                                                       &
!       Field size for flux diagnostics
    , nd_field_rad_diag
!       Field size for radiance diagnostics

! Actual sizes used:
INTEGER ::                                                                     &
          !, intent(in)
    n_points                                                                   &
!       Number of points
    , nozone                                                                   &
!       Number of levels with ozone
    , n_layer                                                                  &
!       Number of layers seen in the radiation scheme
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

! Gaseous mixing ratios:
REAL ::                                                                        &
          !, intent(in)
    h2o(nd_field, model_levels)                                                &
!       Mass mixing ratio of water
    , co2                                                                      &
!       Mass mixing ratio of CO2
    , o3(nd_field, nozone)                                                     &
!       Mass mixing ratios of ozone
    , n2o_mix_ratio                                                            &
!       Mass mixing ratio of nitrous oxide
    , ch4_mix_ratio                                                            &
!       Mass mixing ratio of methane
    , so2_mix_ratio                                                            &
!       Mass mixing ratio of sulphur dioxide
    , cfc11_mix_ratio                                                          &
!       Mass mixing ratio of CFC11
    , cfc12_mix_ratio                                                          &
!       Mass mixing ratio of CFC12
    , cfc113_mix_ratio                                                         &
!       Mass mixing ratio of CFC113
    , cfc114_mix_ratio                                                         &
!       Mass mixing ratio of CFC114
    , hcfc22_mix_ratio                                                         &
!       Mass mixing ratio of HCFC22
    , hfc125_mix_ratio                                                         &
!       Mass mixing ratio of HFC125
    , hfc134a_mix_ratio
!       Mass mixing ratio of HFC134a

INTEGER, INTENT(IN) :: ngrgas
REAL, INTENT(IN) :: grgas_field(nd_field, model_levels, ngrgas)

! General atmospheric properties:
REAL ::                                                                        &
          !, intent(in)
    p_layer_boundaries(nd_field,0:model_levels)                                &
!       Pressure at boundaries of layers
    , p_layer_centres(nd_field,0:model_levels)                                 &
!       Pressure at centres of layers
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
REAL, INTENT(IN) ::                                                            &
    coszin(nd_field)                                                           &
!       Cosines of zenith angle
    , solar_constant                                                           &
!       Total solar irradiance at 1 AU
    , scs                                                                      &
!       Scaling of solar incident field
    , lit(nd_field)                                                            &
!       Fraction of time point is lit
    , sindec
!       sin(solar declination)

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
!       Liquid water contents (these are not used directly in
!       the radiation: the total condensed water content is
!       repartitioned using focwwil).
    , lccwc2(nd_field, nclds+1/(nclds+1))                                      &
!       Ice water contents (these are not used directly in
!       The radiation: the total condensed water content is
!       Repartitioned using focwwil).
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
!       Fraction of grid-box covered by convective cloud

! Aerosols:
LOGICAL ::                                                                     &
          !, intent(in)
    l_climat_aerosol                                                           &
!       Flag for climatological aerosol
    , l_clim_aero_hgt                                                          &
!       Flag to use the depth of the boundary layer to set
!       the climatological aerosol
     , l_hadgem1_clim_aero                                                     &
!       Flag to use HadGEM1 setting for climatological aerosols
    , l_murk_rad
!       Flag for mesoscale model aerosol
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
          !,intent (in)
    sulp_dim1,sulp_dim2                                                        &
!       Dimensions for _sulphate arrays, (P_FIELD,P_LEVELS or 1,1)
    , dust_dim1, dust_dim2                                                     &
!       Dimensions for mineral dust arrays (p_field,p_levels or 1,1)
    , biogenic_dim1, biogenic_dim2                                             &
!       dimensions for biogenic array passed down to
!       r2_set_aerosol_field if direct effect required.
    , soot_dim1, soot_dim2                                                     &
!       dimensions for soot arrays (P_FIELD,P_LEVELS or 1,1)
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
!       Number concentration of film-mode sea-salt aerosol
    , sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                            &
!       Number concentration of jet-mode sea-salt aerosol
    , fresh_soot(soot_dim1, soot_dim2)                                         &
!       Mixing ratios of fresh soot
    , aged_soot(soot_dim1, soot_dim2)                                          &
!       Mixing ratios of aged soot
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
    , aero_meso(nd_field, model_levels)
!       Mixing ratio of 'urban' aerosol of mesoscale model

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
LOGICAL :: l_co2_3d    !controls use of 3d co2 field
INTEGER ::                                                                     &
          !, intent(in)
    co2_dim1, co2_dim2
!       Dimensions for CO2 array, (P_FIELD,P_LEVELS or 1,1)
REAL ::                                                                        &
          !, intent(in)
    co2_3d(co2_dim1, co2_dim2)
!       Mass mixing ratio of carbon dioxide
! Surface fields:
LOGICAL ::                                                                     &
          !, intent(in)
    land(nd_field)
!       Land mask (true if land fraction >0.5)
REAL ::                                                                        &
          !, intent(in)
    flandg(nd_field)
!       land fraction in grid box
REAL ::                                                                        &
          !, intent(in)
    pstar(nd_field)                                                            &
!       Surface pressures
    , t_rad_surf(nd_field)                                                     &
!       Effective radiative temperature over whole grid-box
    , t_rad_land(nd_field)                                                     &
!       Effective radiative temperature over land
    , t_rad_sice(nd_field)                                                     &
!       Effective radiative temperature over sea-ice
    , t_rad_sea(nd_field)                                                      &
!       Radiative temperature over open sea
    , ice_fraction(nd_field)                                                   &
!       Sea ice fraction of sea portion of grid box
    , lying_snow(nd_field)                                                     &
!       Mass loading of lying snow
    , emis_land(nd_field)                                                      &
!       Mean land emissivity in a gridbox if l_um_jules is true
    , bl_depth(nd_field)
!       depth of the boundary layer

REAL, INTENT(IN) :: true_latitude(nd_field)
REAL, INTENT(IN) :: true_longitude(nd_field)
REAL, INTENT(IN) :: seconds_since_midnight
INTEGER, INTENT(IN) :: previous_time(7)

! Solid angle subtended by grid-box for an observer at 1 AU
REAL, INTENT(IN) :: obs_solid_angle(nd_field)

!                   Level of tropopause
INTEGER, INTENT(IN) ::                                                         &
    trindx(nd_field)
!       The layer boundary of the tropopause

! increment of time:
REAL ::                                                                        &
          !, intent(in)
    pts
!       Time increment

INTEGER, INTENT(IN) :: list(nd_field)
!                          list of points where radiation is to be
!                          calculated

! Calculated fluxes:
REAL ::                                                                        &
          !, intent(out)
    olr(nd_field)                                                              &
!       Net outgoing radiation
    , lw_down(nd_field)                                                        &
!       Downwards surface flux
    , top_absorption(nd_field)                                                 &
!       Absorption in the extra radiative layer at the top
!       of the model
    , lwout(nd_field, model_levels+1)                                          &
!       Net downward fluxes or heating rates
    , lwsea(nd_field)
!       Sea-surface components of flux



! Diagnostics:
TYPE (StrDiag) :: lw_diag

INTEGER, INTENT(IN) :: row_list(nd_field)
!                          list of row indices of points
!                          to be treated
INTEGER, INTENT(IN) :: col_list(nd_field)
!                          list of column indices of points
!                          to be treated

! Structures with COSP arguments
TYPE(cosp_gridbox),INTENT(INOUT) :: cosp_gbx
TYPE(cosp_subgrid),INTENT(INOUT) :: cosp_sgx
TYPE(cosp_sghydro),INTENT(INOUT) :: cosp_sgh


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

! Controlling options for diagnostic calls:
TYPE (StrCtrl) :: control_clean
TYPE (StrCtrl) :: control_forc

! Atmospheric properties for diagnostic calls:
TYPE(StrAtm) :: atm_forc

! Output fields from diagnostic calls:
TYPE(StrOut) :: radout_clean
TYPE(StrOut) :: radout_forc

INTEGER :: i, j, k, l, ll, lll, ic
!       Loop variables

! General atmospheric properties:
REAL, PARAMETER :: nullmmr = 0.0
!       Null mass mixing ratio
REAL :: alat(dimen%nd_profile)
!       Latitude in degrees

! Surface fields:
LOGICAL :: land_g(dimen%nd_profile)
!       Gathered land-surface mask
REAL :: ice_fraction_g(dimen%nd_profile)
!       Gathered ice fraction
REAL :: flandg_g(dimen%nd_profile)
!       Gathered land fraction

! Array related to tiling of the surface
INTEGER ::                                                                     &
      index_tile(dimen%nd_tile_type)
!       The indexing number of tiles of the given type

! Cloudy properties:
REAL ::                                                                        &
    condensed_min_dim(dimen%nd_cloud_component)                                &
!       Minimum dimensions of condensed components
    , condensed_max_dim(dimen%nd_cloud_component)                              &
!       Maximum dimensions of condensed components
    , CDNC(dimen%nd_profile, dimen%id_cloud_top : dimen%nd_layer,              &
           dimen%nd_cloud_component)
!       CDNC array for call to r2_set_cloud_field from LW       

! Small real number used in COSP diagnostics calculations
REAL :: epsreal
! Temporary variables used for COSP diagnostics when McICA is used
REAL :: cosp_temp1, cosp_temp2
! Conversion factors for COSP effective radis of ice clouds
REAL :: f_re_lsice, f_re_cvice

! Fluxes:
REAL ::                                                                        &
    flux_net(dimen%nd_flux_profile, 0: dimen%nd_layer, dimen%nd_channel),      &
!       Net flux
    flux_net_clear(dimen%nd_flux_profile, 0: dimen%nd_layer, dimen%nd_channel)
!       Clear-sky net flux

! local arrays to fill diagnostics:
REAL :: total_cloud_cover_g(dimen%nd_profile)
!         cloud fraction at gathered points

! Variables required for compatibility with subroutines:
REAL :: dummy1d(1), dummy2d(1,1), dummy3d(1,1,1)

! To calculate water/ice paths 
REAL :: work_diag(dimen%nd_profile, model_levels)

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'LW_RAD'
CHARACTER (LEN=errormessagelength) :: cmessage

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
CALL set_bound(control, dimen, spectrum, bound, lw_diag,                       &
  nd_field, n_points, n_layer, list, row_list, col_list,                       &
  land, flandg, ice_fraction,                                                  &
  emis_land, t_rad_surf, t_rad_land, t_rad_sice, t_rad_sea,                    &
  dummy3d, dummy2d, dummy2d,                                                   &
  coszin, lit, solar_constant, scs, dummy2d, dummy2d)

! Set the clear-sky atmospheric profiles
CALL set_atm(control, dimen, spectrum, atm,                                  &
! Grid
  n_points, n_layer, nozone, nd_field, list, true_latitude, true_longitude,  &
! Thermodynamic fields
  p_layer_boundaries, p_layer_centres, t_layer_boundaries, t_layer_centres,  &
  p_extra_layer, t_extra_layer, d_mass, density,                             &
  r_layer_centres, r_layer_boundaries,                                       &
! Gas mass mixing ratios
  h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio,                 &
  cfc11_mix_ratio, cfc12_mix_ratio, nullmmr,                                 &
  cfc113_mix_ratio, cfc114_mix_ratio, hcfc22_mix_ratio,                      &
  hfc125_mix_ratio, hfc134a_mix_ratio,                                       &
! 3D CO2
  co2_dim1, co2_dim2, co2_3d, l_co2_3d,                                      &
! Chemical greenhouse gas fields
  ngrgas, grgas_field,                                                       &
! Fields to calculate satellite viewing directions
  sindec, seconds_since_midnight)


! Gathering of fields
DO l=1, n_points
  land_g(l)=land(list(l))
  flandg_g(l)=flandg(list(l))
  alat(l) = recip_pi_over_180 * true_latitude(list(l))
END DO


! Assign the properties of clouds. A dummy array must be passed
! for the microphysical diagnostics since they are not available
! through STASH in the long-wave.
IF (control%l_cloud) THEN

  ! Check the consistency of cloud diagnostics
  IF (lw_diag%l_cloud_absorptivity) THEN
    IF (.NOT. lw_diag%l_cloud_weight_absorptivity) THEN
      cmessage =                                                             &
        '*** Error: The cloud absorptivity ' //                              &
        'may be diagnosed only in conjunction ' //                           &
        'with the corresponding weights.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  END IF

  IF (lw_diag%l_ls_cloud_absorptivity) THEN
    IF (.NOT. lw_diag%l_ls_cloud_weight_absorptivity) THEN
      cmessage =                                                             &
        '*** Error: The layer cloud absorptivity ' //                        &
        'may be diagnosed only in conjunction ' //                           &
        'with the corresponding weights.'
      ierr=i_err_fatal
      GO TO 9999
    END IF
  END IF

  IF (lw_diag%l_cnv_cloud_absorptivity) THEN
    IF (.NOT. lw_diag%l_cnv_cloud_weight_absorptivity) THEN
      cmessage =                                                             &
        '*** Error: The conv. cloud absorptivity ' //                        &
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


  CALL r2_set_cloud_field(n_points, n_layer, nclds, list                     &
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
    , land_g, flandg_g                                                       &
    , control%i_cloud_representation, cld%i_condensed_param                  &
    , condensed_min_dim, condensed_max_dim                                   &
    , cld%n_condensed, cld%type_condensed                                    &
    , cld%w_cloud, cld%n_cloud_type, cld%frac_cloud                          &
    , control%l_local_cnv_partition                                          &
    , cld%condensed_mix_ratio, cld%condensed_dim_char, CDNC                  &
    , lw_diag                                                                &
    , col_list, row_list, row_length, rows                                   &
    , nd_field, dimen%nd_profile, dimen%nd_layer                             &
    , spectrum%dim%nd_aerosol_species, dimen%nd_cloud_component              &
    , dimen%nd_cloud_type, dimen%id_cloud_top, n_cca_lev                     &
    )


  IF (control%i_cloud == ip_cloud_mcica) THEN
    CALL allocate_cld_mcica(cld, dimen, spectrum)
    DO ll=1, spectrum%dim%nd_k_term
      DO l=1, spectrum%basic%n_band
        cld%subcol_k(l,ll)=lw_subcol_k(l,ll)
      END DO
    END DO
    DO l=1,subcol_need
      cld%subcol_reorder(l)=lw_subcol_reorder(l)
    END DO
    CALL mcica_order(control, spectrum, cld)

    DO l=1,n_points
      i=((row_list(l)-1)*row_length)+col_list(l)
      cld%frac_cloudy(l)=frac_cloudy_full(i)
      total_cloud_cover_g(l)=REAL(ncldy(i))/REAL(tot_subcol_gen)
    END DO
    IF (ALLOCATED(cic_sub_full)) THEN
      DO lll=1,MIN(subcol_need, tot_subcol_gen)
        DO ll=dimen%id_cloud_top,n_layer
          DO l=1,n_points
            i=((row_list(l)-1)*row_length)+col_list(l)
            cld%c_sub(l,ll,lll,1)=clw_sub_full(i,ll,lll)
            cld%c_sub(l,ll,lll,2)=cic_sub_full(i,ll,lll)
          END DO
        END DO
      END DO
    ELSE
      DO lll=1,MIN(subcol_need, tot_subcol_gen)
        DO ll=dimen%id_cloud_top,n_layer
          DO l=1,n_points
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
          DO i = 1, n_points
            l=((row_list(i)-1)*row_length)+col_list(i)
            cld%condensed_mix_ratio(i,j,k) = cloud_inhom_param_full(l,j)     &
              * cld%condensed_mix_ratio(i,j,k)
          END DO
        END DO
      END DO
    ELSE
      DO k = 1, dimen%nd_cloud_component
        DO j = dimen%id_cloud_top, n_layer
          DO i = 1, n_points
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
  land, lying_snow, pstar,                                                   &
  p_layer_boundaries, trindx, alat, previous_time,                           &
! Properties of UKCA aerosols
  l_ukca_radaer, l_glomap_clim_radaer, ukca_radaer, ukca_dim1, ukca_dim2,    &
  ukca_mmr, ukca_cvl, ukca_dry, ukca_wet,                                    &
  ukca_rho, ukca_vol, ukca_wtv, ukca_nbr,                                    &
! Properties of EasyAerosol
  l_easyaerosol_rad, easyaerosol_rad,                                        &
! Diagnostics
  lw_diag, row_list, col_list,                                               &
! Dimensions of arrays
  nd_field, n_ukca_cpnt, n_ukca_mode )


CALL socrates_calc(spectrum, control, lw_diag, dimen, atm, cld, aer, bound,    &
  radout, radout_clean, radout_forc)


! Output diagnostics common to SW and LW
CALL set_diag(n_points, n_layer, list, col_list, row_list, model_levels, &
  nd_field, obs_solid_angle, dummy1d, dummy1d, &
  control, atm, spectrum, radout, radout_clean, radout_forc, lw_diag)


! Processing depends on whether the code has been invoked to
! calculate radiances or fluxes.
IF ( (control%i_angular_integration == ip_two_stream) .OR.                     &
     ( (control%i_angular_integration == ip_spherical_harmonic) .AND.          &
        (control%i_sph_mode == ip_sph_mode_flux) ) ) THEN

  ! Convert downward fluxes to net fluxes.
  DO i=0, n_layer
    DO l=1, n_points
      flux_net(l, i, 1)=radout%flux_down(l, i, 1)-radout%flux_up(l, i, 1)
    END DO
  END DO
  IF (control%l_clear) THEN
    DO i=0, n_layer
      DO l=1, n_points
        flux_net_clear(l, i, 1)                                                &
          =radout%flux_down_clear(l, i, 1)-radout%flux_up_clear(l, i, 1)
      END DO
    END DO
  END IF


  ! OLR:
  IF ( control%l_solar_tail_flux ) THEN
    ! Solar tail flux needs to be removed for consistency
    ! with normal OLR
    DO l=1, n_points
      olr(list(l)) = radout%flux_up(l, 0, 1) -                                 &
                     radout%solar_tail_flux(l) / bound%zen_0(l)
    END DO
    IF (lw_diag%l_clear_olr) THEN
      DO l=1, n_points
        lw_diag%clear_olr(col_list(l), row_list(l))                            &
          = radout%flux_up_clear(l, 0, 1) -                                    &
            radout%solar_tail_flux(l) / bound%zen_0(l)
      END DO
    END IF
  ELSE
    DO l=1, n_points
      olr(list(l))=-flux_net(l, 0, 1)
    END DO
    IF (lw_diag%l_clear_olr) THEN
      DO l=1, n_points
        lw_diag%clear_olr(col_list(l), row_list(l))                            &
          =-flux_net_clear(l, 0, 1)
      END DO
    END IF
  END IF


  ! Total cloud cover:
  SELECT CASE (control%i_cloud)
  CASE (ip_cloud_mix_max, ip_cloud_mix_random, ip_cloud_triple,                &
        ip_cloud_part_corr, ip_cloud_part_corr_cnv)
    DO l=1, n_points
      total_cloud_cover_g(l)=radout%tot_cloud_cover(l)
    END DO
  CASE (ip_cloud_mcica)
    DO l=1,n_points
      i=((row_list(l)-1)*row_length)+col_list(l)
      total_cloud_cover_g(l)=REAL(ncldy(i))/REAL(tot_subcol_gen)
    END DO
  CASE DEFAULT
    CALL r2_calc_total_cloud_cover(n_points, nclds, nclds,                     &
      control%i_cloud, cld%w_cloud, total_cloud_cover_g,                       &
      dimen%nd_profile, dimen%nd_layer)
  END SELECT

  IF (lw_diag%l_total_cloud_cover) THEN
    DO l=1, n_points
      lw_diag%total_cloud_cover(col_list(l), row_list(l))                      &
        =total_cloud_cover_g(l)
    END DO
  END IF


  ! Total clear area
  IF (lw_diag%l_total_clear_area) THEN
    DO l=1, n_points
      lw_diag%total_clear_area(col_list(l), row_list(l))                       &
        =1.0-total_cloud_cover_g(l)
    END DO
  END IF


  ! Clear-sky flux at TOA weighted by total clear area
  IF (lw_diag%l_toa_clear_weighted) THEN
    DO l=1, n_points
      lw_diag%toa_clear_weighted(col_list(l), row_list(l))                     &
        =radout%flux_up_clear(l, 0, 1) * (1.0-total_cloud_cover_g(l))
    END DO
  END IF


  ! Net flux at the tropopause:
  IF (lw_diag%l_net_flux_trop) THEN
    DO l=1, n_points
      lw_diag%net_flux_trop(col_list(l), row_list(l))                          &
        =flux_net(l, n_layer+1-trindx(list(l)), 1)
    END DO
  END IF


  ! Downward flux at the tropopause:
  IF (lw_diag%l_down_flux_trop) THEN
    DO l=1, n_points
      lw_diag%down_flux_trop(col_list(l), row_list(l))                         &
        =flux_net(l, n_layer+1-trindx(list(l)), 1)                             &
        +radout%flux_up(l, n_layer+1-trindx(list(l)), 1)
    END DO
  END IF


  ! Downward flux at the surface:
  DO l=1, n_points
    lw_down(list(l)) = flux_net(l, n_layer, 1)                                 &
                     + radout%flux_up(l, n_layer, 1)
  END DO


  ! Clear-sky downward flux at the surface:
  IF (lw_diag%l_surf_down_clr) THEN
    DO l=1, n_points
      lw_diag%surf_down_clr(col_list(l), row_list(l))                          &
        =flux_net_clear(l, n_layer, 1)                                         &
        +radout%flux_up_clear(l, n_layer, 1)
    END DO
  END IF


  ! Cloud absorptivity diagnostics
  IF (lw_diag%l_cloud_absorptivity) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%cloud_absorptivity(col_list(l),row_list(l),i)                  &
           =radout%cloud_absorptivity(l, n_layer+1-i)
        lw_diag%cloud_weight_absorptivity(col_list(l),                         &
          row_list(l), i)                                                      &
          =radout%cloud_weight_absorptivity(l, n_layer+1-i)
      END DO
    END DO
  END IF

  IF (lw_diag%l_ls_cloud_absorptivity) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%ls_cloud_absorptivity(col_list(l),row_list(l),i)               &
           =radout%ls_cloud_absorptivity(l, n_layer+1-i)
        lw_diag%ls_cloud_weight_absorptivity(col_list(l),                      &
           row_list(l), i)                                                     &
           =radout%ls_cloud_weight_absorptivity(l, n_layer+1-i)
      END DO
    END DO
  END IF

  IF (lw_diag%l_cnv_cloud_absorptivity) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%cnv_cloud_absorptivity(col_list(l),row_list(l),i)              &
           =radout%cnv_cloud_absorptivity(l, n_layer+1-i)
        lw_diag%cnv_cloud_weight_absorptivity(col_list(l),                     &
           row_list(l), i)                                                     &
           =radout%cnv_cloud_weight_absorptivity(l, n_layer+1-i)
      END DO
    END DO
  END IF

  ! Total cloud fraction on model levels
  IF (lw_diag%l_total_cloud_on_levels) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%total_cloud_on_levels(col_list(l),row_list(l),i)               &
           =cld%w_cloud(l, n_layer+1-i)
      END DO
    END DO
  END IF


  ! Cloud water mixing ratios
  ! =========================
  
  ! LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
  IF (lw_diag%l_ls_qcl_rad .OR. lw_diag%l_ls_qcl_rad_path) THEN
    work_diag(:,:) = 0.0
    IF (control%i_cloud == ip_cloud_mcica) THEN
      DO i=1, nclds
        DO l=1, n_points
          lll=MIN( subcol_need,                                                &
                   ncldy((row_list(l)-1)*row_length+col_list(l)) )
          IF (lll > 0) THEN
            work_diag(l,i) =                                                   &
              SUM(cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_sw)               &
              *cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_water)        &
              *cld%c_sub(l,n_layer+1-i,1:lll,ip_cloud_type_sw))                &
              *total_cloud_cover_g(l)/REAL(lll)
          END IF
        END DO
      END DO
    ELSE
      DO i=1, nclds
        DO l=1, n_points
          work_diag(l,i) =                                                     &
              cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_water)         &
              *cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_sw)                  &
              *cld%w_cloud(l,n_layer+1-i)
        END DO
      END DO
    END IF

    ! Copy work array to diagnostic when requested 
    IF (lw_diag%l_ls_qcl_rad) THEN
      DO i=1, nclds
        DO l=1, n_points
          lw_diag%ls_qcl_rad(col_list(l),row_list(l),i) = work_diag(l,i)
        END DO
      END DO
    END IF

    ! Large scale liquid water path
    IF (lw_diag%l_ls_qcl_rad_path) THEN
      DO i = 1, nclds
        DO l=1, n_points   
          lw_diag%ls_qcl_rad_path(col_list(l),row_list(l)) =         & 
              lw_diag%ls_qcl_rad_path(col_list(l),row_list(l))       &
              + work_diag(l,i) * atm%mass(l,n_layer+1-i)
        END DO
      END DO
    END IF
  END IF

  ! LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (ICE)
  IF (lw_diag%l_ls_qcf_rad .OR. lw_diag%l_ls_qcf_rad_path) THEN
    work_diag(:,:) = 0.0  
    IF (control%i_cloud == ip_cloud_mcica) THEN
      DO i=1, nclds
        DO l=1, n_points
          lll=MIN( subcol_need,                                              &
                   ncldy((row_list(l)-1)*row_length+col_list(l)) )
          IF (lll > 0) THEN
            work_diag(l,i) =                                                 &
              SUM(cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_si)             &
              *cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_ice)        &
              *cld%c_sub(l,n_layer+1-i,1:lll,ip_cloud_type_si))              &
              *total_cloud_cover_g(l)/REAL(lll)
          END IF
        END DO
      END DO
    ELSE
      DO i=1, nclds
        DO l=1, n_points
          work_diag(l,i) =                                                   &
              cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_ice)         &
             *cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_si)                 &
             *cld%w_cloud(l,n_layer+1-i)
        END DO
      END DO
    END IF

    ! Copy work array to diagnostic when requested
    IF (lw_diag%l_ls_qcf_rad) THEN
      DO i=1, nclds
        DO l=1, n_points
          lw_diag%ls_qcf_rad(col_list(l),row_list(l),i) = work_diag(l,i)
        END DO
      END DO
    END IF

    ! Large scale ice water path
    IF ((lw_diag%l_ls_qcf_rad_path)) THEN
      DO i = 1, nclds
        DO l=1, n_points
          lw_diag%ls_qcf_rad_path(col_list(l),row_list(l)) =                 & 
                  lw_diag%ls_qcf_rad_path(col_list(l),row_list(l))           &
                + work_diag(l,i) * atm%mass(l,n_layer+1-i)
        END DO
      END DO
    END IF
  END IF

  ! CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
  IF (lw_diag%l_cc_qcl_rad .OR. lw_diag%l_cc_qcl_rad_path) THEN
    work_diag(:,:) = 0.0
    DO i=1, nclds
      DO l=1, n_points
        work_diag(l,i) =                                                   &
           cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_cnv_water)       &
           * cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_cw)                &
           * cld%w_cloud(l,n_layer+1-i)
      END DO
    END DO

    ! Copy work array to diagnostic when requested
    IF (lw_diag%l_cc_qcl_rad) THEN
      DO i=1, nclds
        DO l=1, n_points
          lw_diag%cc_qcl_rad(col_list(l),row_list(l),i) = work_diag(l,i)
        END DO
      END DO
    END IF
  
    ! Convective  liquid water path
    IF ((lw_diag%l_cc_qcl_rad_path)) THEN
      DO i = 1, nclds
        DO l=1, n_points
          lw_diag%cc_qcl_rad_path(col_list(l),row_list(l)) =               & 
            lw_diag%cc_qcl_rad_path(col_list(l),row_list(l)) +             &
            work_diag(l,i) * atm%mass(l,n_layer+1-i)
        END DO
      END DO
    END IF
  END IF

  ! CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (ICE)
  IF (lw_diag%l_cc_qcf_rad .OR. lw_diag%l_cc_qcf_rad_path) THEN
    work_diag(:,:) = 0.0
    DO i=1, nclds
      DO l=1, n_points
        work_diag(l,i) =                                                   &
           cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_cnv_ice)         &
           * cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_ci)                &
           * cld%w_cloud(l,n_layer+1-i)
      END DO
    END DO

    ! Copy work array to diagnostic when requested
    IF (lw_diag%l_cc_qcf_rad) THEN
      DO i=1, nclds
        DO l=1, n_points
          lw_diag%cc_qcf_rad(col_list(l),row_list(l),i) = work_diag(l,i)
        END DO
      END DO
    END IF

    ! Convective ice water path
    IF ((lw_diag%l_cc_qcf_rad_path)) THEN
      DO i = 1, nclds
        DO l=1, n_points
          lw_diag%cc_qcf_rad_path(col_list(l),row_list(l)) =               & 
            lw_diag%cc_qcf_rad_path(col_list(l),row_list(l)) +             &
            work_diag(l,i) * atm%mass(l,n_layer+1-i)
        END DO
      END DO
    END IF
  END IF


  ! Cloud fractions
  ! ===============

  ! LARGE-SCALE cloud GRIDBOX FRACTION seen by radiation. (LIQUID)
  IF (lw_diag%l_ls_cl_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%ls_cl_rad(col_list(l),row_list(l),i)                           &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_sw)                    &
           * cld%w_cloud(l,n_layer+1-i)
      END DO
    END DO
  END IF

  ! LARGE-SCALE cloud GRIDBOX fraction seen by radiation. (ICE)
  IF (lw_diag%l_ls_cf_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%ls_cf_rad(col_list(l),row_list(l),i)                           &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_si)                    &
           * cld%w_cloud(l,n_layer+1-i)
      END DO
    END DO
  END IF

  ! CONVECTIVE cloud GRIDBOX fraction seen by radiation. (LIQUID)
  IF (lw_diag%l_cc_cl_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%cc_cl_rad(col_list(l),row_list(l),i)                           &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_cw)                    &
           * cld%w_cloud(l,n_layer+1-i)
      END DO
    END DO
  END IF

  ! CONVECTIVE cloud GRIDBOX FRACTION seen by radiation. (ICE)
  IF (lw_diag%l_cc_cf_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%cc_cf_rad(col_list(l),row_list(l), i)                          &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_ci)                    &
           * cld%w_cloud(l,n_layer+1-i)
      END DO
    END DO
  END IF

  ! Weighted cloud particles effective dimensions
  ! =============================================
  ! Large-scale liquid
  IF (lw_diag%l_ls_del_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%ls_del_rad(col_list(l),row_list(l),i)                          &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_sw)                    &
           * cld%w_cloud(l,n_layer+1-i)                                        &
           * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_water)
      END DO
    END DO
  END IF
  ! Large-scale ice
  IF (lw_diag%l_ls_def_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%ls_def_rad(col_list(l),row_list(l),i)                          &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_si)                    &
           * cld%w_cloud(l,n_layer+1-i)                                        &
           * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_ice)
      END DO
    END DO
  END IF
  ! Convective liquid
  IF (lw_diag%l_cc_del_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%cc_del_rad(col_list(l),row_list(l),i)                          &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_cw)                   &
           * cld%w_cloud(l,n_layer+1-i)                                        &
           * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_cnv_water)
      END DO
    END DO
  END IF
  ! Convective ice
  IF (lw_diag%l_cc_def_rad) THEN
    DO i=1, nclds
      DO l=1, n_points
        lw_diag%cc_def_rad(col_list(l),row_list(l),i)                          &
           = cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_ci)                    &
           * cld%w_cloud(l,n_layer+1-i)                                        &
           * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_cnv_ice)
      END DO
    END DO
  END IF

  ! COSP arguments
  IF (l_cosp) THEN
    epsreal = EPSILON(1.0)
    ! The ice effective dimension for most parametrisations is equivalent
    ! to the diameter rather than the radius. Here we choose the adequate
    ! multiplicative factor to calculate the effective radius.
    f_re_lsice = 0.5
    f_re_cvice = 0.5
    IF (cld%i_condensed_param(ip_clcmp_st_ice) == ip_slingo_schrecker_ice)     &
      f_re_lsice = 1.0
    IF (cld%i_condensed_param(ip_clcmp_cnv_ice) == ip_slingo_schrecker_ice)    &
      f_re_cvice = 1.0
    IF (control%i_cloud == IP_cloud_mcica .AND.                             &
     cosp_sgx%Ncolumns == tot_subcol_gen) THEN
     DO i=1, nclds
       DO l=1, n_points
         j = (row_list(l)-1)*row_length + col_list(l)
         cosp_temp1 =                                                          &
             cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_water)*         &
             cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_sw)
         cosp_temp2 =                                                          &
             cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_ice)*           &
             cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_si)
         DO lll=1, ncldy(j)
           ! Cloud water SUBGRID mixing ratio (LIQUID)
           cosp_sgh%mr_hydro(j,lll,i,i_lscliq) =                               &
             clw_sub_full(j,n_layer+1-i,lll)*cosp_temp1
           ! Cloud water SUBGRID mixing ratio (ICE)
           cosp_sgh%mr_hydro(j,lll,i,i_lscice) =                               &
             clw_sub_full(j,n_layer+1-i,lll)*cosp_temp2
           ! Cloud water effective radius (LIQUID)
           cosp_sgh%Reff(j,lll,i,i_lscliq) =                                   &
             cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_water)
           ! Cloud water effective dimension (ICE)
           cosp_sgh%Reff(j,lll,i,i_lscice) =                                   &
             f_re_lsice * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_ice)
         END DO
         ! Sub-grid cloud emissivity for COSP assumes
         ! same generated sub-columns for liquid and ice.
         IF (radout%ls_cloud_weight_absorptivity(l, n_layer+1-i)               &
               > epsreal) THEN
           cosp_temp1 = atm%mass(l,n_layer+1-i)*                               &
             radout%ls_cloud_absorptivity(l, n_layer+1-i) /                    &
             radout%ls_cloud_weight_absorptivity(l, n_layer+1-i)
           DO lll=1, ncldy(j)
             cosp_sgx%dem(j,lll,i) =  1.0-EXP(-1.666*cosp_temp1*               &
               clw_sub_full(j,n_layer+1-i,lll))
           END DO
         END IF
         ! TOTAL cloud GRIDBOX fraction seen by radiation
         cosp_gbx%tca(j,i) = cld%w_cloud(l,n_layer+1-i)
         ! GRIDBOX effective radii
         ! Cloud water effective radius (LIQUID)
         cosp_gbx%reff(j,i,i_lscliq) =                                         &
           cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_water)
         ! Cloud water effective dimension (ICE)
         cosp_gbx%reff(j,i,i_lscice) =                                         &
           f_re_lsice * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_ice)
       END DO
     END DO
    ELSE
     DO i=1, nclds
      DO l=1, n_points
        j = (row_list(l)-1)*row_length + col_list(l)
        ! LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
        cosp_gbx%mr_hydro(j,i,i_lscliq) =                                      &
           cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_water)            &
           * cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_sw)                    &
           * cld%w_cloud(l,n_layer+1-i)
        ! LARGE-SCALE cloud water GRIDBOX MEAN mixing ratio (ICE)
        cosp_gbx%mr_hydro(j,i,i_lscice) =                                      &
           cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_st_ice)              &
           * cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_si)                    &
           * cld%w_cloud(l,n_layer+1-i)
        ! CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (LIQUID)
        cosp_gbx%mr_hydro(j,i,i_cvcliq) =                                      &
           cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_cnv_water)           &
           * cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_cw)                    &
           * cld%w_cloud(l,n_layer+1-i)
        ! CONVECTIVE cloud water GRIDBOX MEAN mixing ratio (ICE)
        cosp_gbx%mr_hydro(j,i,i_cvcice) =                                      &
           cld%condensed_mix_ratio(l,n_layer+1-i,ip_clcmp_cnv_ice)             &
           * cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_ci)                    &
           * cld%w_cloud(l,n_layer+1-i)
        ! TOTAL cloud GRIDBOX fraction seen by radiation
        cosp_gbx%tca(j,i) = cld%w_cloud(l,n_layer+1-i)
        ! CONVECTIVE cloud GRIDBOX fraction seen by radiation
        cosp_gbx%cca(j,i) =                                                    &
           (cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_cw) +                   &
            cld%frac_cloud(l,n_layer+1-i,ip_cloud_type_ci))                    &
           * cld%w_cloud(l,n_layer+1-i)
        ! LARGE-SCALE cloud water effective radius
        cosp_gbx%reff(j,i,i_lscliq) =                                          &
           cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_water)
        ! LARGE-SCALE cloud ice effective radius
        cosp_gbx%reff(j,i,i_lscice) =                                          &
           f_re_lsice * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_st_ice)
        ! CONVECTIVE cloud water effective radius
        cosp_gbx%reff(j,i,i_cvcliq) =                                          &
           cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_cnv_water)
        ! CONVECTIVE cloud ice effective radius
        cosp_gbx%reff(j,i,i_cvcice) =                                          &
           f_re_lsice * cld%condensed_dim_char(l,n_layer+1-i,ip_clcmp_cnv_ice)
        ! Large-Scale and convective cloud emissivities
        IF (radout%ls_cloud_weight_absorptivity(l, n_layer+1-i)                &
            > epsreal) THEN
          cosp_gbx%dem_s(j,i) = 1.0-EXP(-1.666*atm%mass(l,n_layer+1-i)         &
                * radout%ls_cloud_absorptivity(l, n_layer+1-i)                 &
                / radout%ls_cloud_weight_absorptivity(l, n_layer+1-i))
        END IF
        IF (radout%cnv_cloud_weight_absorptivity(l, n_layer+1-i)               &
            > epsreal) THEN
          cosp_gbx%dem_c(j,i) = 1.0-EXP(-1.666*atm%mass(l,n_layer+1-i)         &
                * radout%cnv_cloud_absorptivity(l, n_layer+1-i)                &
                / radout%cnv_cloud_weight_absorptivity(l, n_layer+1-i))
        END IF
      END DO
     END DO
    END IF
  END IF


  ! Output arrays:

  ! Convert the fluxes to increments in the heating rate except at
  ! the surface: there, the net downward flux is assigned to LWOUT.
  DO k=model_levels, 1, -1
    DO l=1, n_points
      lwout(list(l), k+1)=(flux_net(l, n_layer-k, 1)                           &
        -flux_net(l, n_layer+1-k, 1))                                          &
        *pts/layer_heat_capacity(list(l), k)
    END DO
    IF (lw_diag%l_clear_hr) THEN
      ! The factor of PTS is included here to yield a rate from an
      ! increment.
      DO l=1, n_points
        lw_diag%clear_hr(col_list(l), row_list(l), k)                          &
          =(flux_net_clear(l, n_layer-k, 1)                                    &
          -flux_net_clear(l, n_layer+1-k, 1))                                  &
          /layer_heat_capacity(list(l), k)
      END DO
    END IF
  END DO

  IF (l_extra_top) THEN
    ! calculate the radiation absorbed in the extra layer
    ! above the top of the rest of the model.
    DO l=1, n_points
      top_absorption(list(l))=flux_net(l, 0, 1)                                &
        -flux_net(l, n_layer-model_levels, 1)
    END DO
  END IF

  DO l=1, n_points
    lwout(list(l), 1)=flux_net(l, n_layer, 1)
  END DO

  ! Separate the contributions over open sea and sea-ice.
  ! Fluxes returned from the radiation code itself are not
  ! weighted by the fraction of the tile, but here are converted
  ! to grid-box mean values. This split is possible only if the
  ! ocean surface has been tiled.

  IF (control%l_tile) THEN

    ! The variable flandg is set even if coastal tiling is not
    ! used, so fairly generic code can be written.

    ! It is simplest to zero LWsea at all points and reset
    ! elsewhere.
    lwsea(list(1:n_points)) = 0.0

    DO ll=1, n_points
      l=list(ll)
      IF ( (flandg(l) < TINY(flandg)) .AND.                                    &
           (ice_fraction(l) < TINY(ice_fraction) )                             &
         ) THEN
        ! This point is open sea with no sea ice.
        lwsea(l)=lwout(l, 1)
        lwout(l, 1)=0.0
      END IF
    END DO

    index_tile(ip_ocean_tile) = 1

    ! Tiled points will have both land and sea. Note that the
    ! channel index of flux_up_tile is hard-wired to 1 because
    ! we don't envisage calling the code in other cases.
    DO lll=1, bound%n_point_tile
      ll=bound%list_tile(lll)
      l=list(ll)
      lwsea(l)=(1.0-ice_fraction(l))*(flux_net(ll, n_layer, 1)                 &
        +radout%flux_up(ll, n_layer, 1)                                        &
        -radout%flux_up_tile(lll, index_tile(ip_ocean_tile), 1))
      lwout(l, 1)=lwout(l, 1)-(1.0-flandg(l))*lwsea(l)
    END DO

    ! The remaining points are entirely land points and lwout
    ! need not be altered.

  ELSE

    ! Without radiative tiling we must assume that fluxes are
    ! uniform across the grid-box.
    WHERE (flandg(list(1:n_points)) < 1.0-TINY(flandg))
      lwsea(list(1:n_points))=(1.0-ice_fraction(list(1:n_points)))             &
        *flux_net(1:n_points, n_layer, 1)
      lwout(list(1:n_points), 1)=lwout(list(1:n_points), 1)                    &
        -(1.0-flandg(list(1:n_points)))*lwsea(list(1:n_points))
    END WHERE

  END IF

END IF

CALL deallocate_out(radout)
CALL deallocate_aer_prsc(aer)
CALL deallocate_aer(aer)
CALL deallocate_cld_mcica(cld)
CALL deallocate_bound(bound)
CALL deallocate_cld_prsc(cld)
CALL deallocate_cld(cld)
CALL deallocate_atm(atm)

9999 CONTINUE
! Check error condition
IF (ierr /= i_normal) THEN
  CALL ereport(RoutineName, ierr, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE lw_rad
END MODULE lw_rad_mod
