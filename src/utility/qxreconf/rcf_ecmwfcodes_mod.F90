! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  List of ECMWF magic numbers

MODULE Rcf_ECMWFcodes_Mod

! Description:
!    Magic numbers for ECMWF [dia|pro]gnostic codes, as defined by
!    the ECMWF version of 'Table 2' used in the GRIB description.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: Fortran 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE
PUBLIC

! ECMWF codes used in the reconfiguration

! Sea-ice Fraction
INTEGER, PARAMETER      :: ECMWFcode_icefrac             = 31

! Snow Density
INTEGER, PARAMETER      :: ECMWFcode_snow_density        = 33

! Volumetric Soil Water Layer 1
INTEGER, PARAMETER      :: ECMWFcode_soil_moist_1        = 39

! Volumetric Soil Water Layer 2
INTEGER, PARAMETER      :: ECMWFcode_soil_moist_2        = 40

! Volumetric Soil Water Layer 3
INTEGER, PARAMETER      :: ECMWFcode_soil_moist_3        = 41

! Volumetric Soil Water Layer 4
INTEGER, PARAMETER      :: ECMWFcode_soil_moist_4        = 42

! Specific rain water content
INTEGER, PARAMETER      :: ECMWFcode_qrain               = 75

! Specific snow water content
INTEGER, PARAMETER      :: ECMWFcode_qsnow               = 76

! Atmospheric Tide
INTEGER, PARAMETER      :: ECMWFcode_atm_tide            = 127

! Budget Values
INTEGER, PARAMETER      :: ECMWFcode_budget              = 128

! Geopotential
INTEGER, PARAMETER      :: ECMWFcode_geopot              = 129

! Temperature
INTEGER, PARAMETER      :: ECMWFcode_T                   = 130

! U-component of Wind
INTEGER, PARAMETER      :: ECMWFcode_u                   = 131

! V-component of Wind
INTEGER, PARAMETER      :: ECMWFcode_v                   = 132

! Specific Humidity
INTEGER, PARAMETER      :: ECMWFcode_spec_hum            = 133

! Surface Pressure
INTEGER, PARAMETER      :: ECMWFcode_pstar               = 134

! Vertical Velocity
INTEGER, PARAMETER      :: ECMWFcode_w                   = 135

! Precipitable Water Content
INTEGER, PARAMETER      :: ECMWFcode_precip_water        = 137

! Vorticity
INTEGER, PARAMETER      :: ECMWFcode_vort                = 138

! Soil Temperature level 1
INTEGER, PARAMETER      :: ECMWFcode_soil_temp_1         = 139

! Soil Wetness level 1
INTEGER, PARAMETER      :: ECMWFcode_soil_wet_1          = 140

! Snow Depth
INTEGER, PARAMETER      :: ECMWFcode_snow_depth          = 141

! Large Scale Precipitation
INTEGER, PARAMETER      :: ECMWFcode_lge_precip          = 142

! Convective Precipitation
INTEGER, PARAMETER      :: ECMWFcode_conv_precip         = 143

! Snow Fall
INTEGER, PARAMETER      :: ECMWFcode_snowfall            = 144

! Boundary Layer Dissipation
INTEGER, PARAMETER      :: ECMWFcode_b_layer_diss        = 145

! Surface Flux of Sensible Heat
INTEGER, PARAMETER      :: ECMWFcode_s_flux_s_heat       = 146

! Surface Flux of Latent Heat
INTEGER, PARAMETER      :: ECMWFcode_s_flux_l_heat       = 147

! Mean Sea Level Pressure
INTEGER, PARAMETER      :: ECMWFcode_msl_press           = 148

! Log Surface Pressure
INTEGER, PARAMETER      :: ECMWFcode_log_p               = 152

! Divergence
INTEGER, PARAMETER      :: ECMWFcode_div                 = 153

! Height (Geopotential)
INTEGER, PARAMETER      :: ECMWFcode_h_geopot            = 156

! Relative Humidity
INTEGER, PARAMETER      :: ECMWFcode_rel_hum             = 157

! Tendency of Surface Pressure
INTEGER, PARAMETER      :: ECMWFcode_tend_p              = 158

! Total Cloud Cover
INTEGER, PARAMETER      :: ECMWFcode_t_cloud_cover       = 159

! U-wind at 10m
INTEGER, PARAMETER      :: ECMWFcode_u_10                = 160

! V-wind at 10m
INTEGER, PARAMETER      :: ECMWFcode_v_10                = 161

! Temperature at 2m
INTEGER, PARAMETER      :: ECMWFcode_temp_2              = 167

! Dewpoint at 2m
INTEGER, PARAMETER      :: ECMWFcode_dew_p_2             = 168

! Soil Temperature level 2
INTEGER, PARAMETER      :: ECMWFcode_soil_temp_2         = 170

! Soil Wetness level 2
INTEGER, PARAMETER      :: ECMWFcode_soil_wet_2          = 171

! Land-Sea Mask
INTEGER, PARAMETER      :: ECMWFcode_lsm                 = 172

! Surface Roughness
INTEGER, PARAMETER      :: ECMWFcode_surface_rough       = 173

! Albedo
INTEGER, PARAMETER      :: ECMWFcode_albedo              = 174

! Downwards Shortwave Radiation (surface)
INTEGER, PARAMETER      :: ECMWFcode_down_swave_surface  = 175

! Net Shortwave Radiation (surface)
INTEGER, PARAMETER      :: ECMWFcode_net_swave_surface   = 176

! Net Longwave Radiation (surface)
INTEGER, PARAMETER      :: ECMWFcode_net_lwave_surface   = 177

! Net Shortwave Radiation (top of atmosphere)
INTEGER, PARAMETER      :: ECMWFcode_net_swave_top       = 178

! Net Longwave Radiation (top of atmosphere)
INTEGER, PARAMETER      :: ECMWFcode_net_lwave_top       = 179

! U-component of Surface Wind Stress
INTEGER, PARAMETER      :: ECMWFcode_u_stress            = 180

! V-component of Surface Wind Stress
INTEGER, PARAMETER      :: ECMWFcode_v_stress            = 181

!Evaporation
INTEGER, PARAMETER      :: ECMWFcode_evap                = 182

! Soil Temperature level 3
INTEGER, PARAMETER      :: ECMWFcode_soil_temp_3         = 183

! Soil Wetness level 3
INTEGER, PARAMETER      :: ECMWFcode_soil_wet_3          = 184

! Convective Cloud Cover
INTEGER, PARAMETER      :: ECMWFcode_conv_cloud_cover    = 185

! Low Cloud Cover
INTEGER, PARAMETER      :: ECMWFcode_low_cloud_cover     = 186

! Medium Cloud Cover
INTEGER, PARAMETER      :: ECMWFcode_med_cloud_cover     = 187

! High Cloud Cover
INTEGER, PARAMETER      :: ECMWFcode_high_cloud_cover    = 188

! Sunshine Duration
INTEGER, PARAMETER      :: ECMWFcode_sun_duration        = 189

! Ozone mass mixing ratio
INTEGER, PARAMETER      :: ECMWFcode_ozone               = 203

! Skin Temperature
INTEGER, PARAMETER      :: ECMWFcode_skin_temp           = 235

! Soil Temperature level 4
INTEGER, PARAMETER      :: ECMWFcode_soil_temp_4         = 236

! Soil Wetness level 4
INTEGER, PARAMETER      :: ECMWFcode_soil_wet_4          = 237

! Cloud Liquid Water Content
INTEGER, PARAMETER      :: ECMWFcode_qcl                 = 246

! Cloud Ice Water Content
INTEGER, PARAMETER      :: ECMWFcode_qcf                 = 247

! Cloud Cover
INTEGER, PARAMETER      :: ECMWFcode_cc                  = 248

! GEMS Hydrophobic Organic Matter Aerosol (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_OMFRSH              = 7

! GEMS Hydrophilic Organic Matter Aerosol (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_OMAGD               = 8

! GEMS Hydrophobic Black Carbon Aerosol (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_BCFRSH              = 9

! GEMS Hydrophilic Black Carbon Aerosol (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_BCAGD               = 10

! Methane (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_CH4                 = 62

! Nitrogen Oxides (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_NOX                 = 129

! GEMS NO2
INTEGER, PARAMETER      :: ECMWFcode_NO2                 = 121

! Carbon monoxide (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_CO                  = 123

! Formaldehyde (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_HCHO                = 124

! GEMS Ozone (GEMS Table 210)
INTEGER, PARAMETER      :: ECMWFcode_GO3                 = 203

END MODULE Rcf_ECMWFcodes_Mod
