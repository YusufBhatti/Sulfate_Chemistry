! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control

!- ---------------------------------------------------------------------
MODULE r2_set_cloud_field_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'R2_SET_CLOUD_FIELD_MOD'
CONTAINS

SUBROUTINE r2_set_cloud_field(n_profile, n_layer, nclds, i_gather &
   , p, t, d_mass, alat                                           &
   , ccb_in, cct_in, cca, cccwp, ccw, lcbase                      &
   , lccwc1, lccwc2, lca_area, lca_bulk, n_drop_pot               &
   , l_microphysics, l_aerosol_ccn                                &
   , sea_salt_film, sea_salt_jet                                  &
   , l_seasalt_ccn, salt_dim1, salt_dim2, salt_dim3               &
   , l_use_biogenic, biogenic, biogenic_dim1, biogenic_dim2       &
   , sulp_dim1, sulp_dim2, accum_sulphate, diss_sulphate          &
   , aitken_sulphate, l_biomass_ccn                               &
   , bmass_dim1, bmass_dim2, aged_bmass, cloud_bmass              &
   , l_ocff_ccn, ocff_dim1, ocff_dim2, aged_ocff, cloud_ocff      &
   , l_nitrate_ccn, nitrate_dim1, nitrate_dim2, accum_nitrate     &
   , diss_nitrate                                                 &
   , l_easyaerosol_cdnc, easyaerosol_cdnc                         &
   , lying_snow                                                   &
   , land_g, flandg_g                                             &
   , i_cloud_representation, i_condensed_param                    &
   , condensed_min_dim, condensed_max_dim                         &
   , n_condensed, type_condensed                                  &
   , w_cloud, n_cloud_type, frac_cloud, l_local_cnv_partition     &
   , condensed_mix_rat_area, condensed_dim_char, CDNC             &
   , diag                                                         &
   , col_list, row_list, row_length, rows                         &
   , npd_field, npd_profile, npd_layer, npd_aerosol_species       &
   , npd_cloud_component, npd_cloud_type, id_ct, n_cca_lev        &
   )

!  Subroutine to assign Properties of Clouds.

! Purpose:
!   The fractions of different types of clouds and their microphysical
!   properties are set.

! Method:
!   Straightforward.

USE rad_pcf
USE planet_constants_mod, ONLY: r
USE water_constants_mod, ONLY: tm
USE def_diag,   ONLY: StrDiag
USE cv_run_mod, ONLY: l_3d_cca

USE rad_input_mod, ONLY: l_rad_ovrlap, l_rad_ccw_scav, lrad_ccrad
USE bl_option_mod, ONLY: blending_option, off

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE cloud_inputs_mod, ONLY: l_add_cca_to_mcica,                  &
                             starticeTKelvin, alliceTdegC

USE nlsizes_namelist_mod, ONLY: model_levels
USE def_easyaerosol, ONLY: t_easyaerosol_cdnc

USE lsp_focwwil_mod, ONLY: lsp_focwwil
USE r2_re_mrf_umist_mod, ONLY: r2_re_mrf_umist
IMPLICIT NONE

!     DIMENSIONS OF ARRAYS:
INTEGER, INTENT(IN) :: row_length
!                              Number of grid-points in EW-direction
!                              in the local domain
INTEGER, INTENT(IN) :: rows
!                              Number of grid-points in NS-direction
!                              in the local domain
INTEGER ::                                                        &
          !, INTENT(IN)
     npd_field                                                    &
!             FIELD SIZE IN CALLING PROGRAM
         , npd_profile                                                  &
!             SIZE OF ARRAY OF PROFILES
         , npd_layer                                                    &
!             MAXIMUM NUMBER OF LAYERS
         , npd_aerosol_species                                          &
!             MAXIMUM NUMBER OF AEROSOL_SPECIES
         , npd_cloud_component                                          &
!             Number of components of clouds
         , npd_cloud_type                                               &
!             Number of permitted types of clouds.
         , id_ct                                                        &
!             Topmost declared cloudy layer
         , sulp_dim1                                                    &
!             1ST DIMENSION OF ARRAYS OF SULPHATE
         , sulp_dim2                                                    &
!             2ND DIMENSION OF ARRAYS OF SULPHATE
         , bmass_dim1                                                   &
!             1ST DIMENSION OF ARRAYS OF BIOMASS SMOKE
         , bmass_dim2                                                   &
!             2ND DIMENSION OF ARRAYS OF BIOMASS SMOKE
         , ocff_dim1                                                    &
!             1ST DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
         , ocff_dim2                                                    &
!             2ND DIMENSION OF ARRAYS OF FOSSIL-FUEL ORGANIC CARBON
         , salt_dim1                                                    &
!             1ST DIMENSION OF ARRAYS OF SEA-SALT
         , salt_dim2                                                    &
!             2ND DIMENSION OF ARRAYS OF SEA-SALT
         , salt_dim3                                                    &
!             3RD DIMENSION OF ARRAYS OF SEA-SALT
         , biogenic_dim1                                                &
!             1ST DIMENSION OF BIOGENIC AEROSOL ARRAY
         , biogenic_dim2                                                &
!             2ND DIMENSION OF BIOGENIC AEROSOL ARRAY
         , nitrate_dim1                                                 &
!             1ST DIMENSION OF NITRATE AEROSOL ARRAY
         , nitrate_dim2                                                 &
!             2ND DIMENSION OF NITRATE AEROSOL ARRAY
         , n_cca_lev
!             NUMBER OF LEVELS FOR CONVECTIVE CLOUD AMOUNT

!     ACTUAL SIZES USED:
INTEGER ::                                                        &
          !, INTENT(IN)
     n_profile                                                    &
!             NUMBER OF PROFILES
         , n_layer                                                      &
!             Number of layers seen by the radiation code
         , nclds
!             NUMBER OF CLOUDY LEVELS

!     GATHERING ARRAY:
INTEGER ::                                                        &
          !, INTENT(IN)
     i_gather(npd_field)
!             LIST OF POINTS TO BE GATHERED
INTEGER, INTENT(IN) :: col_list(npd_field)
!                              EW indices of gathered points in the 2-D
!                              domain
INTEGER, INTENT(IN) :: row_list(npd_field)
!                              NS indices of gathered points in the 2-D
!                              domain

!     THERMODYNAMIC FIELDS:
REAL ::                                                           &
          !, INTENT(IN)
     p(npd_profile, npd_layer)                                    &
!             PRESSURES
         , t(npd_profile, npd_layer)                              &
!             TEMPERATURES
         , d_mass(npd_profile, npd_layer)                         &
!             MASS THICKNESSES OF LAYERS
         , alat(npd_profile)
!             Latitude in degree

!     CONVECTIVE CLOUDS:
INTEGER ::                                                        &
          !, INTENT(IN)
     ccb_in(npd_field)                                            &
!             BASE OF CONVECTIVE CLOUD
         , cct_in(npd_field)
!             TOP OF CONVECTIVE CLOUD

INTEGER, INTENT(IN):: lcbase(npd_field)
!             Lowest cloud base level in vertical profile.


REAL ::                                                           &
          !, INTENT(IN)
     cca(npd_field,n_cca_lev)                                     &
!             FRACTION OF CONVECTIVE CLOUD
         , cccwp(npd_field)
!             WATER PATH OF CONVECTIVE CLOUD

REAL, INTENT(IN) :: ccw(npd_field, model_levels)
!             Convective Cloud Water (kg/kg) (If CCRad =.true.)

LOGICAL ::                                                        &
          !, INTENT(IN)
     l_local_cnv_partition                                        &
!             FLAG TO CARRY OUT THE PARTITIONING BETWEEN ICE
!             AND WATER IN CONVECTIVE CLOUDS AS A FUNCTION OF
!             THE LOCAL TEMPERATURE
         , l_seasalt_ccn                                                &
!              FLAG FOR SEA-SALT PARAMETRIZATION FOR CCN
         , l_use_biogenic                                               &
!              FLAG TO USE BIOGENIC AEROSOLS AS CCN
         , l_biomass_ccn                                                &
!              FLAG FOR BIOMASS PARAMETRIZATION FOR CCN
         , l_ocff_ccn                                                   &
!              FLAG FOR FOSSIL-FUEL ORG CARB PARAMETRIZATION FOR CCN
         , l_nitrate_ccn                                                &
!              FLAG FOR NITRATE AEROSOL AS CCN
         , l_easyaerosol_cdnc
!              Flag to use EasyAerosol CDNC

!     LAYER CLOUDS:
REAL ::                                                           &
          !, INTENT(IN)
     lccwc1(npd_field, nclds+1/(nclds+1))                         &
!             LIQUID WATER CONTENTS
         , lccwc2(npd_field, nclds+1/(nclds+1))                         &
!             ICE WATER CONTENTS
         , lca_area(npd_field, nclds+1/(nclds+1))                       &
!             AREA COVERAGE FRACTIONS OF LAYER CLOUDS
         , lca_bulk(npd_field, nclds+1/(nclds+1))
!             BULK COVERAGE FRACTIONS OF LAYER CLOUDS

REAL, INTENT(IN) :: n_drop_pot(npd_field, nclds)

!     ARRAYS FOR MICROPHYSICS:
LOGICAL ::                                                        &
          !, INTENT(IN)
     l_microphysics                                               &
!             MICROPHYSICAL FLAG
         , l_aerosol_ccn                                                &
!             FLAG TO USE AEROSOLS TO FIND CCN
         , land_g(npd_profile)
!             FLAG FOR LAND POINTS
REAL ::                                                           &
          !, INTENT(IN)
     accum_sulphate(sulp_dim1, sulp_dim2)                         &
!             MIXING RATIOS OF ACCUMULATION-MODE SULPHATE
         , aitken_sulphate(sulp_dim1, sulp_dim2)                        &
!             Mixing ratios of Aitken-mode sulphate
         , diss_sulphate(sulp_dim1, sulp_dim2)                          &
!             MIXING RATIOS OF DISSOLVED SULPHATE
         , aged_bmass(bmass_dim1, bmass_dim2)                           &
!             MIXING RATIOS OF AGED BIOMASS SMOKE
         , cloud_bmass(bmass_dim1, bmass_dim2)                          &
!             MIXING RATIOS OF IN-CLOUD BIOMASS SMOKE
         , aged_ocff(ocff_dim1, ocff_dim2)                              &
!             MIXING RATIOS OF AGED FOSSIL-FUEL ORGANIC CARBON
         , cloud_ocff(ocff_dim1, ocff_dim2)                             &
!             MIXING RATIOS OF IN-CLOUD FOSSIL-FUEL ORGANIC CARBON
         , sea_salt_film(salt_dim1, salt_dim2, salt_dim3)               &
!             NUMBER CONCENTRATION OF FILM-MODE SEA-SALT AEROSOL
         , sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)                &
!             NUMBER CONCENTRATION OF JET-MODE SEA-SALT AEROSOL
         , biogenic(biogenic_dim1, biogenic_dim2)                       &
!             M.M.R. OF BIOGENIC AEROSOL
         , accum_nitrate(nitrate_dim1, nitrate_dim2)                    &
!             MIXING RATIOS OF ACCUMULATION NITRATE AEROSOL
         , diss_nitrate(nitrate_dim1, nitrate_dim2)
!             MIXING RATIOS OF DISSOLVED NITRATE AEROSOL

! EasyAerosol cloud droplet number concentrations
TYPE (t_easyaerosol_cdnc), INTENT(IN) :: easyaerosol_cdnc

!     REPRESENTATION OF CLOUDS
INTEGER ::                                                        &
          !, INTENT(IN)
     i_cloud_representation
!             REPRESENTATION OF CLOUDS

!     PARAMETRIZATIONS FOR CLOUDS:
INTEGER ::                                                        &
          !, INTENT(IN)
     i_condensed_param(npd_cloud_component)
!             TYPES OF PARAMETRIZATION USED FOR CONDENSED
!             COMPONENTS IN CLOUDS
!     LIMITS ON SIZES OF PARTICLES
REAL ::                                                           &
          !, INTENT(IN)
     condensed_min_dim(npd_cloud_component)                       &
!             MINIMUM DIMENSION OF EACH CONDENSED COMPONENT
         , condensed_max_dim(npd_cloud_component)
!             MAXIMUM DIMENSION OF EACH CONDENSED COMPONENT

!     ASSIGNED CLOUD FIELDS:
INTEGER ::                                                        &
          !, INTENT(OUT)
     n_condensed                                                  &
!             NUMBER OF CONDENSED COMPONENTS
         , type_condensed(npd_cloud_component)                          &
!             TYPES OF CONDENSED COMPONENTS
         , n_cloud_type
!             Number of types of clouds
REAL ::                                                           &
          !, INTENT(OUT)
     w_cloud(npd_profile, id_ct:npd_layer)                        &
!             TOTAL AMOUNTS OF CLOUD
         , frac_cloud(npd_profile, id_ct:npd_layer, npd_cloud_type)     &
!             FRACTION OF EACH TYPE OF CLOUD
         , condensed_dim_char(npd_profile, id_ct:npd_layer              &
            , npd_cloud_component)                                      &
!             CHARACTERISTIC DIMENSIONS OF CLOUDY COMPONENTS
         , CDNC(npd_profile, id_ct:npd_layer, npd_cloud_component)      &
!             CDNC ARRAY FOR GETTING CDNC AT CLOUD TOP
         , condensed_mix_rat_area(npd_profile, id_ct:npd_layer          &
            , npd_cloud_component)                                      &
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING AREA CLD
         , ntot_diag_g(npd_profile, npd_layer)                          &
!             DIAGNOSTIC ARRAY FOR NTOT (GATHERED)
         , strat_lwc_diag_g(npd_profile, npd_layer)                     &
!             DIAGNOSTIC ARRAY FOR STRATIFORM LWC (GATHERED)
         , so4_ccn_diag_g(npd_profile, npd_layer)
!             DIAGNOSTIC ARRAY FOR SO4 CCN MASS CONC (GATHERED)


REAL ::                                                           &
     lying_snow(npd_field)                                        &
!            SNOW DEPTH (>5000m = LAND ICE SHEET)
         , lying_snow_g(npd_profile)                                    &
!            GATHERED VERSION OF THE ABOVE
         , flandg_g(npd_profile)
!            GATHERED GLOBAL LAND FRACTION FIELD

! Diagnostics:
TYPE (StrDiag), INTENT(INOUT) :: diag



!     LOCAL VARIABLES:
INTEGER ::                                                        &
     i                                                            &
!             LOOP VARIABLE
         , j                                                            &
!             LOOP VARIABLE
         , k                                                            &
!             LOOP VARIABLE
         , l                                                            &
!             LOOP VARIABLE
         , lg
!             INDEX TO GATHER

LOGICAL ::                                                        &
     l_glaciated_top(npd_profile)
!             LOGICAL FOR GLACIATED TOPS IN CONVECTIVE CLOUD.

REAL ::                                                           &
     liq_frac(npd_profile)                                        &
!             FRACTION OF LIQUID CLOUD WATER
         , liq_frac_conv(npd_profile)                                   &
!             FRACTION OF LIQUID WATER IN CONVECTIVE CLOUD
         , t_gather(npd_profile)                                        &
!             GATHERED TEMPERATURE FOR LSP_FOCWWIL
         , t_list(npd_profile)                                          &
!             LIST OF TEMPERATURES
         , total_mass(npd_profile)                                      &
!             TOTAL MASS IN CONVECTIVE CLOUD
         , cc_depth(npd_profile)                                        &
!             DEPTH OF CONVECTIVE CLOUD
         , condensed_mix_rat_bulk(npd_profile, id_ct:npd_layer          &
            , npd_cloud_component)                                      &
!             MASS MIXING RATIOS OF CONDENSED COMPONENTS USING BULK CLD
         , density_air(npd_profile, npd_layer)                          &
!             DENSITY OF AIR
         , convective_cloud_layer(npd_profile)                          &
!             AMOUNT OF CONVECTIVE CLOUD IN TH CURRENT LAYER
         , strat_liq_cloud_fraction(npd_profile, npd_layer)             &
!             STRATIFORM LIQUID CLOUD FRACTION (T>273K)
         , conv_liq_cloud_fraction(npd_profile, npd_layer)              &
!             CONVECTIVE LIQUID CLOUD FRACTION (T>273K)
         , nc_diag_g(npd_profile)                                       &
!             DIAGNOSTIC ARRAY FOR COLUMN DROPLET NUMBER (GATHERED)
         , nc_weight_g(npd_profile)
!             DIAGNOSTIC ARRAY FOR COL. DROP NO. SAMPLING WGT (GATHERED)

REAL :: work_qcl(n_profile, n_layer), work_qcf(n_profile, n_layer)

      !-----------------------------------------------------------------
      ! Convection variables - CCRad
      !-----------------------------------------------------------------
INTEGER :: ccb(npd_field)
!             CONVECTIVE CLOUD BASE TO BE USED BY RADIATION.

INTEGER :: cct(npd_field)
!             CONVECTIVE CLOUD TOP TO BE USED BY RADIATION.

REAL    :: lca_of_grdbx (npd_field, MAX(nclds,1))
!             LARGE-SCALE CLOUD AREA OF GRIDBOX.

REAL    :: iwc
!             Ice water content


!     Parameters for the aggregate parametrization.
REAL, PARAMETER :: a0_agg_cold = 7.5094588e-04
REAL, PARAMETER :: b0_agg_cold = 5.0830326e-07
REAL, PARAMETER :: a0_agg_warm = 1.3505403e-04
REAL, PARAMETER :: b0_agg_warm = 2.6517429e-05
REAL, PARAMETER :: t_switch    = 216.208
REAL, PARAMETER :: t0_agg      = 279.5
REAL, PARAMETER :: s0_agg      = 0.05

!     Coefficients for effective size of ice crystal from Sun scheme
REAL, PARAMETER :: a1_s=45.9866
REAL, PARAMETER :: a2_s=0.2214
REAL, PARAMETER :: a3_s=0.7957
REAL, PARAMETER :: a4_s=0.2535
REAL, PARAMETER :: a5_s=83.15
REAL, PARAMETER :: a6_s=1.2351
REAL, PARAMETER :: a7_s=0.0105
!     Latitude threshold for division of ice crystal correction
REAL, PARAMETER :: a8_s=30.0

!     Parameters for Baran's diagnostic parameterisation of effective size.
!     m and n are the slope and intercept (positive sign), and max and min the 
!     maximum and minimum values. Values for parameterisation expressed in
!     units of metres.
REAL, PARAMETER :: m_ice_baran = 1.868e-6
REAL, PARAMETER :: n_ice_baran = 353.613e-6
REAL, PARAMETER :: min_ice_baran = 7.0e-6
REAL, PARAMETER :: max_ice_baran = 156.631e-6

!     Tolerance for cloud fraction to prevent unrealistic large LWC
REAL :: lca_tol

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='R2_SET_CLOUD_FIELD'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

lca_tol = EPSILON(lca_tol) ! Was 1.e-6 for SES2

DO l=1, n_profile
  l_glaciated_top(l) = .FALSE.
  t_gather(l)        =  0.0e+00
END DO

!     SET THE TOTAL AMOUNTS OF CLOUD AND THE FRACTIONS COMPRISED BY
!     CONVECTIVE AND STRATIFORM COMPONENTS.

!     ZERO THE AMOUNTS OF CLOUD IN THE UPPER LAYERS.
DO i=id_ct, n_layer-nclds
  DO l=1, n_profile
    w_cloud(l, i)=0.0e+00
  END DO
END DO

! Initialise the characterisitc dimension
condensed_dim_char = 0.0

IF (lrad_ccrad) THEN
  !---------------------------------------------------------------------
  ! Convection passes in ccb_in and cct_in on model layers, and these
  ! only to refer the top most contigiuos cloud bank.  Set local
  ! convective cloud base and top to:
  !
  !   ccb = LOWER BOUNDARY of LOWEST  LAYER with convective cloud
  !   cct = UPPER BOUNDARY of HIGHEST LAYER with convective cloud
  !
  ! as expected by the radiation scheme.
  !---------------------------------------------------------------------

  DO l=1, n_profile
    lg=i_gather(l)

    ccb(lg) = lcbase(lg)     ! Set ccb to lowest cloud boundary
                             ! of lowest cloud
    cct(lg) = cct_in(lg) + 1 ! Set cct to upper boundary of highest
                             ! Cloud
  END DO

ELSE ! Original

  DO l=1, n_profile
    lg=i_gather(l)
    ccb(lg) = ccb_in(lg)
    cct(lg) = cct_in(lg)
  END DO

END IF      ! lrad_ccrad



IF (i_cloud_representation == ip_cloud_conv_strat) THEN

  n_cloud_type=2

  !  This cloud representation not available with new cloud microphysics

  !        THE CLOUDS ARE DIVIDED INTO MIXED-PHASE STRATIFORM AND
  !        CONVECTIVE CLOUDS: LSP_FOCWWIL GIVES THE PARTITIONING BETWEEN
  !        ICE AND WATER IN STRATIFORM CLOUDS AND IN CONVECTIVE CLOUD,
  !        UNLESS THE OPTION TO PARTITION AS A FUNCTION OF THE LOCAL
  !        TEMPERATURE IS SELECTED. WITHIN CONVECTIVE CLOUD THE LIQUID
  !        WATER CONTENT IS DISTRIBUTED UNIFORMLY THROUGHOUT THE CLOUD.


  !        Set the components within the clouds. Here we have four
  !        components: stratiform ice and water and convective
  !        ice and water.
  n_condensed=4
  type_condensed(1)=ip_clcmp_st_water
  type_condensed(2)=ip_clcmp_st_ice
  type_condensed(3)=ip_clcmp_cnv_water
  type_condensed(4)=ip_clcmp_cnv_ice


  !        CONVECTIVE CLOUD:

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)=0.0e+00
      condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)=0.0e+00
      condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)=0.0e+00
      condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)=0.0e+00
    END DO
  END DO


  IF (l_local_cnv_partition) THEN

    !           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
    !           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
    !           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT
    !           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.

    DO l=1, n_profile
      IF ((n_layer+2-cct(i_gather(l))) <= n_layer) THEN
        l_glaciated_top(l) =                                    &
          (t(l, MAX( n_layer+2-cct(i_gather(l))                 &
                   , n_layer-nclds+1)) <  tm)
      END IF
    END DO

  ELSE

    !           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
    !           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
    !           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.

    DO l=1, n_profile
      IF ((n_layer+2-cct(i_gather(l))) <= n_layer) THEN

        t_gather(l) = t(l, MAX( n_layer+2-cct(i_gather(l))      &
                             , n_layer-nclds+1))
      END IF
    END DO

    CALL lsp_focwwil(t_gather, n_profile, liq_frac_conv)

  END IF ! Test on local_cnv_partition


  DO l=1, n_profile
    total_mass(l)=0.0e+00
  END DO

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      lg=i_gather(l)
      IF ( (cct(lg) >= n_layer+2-i) .AND.                       &
           (ccb(lg) <= n_layer+1-i) ) THEN
        total_mass(l)=total_mass(l)+d_mass(l, i)
      END IF
    END DO
  END DO

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      lg=i_gather(l)
      IF ( (cct(lg) >= n_layer+2-i) .AND.                       &
           (ccb(lg) <= n_layer+1-i) ) THEN
        IF (l_local_cnv_partition) THEN
          !                    THE PARTITIONING IS RECALCULATED FOR EACH LAYER
          !                    OTHERWISE A GENERIC VALUE IS USED.
          liq_frac_conv(l)=MAX(0.0e+00, MIN(1.0e+00          &
             , 1.61e-02*(t(l, i)-tm)+8.9e-01))
          !                    Do not allow ice above 0 Celsius unless the top
          !                    of the cloud is glaciated and force homogeneous
          !                    nucleation at -40 Celsius.
          IF ( (t(l, i) >  tm) .AND.                          &
               (.NOT. l_glaciated_top(l)) ) THEN
            liq_frac_conv(l)=1.0e+00
          ELSE IF (t(l, i) <  tm-4.0e+01) THEN
            liq_frac_conv(l)=0.0e+00
          END IF
        END IF
        condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)      &
           =cccwp(lg)*liq_frac_conv(l)                        &
           /(total_mass(l)+TINY(cccwp))
        condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)        &
           =cccwp(lg)*(1.0e+00-liq_frac_conv(l))              &
           /(total_mass(l)+TINY(cccwp))
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)      &
           =cccwp(lg)*liq_frac_conv(l)                        &
           /(total_mass(l)+TINY(cccwp))
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)        &
           =cccwp(lg)*(1.0e+00-liq_frac_conv(l))              &
           /(total_mass(l)+TINY(cccwp))
      END IF
    END DO
  END DO


  !        STRATIFORM CLOUDS:

  !        PARTITION BETWEEN ICE AND WATER DEPENDING ON THE
  !        LOCAL TEMPERATURE.

  DO i=1, nclds
    DO l=1, n_profile
      lg=i_gather(l)
      IF (lca_area(lg, i) > lca_tol) THEN
        liq_frac(l) = lccwc1(lg, i)                            &
          / (lccwc1(lg, i)+lccwc2(lg, i))
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_water)                                &
            =(lccwc1(lg, i)+lccwc2(lg, i))                     &
            *liq_frac(l)/lca_area(lg, i)
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)                                  &
            =(lccwc1(lg, i)+lccwc2(lg, i))                     &
            *(1.0e+00-liq_frac(l))/lca_area(lg, i)
      ELSE
        liq_frac(l) = 0.0
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_water)=0.0e+00
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)=0.0e+00
      END IF

      IF (lca_bulk(lg, i) > lca_tol) THEN
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_water)                                &
            =(lccwc1(lg, i)+lccwc2(lg, i))                     &
            *liq_frac(l)/lca_bulk(lg, i)
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)                                  &
            =(lccwc1(lg, i)+lccwc2(lg, i))                     &
            *(1.0e+00-liq_frac(l))/lca_bulk(lg, i)
      ELSE
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_water)=0.0e+00
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)=0.0e+00
      END IF
    END DO
  END DO


  !        CLOUD FRACTIONS:

  IF (l_3d_cca) THEN
    DO i=1, nclds
      DO l=1, n_profile
        lg=i_gather(l)
        w_cloud(l, n_layer+1-i)                                  &
           =cca(lg,i)+(1.0e+00-cca(lg,i))*lca_area(lg, i)
        frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)           &
           =cca(lg,i)/(w_cloud(l, n_layer+1-i)+TINY(cca))
        frac_cloud(l, n_layer+1-i, ip_cloud_type_strat)          &
           =1.0e+00-frac_cloud(l, n_layer+1-i                    &
           , ip_cloud_type_conv)
      END DO
    END DO
  ELSE
    DO i=1, nclds
      DO l=1, n_profile
        lg=i_gather(l)
        IF ( (i <= cct(lg)-1) .AND. (i >= ccb(lg)) ) THEN
          w_cloud(l, n_layer+1-i)                               &
             =cca(lg,1)+(1.0e+00-cca(lg,1))*lca_area(lg, i)
          frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)        &
             =cca(lg,1)/(w_cloud(l, n_layer+1-i)+TINY(cca))
        ELSE
          w_cloud(l, n_layer+1-i)=lca_area(lg, i)
          frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)        &
             =0.0e+00
        END IF
        frac_cloud(l, n_layer+1-i, ip_cloud_type_strat)          &
           =1.0e+00-frac_cloud(l, n_layer+1-i                    &
           , ip_cloud_type_conv)
      END DO
    END DO
  END IF




ELSE IF (i_cloud_representation == ip_cloud_csiw) THEN

  n_cloud_type=4

  !        HERE THE CLOUDS ARE SPLIT INTO FOUR SEPARATE TYPES.
  !        THE PARTITIONING BETWEEN ICE AND WATER IS REGARDED AS
  !        DETERMINING THE AREAS WITHIN THE GRID_BOX COVERED BY
  !        ICE OR WATER CLOUD, RATHER THAN AS DETERMINING THE IN-CLOUD
  !        MIXING RATIOS. THE GRID-BOX MEAN ICE WATER CONTENTS IN
  !        STRATIFORM CLOUDS MAY BE PREDICTED BY THE ICE MICROPHYSICS
  !        SCHEME OR MAY BE DETERMINED AS A FUNCTION OF THE TEMPERATURE
  !        (LSP_FOCWWIL). IN CONVECTIVE CLOUDS THE PARTITIONING MAY BE
  !        DONE USING THE SAME FUNCTION, LSP_FOCWWIL, BASED ON A SINGLE
  !        TEMPERATURE, OR USING A PARTITION BASED ON THE LOCAL
  !        TEMPERATURE.


  !        Set the components within the clouds. Here we have four
  !        components: stratiform ice and water and convective
  !        ice and water.
  n_condensed=4
  type_condensed(1)=ip_clcmp_st_water
  type_condensed(2)=ip_clcmp_st_ice
  type_condensed(3)=ip_clcmp_cnv_water
  type_condensed(4)=ip_clcmp_cnv_ice


  !        CONVECTIVE CLOUD:


  IF (lrad_ccrad) THEN

    ! Initialise Convection mixing ratio arrays on ALL levels, not just
    ! the levels that radiation is using (may affect diagnostics
    ! otherwise)

    DO i=id_ct, npd_layer
      DO l=1, npd_profile
        condensed_mix_rat_area(l, i, ip_clcmp_cnv_water) = 0.0e+00
        condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)   = 0.0e+00
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water) = 0.0e+00
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)   = 0.0e+00
      END DO      ! L
    END DO      ! I

  ELSE      ! Original Code

    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)  = 0.0e+00
        condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)    = 0.0e+00
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)  = 0.0e+00
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)    = 0.0e+00
      END DO
    END DO
  END IF      ! lrad_ccrad



  IF (lrad_ccrad) THEN
    !-------------------------------------------------------------------
    ! USE CCW PROFILE WHICH WAS PASSED IN FROM CONVECTION
    !-------------------------------------------------------------------

    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        lg = i_gather(l)

        condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)                &
                                                   = ccw(lg,n_layer+1-i)
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)                &
                                                   = ccw(lg,n_layer+1-i)

        condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)                  &
                      = condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)                  &
                      = condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)

      END DO      ! L (N_LAYER)
    END DO      ! I (N_PROFILE)

  ELSE ! Original Code, use cclwp spread from cloud bottom to top

    DO l=1, n_profile
      total_mass(l) = 0.0e+00
    END DO


    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        lg = i_gather(l)
        IF ( (cct(lg) >= n_layer+2-i) .AND.                              &
             (ccb(lg) <= n_layer+1-i) ) THEN
          total_mass(l) = total_mass(l) + d_mass(l, i)
        END IF
      END DO      ! L
    END DO      ! I


    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        lg = i_gather(l)

        IF ( (cct(lg) >= n_layer+2-i) .AND.                              &
             (ccb(lg) <= n_layer+1-i) ) THEN

          condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)              &
                                 = cccwp(lg)/(total_mass(l)+TINY(cccwp))
          condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)              &
                                 = cccwp(lg)/(total_mass(l)+TINY(cccwp))

          condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)                &
                       =condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)
          condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)                &
                       =condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)

        END IF

      END DO      ! L
    END DO      ! I
  END IF      ! lrad_ccrad

  !        STRATIFORM CLOUDS:

  DO i=1, nclds
    DO l=1, n_profile
      lg=i_gather(l)
      IF (lca_area(lg, i) > lca_tol) THEN
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_water)                                &
           =(lccwc1(lg, i)+lccwc2(lg, i))/lca_area(lg, i)
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)                                  &
           =condensed_mix_rat_area(l, n_layer+1-i              &
           , ip_clcmp_st_water)
      ELSE
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_water)=0.0e+00
        condensed_mix_rat_area(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)=0.0e+00
      END IF

      IF (lca_bulk(lg, i) > lca_tol) THEN
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_water)                                &
           =(lccwc1(lg, i)+lccwc2(lg, i))/lca_bulk(lg, i)
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)                                  &
           =condensed_mix_rat_bulk(l, n_layer+1-i              &
           , ip_clcmp_st_water)
      ELSE
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_water)=0.0e+00
        condensed_mix_rat_bulk(l, n_layer+1-i                  &
           , ip_clcmp_st_ice)=0.0e+00
      END IF
    END DO
  END DO



  !        CLOUD FRACTIONS:

  IF (l_local_cnv_partition) THEN

    !           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
    !           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
    !           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT
    !           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.

    DO l=1, n_profile
      IF ((n_layer+2-cct(i_gather(l))) <= n_layer) THEN
        l_glaciated_top(l) =                                    &
          (t(l, MAX(n_layer+2-cct(i_gather(l))                  &
                  , n_layer-nclds+1)) <  tm)
      END IF
    END DO

  ELSE

    !           PARTITION BETWEEN ICE AND WATER AS DIRECTED BY THE
    !           TEMPERATURE IN THE MIDDLE OF THE TOP LAYER OF THE CLOUD.
    !           THE PARTITIONING MAY BE PRECALCULATED IN THIS CASE.

    DO l=1, n_profile
      IF ((n_layer+2-cct(i_gather(l))) <= n_layer) THEN
        t_gather(l) = t(l, MAX(n_layer+2-cct(i_gather(l))       &
                             , n_layer-nclds+1))
      END IF
    END DO

    CALL lsp_focwwil(t_gather, n_profile, liq_frac_conv)

  END IF



  IF (lrad_ccrad .AND. l_rad_ovrlap) THEN

    DO i=n_layer+1-nclds, n_layer

      !----------------------------------------------------------------
      ! Calculate cloud fractions of respective cloud types
      !----------------------------------------------------------------
      DO l=1, n_profile

        ! Modify Convective/Large-scale cloud fractions to account. For
        ! overlap, CCA overlaps LCA and takes precedence, i.e. LCA must
        ! be > CCA in order to have a non-zero cloud fraction.

        lg                            = i_gather(l)
        w_cloud(l, i)                 = 0.0
        lca_of_grdbx(lg, n_layer+1-i) = 0.0
        convective_cloud_layer(l)     = cca(lg,n_layer+1-i)

        ! Calculate Total Cloud Cover in gridbox
        w_cloud(l, i) = MAX(  convective_cloud_layer(l)                &
                            , lca_area(lg, n_layer+1-i))


        !--------------------------------------------------------------
        ! Scavenge LS Cloud water by CCW
        !--------------------------------------------------------------
        IF (l_rad_ccw_scav) THEN

          ! If Convective cloud and Large-scale cloud are non-zero and
          ! occur on the same model level. Recalc. CCW to allow for the
          ! additional condensate from the large-scale cloud

          IF ((lca_area(lg, n_layer+1-i) > 1.0e-10) .AND.              &
              (convective_cloud_layer(l) > 1.0e-10)) THEN

            condensed_mix_rat_area(l,i,ip_clcmp_cnv_water)             &
              = ccw(lg,n_layer+1-i)                                    &
              + (  (lccwc1(lg, i)+lccwc2(lg, i))                       &
                 * MIN(  lca_area(lg, n_layer+1-i)                     &
                       , convective_cloud_layer(l))                    &
                 / convective_cloud_layer(l) )

            condensed_mix_rat_area(l,i,ip_clcmp_cnv_ice)               &
              = condensed_mix_rat_area(l,i,ip_clcmp_cnv_water)

          END IF
        END IF      ! l_rad_ccw_scav


        !--------------------------------------------------------------
        ! Recalc LCA_OF_GRDBX if there is more large-scale cloud
        ! fraction than cca
        !--------------------------------------------------------------
        IF (  (lca_area(lg, n_layer+1-i) - convective_cloud_layer(l))  &
            > 1.0e-10) THEN

          lca_of_grdbx(lg, n_layer+1-i) = lca_area(lg, n_layer+1-i)    &
                                        - convective_cloud_layer(l)
        END IF


        !--------------------------------------------------------------
        ! Partition stratiform clouds using the ratio of cloud water
        ! contents.
        !--------------------------------------------------------------
        IF (( lca_of_grdbx(lg, n_layer+1-i) > EPSILON(lca_of_grdbx))   &
           .AND.                                                       &
           (( lccwc1(lg,n_layer+1-i) + lccwc2(lg,n_layer+1-i))         &
             > 0.0)) THEN

          liq_frac(l) = lccwc1(lg, n_layer+1-i)                        &
                      / (  lccwc1(lg, n_layer+1-i)                     &
                         + lccwc2(lg, n_layer+1-i))
        ELSE
          liq_frac(l) = 0.0e+00
        END IF


        !--------------------------------------------------------------
        ! Partition Convective liquid/ice cloud fraction based on local
        ! temperature.
        !--------------------------------------------------------------
        IF (l_local_cnv_partition) THEN

          ! NOTE: Ice is allowed above the freezing point ONLY if the
          !       TOP of the cloud is glaciated.

          liq_frac_conv(l) = MAX(  0.0e+00                             &
                                 , MIN(  1.0e+00                       &
                                       , 1.61e-02*(t(l, i)-tm)         &
                                         + 8.9e-01))

          ! Do not allow ice above 0 Celsius unless the top of the
          ! cloud is glaciated and force homogeneous nucleation at -40
          ! Celsius

          IF ( (t(l, i) > tm) .AND.                                    &
               (.NOT. l_glaciated_top(l)) ) THEN
            liq_frac_conv(l) = 1.0e+00
          ELSE IF (t(l, i) < tm-4.0e+01) THEN
            liq_frac_conv(l) = 0.0e+00
          END IF

        END IF      ! L_LOCAL_CNV_PARTITION


        !--------------------------------------------------------------
        ! Split cloudly fraction of gridbox into ls/convective
        ! liquid/ice cloud fractions. i.e. Fraction of (fraction of
        ! gridbox).
        !--------------------------------------------------------------

        frac_cloud(l,i,ip_cloud_type_sw)                               &
                         = (liq_frac(l) * lca_of_grdbx(lg,n_layer+1-i))&
                         / (w_cloud(l,i) + TINY(w_cloud))

        frac_cloud(l,i,ip_cloud_type_si)                               &
                         = ((1.0e+00-liq_frac(l))                      &
                         * lca_of_grdbx(lg, n_layer+1-i))              &
                         / (w_cloud(l, i) + TINY(w_cloud))

        IF (.NOT. l_local_cnv_partition         &
            .AND. (t(l,i) >= tm)) THEN
            !--------------------------------------------------------
            ! Convective cloud gridbox fraction (liquid)
            !--------------------------------------------------------
          frac_cloud(l,i,ip_cloud_type_cw)                         &
                      = convective_cloud_layer(l)                  &
                      / (w_cloud(l,i) + TINY(w_cloud))

          !--------------------------------------------------------
          ! Convective cloud gridbox fraction (ice)
          !--------------------------------------------------------
          frac_cloud(l,i,ip_cloud_type_ci) = 0.0

        ELSE
          frac_cloud(l, i, ip_cloud_type_cw)                             &
                           = liq_frac_conv(l) * convective_cloud_layer(l)&
                           / (w_cloud(l, i) + TINY(w_cloud))

          frac_cloud(l, i, ip_cloud_type_ci)                             &
                           = (1.0e+00-liq_frac_conv(l))                  &
                           * convective_cloud_layer(l)                   &
                           / (w_cloud(l, i) + TINY(w_cloud))
        END IF

      END DO      ! L (N_PROFILE)
    END DO      ! I (N_LAYER)

  ELSE     ! Original Code

    DO i=n_layer+1-nclds, n_layer

      IF (l_3d_cca) THEN
        DO l=1, n_profile
          lg = i_gather(l)
          convective_cloud_layer(l)=cca(lg,n_layer+1-i)
        END DO
      ELSE
        DO l=1, n_profile
          lg = i_gather(l)
          IF ( (cct(lg) >= n_layer+2-i) .AND.                           &
               (ccb(lg) <= n_layer+1-i) ) THEN
            convective_cloud_layer(l)=cca(lg,1)
          ELSE
            convective_cloud_layer(l)=0.0e+00
          END IF
        END DO
      END IF


      DO l=1, n_profile
        lg = i_gather(l)

        w_cloud(l, i)                                                  &
          = convective_cloud_layer(l)                                  &
          + (1.0e+00 - convective_cloud_layer(l))                      &
          * lca_area(lg, n_layer+1-i)



        ! Partition stratiform clouds using the ratio of cloud water
        ! contents.

        IF (lca_area(lg, n_layer+1-i) > EPSILON(lca_area)              &
           .AND.                                                       &
           (lccwc1(lg, n_layer+1-i) + lccwc2(lg, n_layer+1-i))         &
            > 0.0) THEN

          liq_frac(l) =  lccwc1(lg, n_layer+1-i)                       &
                      / (lccwc1(lg, n_layer+1-i)                       &
                      +  lccwc2(lg, n_layer+1-i))
        ELSE
          liq_frac(l) = 0.0e+00
        END IF



        IF (l_local_cnv_partition) THEN

          ! The partitioning between ice and water must be recalculated
          ! for this layer as a function of the local temperature, but
          ! ice is allowed above the freezing point only if the top of
          ! the cloud is glaciated.

          liq_frac_conv(l)=MAX(  0.0e+00                               &
                               , MIN(  1.0e+00                         &
                                     , 1.61e-02*(t(l, i)-tm)+8.9e-01))

          ! Do not allow ice above 0 Celsius unless the top of the
          ! cloud is glaciated and force homogeneous nucleation at -40
          ! Celsius.

          IF ( (t(l, i) > tm) .AND. (.NOT. l_glaciated_top(l)) ) THEN
            liq_frac_conv(l) = 1.0e+00
          ELSE IF (t(l, i) <  tm-4.0e+01) THEN
            liq_frac_conv(l) = 0.0e+00
          END IF

        END IF

        frac_cloud(l, i, ip_cloud_type_sw)                             &
                  = liq_frac(l) * (1.0e+00-convective_cloud_layer(l))  &
                  * lca_area(lg, n_layer+1-i)                          &
                  / (w_cloud(l, i) + TINY(w_cloud))

        frac_cloud(l, i, ip_cloud_type_si)                             &
                  = (1.0e+00 - liq_frac(l))                            &
                  * (1.0e+00 - convective_cloud_layer(l))              &
                  * lca_area(lg, n_layer+1-i)                          &
                  / (w_cloud(l, i) + TINY(w_cloud))



        IF (.NOT. l_local_cnv_partition                                &
            .AND. (t(l,i) >= tm)) THEN
            !--------------------------------------------------------
            ! Convective cloud gridbox fraction (liquid)
            !--------------------------------------------------------
          frac_cloud(l,i,ip_cloud_type_cw)                           &
                      = convective_cloud_layer(l)                    &
                      / (w_cloud(l,i) + TINY(w_cloud))

          !--------------------------------------------------------
          ! Convective cloud gridbox fraction (ice)
          !--------------------------------------------------------
          frac_cloud(l,i,ip_cloud_type_ci) = 0.0

        ELSE


          frac_cloud(l, i, ip_cloud_type_cw)                             &
                    = liq_frac_conv(l) * convective_cloud_layer(l)       &
                    / (w_cloud(l, i) + TINY(w_cloud))

          frac_cloud(l, i, ip_cloud_type_ci)                             &
                    = (1.0e+00 - liq_frac_conv(l))                       &
                    * convective_cloud_layer(l)                          &
                    / (w_cloud(l, i) + TINY(w_cloud))
        END IF
      END DO      ! L (N_PROFILE)
    END DO      ! I (N_LAYER)
  END IF      ! lrad_ccrad .AND. l_overlap



  !--------------------------------------------------------------------
  ! Cloud fraction/Mixing Ratio consistency check
  !--------------------------------------------------------------------
  IF (lrad_ccrad) THEN

    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        lg = i_gather(l)
        IF (frac_cloud(l, i, ip_cloud_type_sw) <= 1.0e-10)             &
                 condensed_mix_rat_area(l, i, ip_clcmp_st_water)  = 0.0
        IF (frac_cloud(l, i, ip_cloud_type_si) <= 1.0e-10)             &
                 condensed_mix_rat_area(l, i, ip_clcmp_st_ice)    = 0.0
        IF (frac_cloud(l, i, ip_cloud_type_cw) <= 1.0e-10)             &
                 condensed_mix_rat_area(l, i, ip_clcmp_cnv_water) = 0.0
        IF (frac_cloud(l, i, ip_cloud_type_ci) <= 1.0e-10)             &
                 condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)   = 0.0
      END DO      ! L (N_PROFILE)
    END DO      ! I (N_LAYER)


    ! For consisitency with original code, bulk ratios
    ! are set to be the same as the area mixing ratios
    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        lg = i_gather(l)
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water)               &
               = condensed_mix_rat_area(l, i, ip_clcmp_cnv_water)
        condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)                 &
               = condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)
      END DO      ! L (N_PROFILE)
    END DO      ! I (N_LAYER)

  END IF      ! lrad_ccrad


ELSE IF (i_cloud_representation == ip_cloud_ice_water) THEN

  n_cloud_type=2

  ! Here the clouds are split into two separate types.
  ! The partitioning between ice and water is regarded as
  ! determining the areas within the grid_box covered by
  ! ice or water cloud, rather than as determining the in-cloud
  ! mixing ratios. The grid-box mean ice water contents may
  ! be predicted by the ice microphysics scheme or may be
  ! determined as a function of the temperature (LSP_FOCWWIL).

  ! Set the components within the clouds. Here we have two
  ! components: stratiform ice and water.
  n_condensed=2
  type_condensed(1)=ip_clcmp_st_water
  type_condensed(2)=ip_clcmp_st_ice


  ! Convective cloud is not considered here, so we set to zero:
  DO i=id_ct, npd_layer
    DO l=1, npd_profile
      condensed_mix_rat_area(l, i, ip_clcmp_cnv_water) = 0.0e+00
      condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)   = 0.0e+00
      condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_water) = 0.0e+00
      condensed_mix_rat_bulk(l, i, ip_clcmp_cnv_ice)   = 0.0e+00
    END DO      ! L
  END DO      ! I

  ! Stratiform clouds:
  DO i=1, nclds
    j = n_layer+1-i
    DO l=1, n_profile
      lg=i_gather(l)
      IF (lca_area(lg, i) > EPSILON(lca_area)) THEN
        IF (l_add_cca_to_mcica) THEN
          IF (cca(lg, i) > 0.0 .AND. ccw(lg, i) > 0.0) THEN
            condensed_mix_rat_area(l, j, ip_clcmp_st_water) =                  &
              (lccwc1(lg, i)+lccwc2(lg, i)+cca(lg, i)*ccw(lg, i))              &
              /( cca(lg, i) + (1.0-cca(lg, i)) * lca_area(lg, i) )
          ELSE
            condensed_mix_rat_area(l, j, ip_clcmp_st_water) =                  &
              (lccwc1(lg, i)+lccwc2(lg, i))/lca_area(lg, i)
          END IF
        ELSE
          condensed_mix_rat_area(l, j, ip_clcmp_st_water) =                    &
            (lccwc1(lg, i)+lccwc2(lg, i))/lca_area(lg, i)
        END IF
      ELSE
        IF (l_add_cca_to_mcica) THEN
          IF (cca(lg, i) > 0.0 .AND. ccw(lg, i) > 0.0) THEN
            condensed_mix_rat_area(l, j, ip_clcmp_st_water) = ccw(lg, i)
          ELSE
            condensed_mix_rat_area(l, j, ip_clcmp_st_water) = 0.0
          END IF
        ELSE
          condensed_mix_rat_area(l, j, ip_clcmp_st_water) = 0.0
        END IF
      END IF
      condensed_mix_rat_area(l, j, ip_clcmp_st_ice) =                          &
        condensed_mix_rat_area(l, j, ip_clcmp_st_water)

      IF (lca_bulk(lg, i) > EPSILON(lca_bulk)) THEN
        condensed_mix_rat_bulk(l, j, ip_clcmp_st_water) =                      &
          (lccwc1(lg, i)+lccwc2(lg, i))/lca_bulk(lg, i)
      ELSE
        IF (l_add_cca_to_mcica) THEN
          IF (cca(lg, i) > 0.0 .AND. ccw(lg, i) > 0.0) THEN
            condensed_mix_rat_bulk(l, j, ip_clcmp_st_water) = ccw(lg, i)
          ELSE
            condensed_mix_rat_bulk(l, j, ip_clcmp_st_water) = 0.0
          END IF
        ELSE
          condensed_mix_rat_bulk(l, j, ip_clcmp_st_water) = 0.0
        END IF
      END IF
      condensed_mix_rat_bulk(l, j, ip_clcmp_st_ice) =                          &
        condensed_mix_rat_bulk(l, j, ip_clcmp_st_water)
    END DO
  END DO

  !
  ! Cloud fractions:
  !
  DO i=n_layer+1-nclds, n_layer
    j = n_layer+1-i
    DO l=1, n_profile
      lg = i_gather(l)

      work_qcl(l,j) = 0.0
      work_qcf(l,j) = 0.0

      w_cloud(l, i) = lca_area(lg, j)
      IF (l_add_cca_to_mcica) THEN
        IF (cca(lg, j) > 0.0 .AND. ccw(lg, j) > 0.0) THEN
          w_cloud(l, i) = cca(lg, j) + (1.0-cca(lg, j)) * w_cloud(l, i)
        END IF
      END IF

      ! Partition clouds using the ratio of cloud water contents.
      IF (l_add_cca_to_mcica) THEN
        IF (cca(lg, j) > 0.0 .AND. ccw(lg, j) > 0.0 ) THEN
          liq_frac(l)  = 1.0 - ( (t(l,i)-starticeTKelvin) /                    &
                         (alliceTdegC-(starticeTKelvin-tm)) )
          liq_frac(l)=MIN(MAX(0.0,liq_frac(l)),1.0)
          ! Convective core diagnostics, before stratiform contribution
          ! is added to liq_frac
          IF (diag%l_ccore_clt_rad) THEN
            diag%ccore_clt_rad(col_list(l),row_list(l),j) = cca(lg, j)
          END IF

          IF (diag%l_ccore_qcl_rad .OR. diag%l_ccore_qcl_rad_path) THEN
            work_qcl(l,j) = cca(lg,j)*ccw(lg,j)*liq_frac(l)
          END IF

          IF (diag%l_ccore_qcf_rad .OR. diag%l_ccore_qcf_rad_path) THEN
            work_qcf(l,j) = cca(lg,j)*ccw(lg,j)*(1.0 - liq_frac(l))
          END IF

          liq_frac(l) =                                                        &
            ( lccwc1(lg,j) + cca(lg,j)*ccw(lg,j)*liq_frac(l) ) /               &
            ( lccwc1(lg,j) + cca(lg,j)*ccw(lg,j) + lccwc2(lg,j) )
        ELSE IF ( lca_area(lg, j) > EPSILON(lca_area) .AND.                    &
                  (lccwc1(lg, j) + lccwc2(lg, j)) > 0.0 ) THEN
          liq_frac(l) =  lccwc1(lg, j) / (lccwc1(lg, j) + lccwc2(lg, j))
        ELSE
          liq_frac(l) = 0.0e+00
        END IF
      ELSE IF ( lca_area(lg, j) > EPSILON(lca_area) .AND.                      &
                (lccwc1(lg, j) + lccwc2(lg, j)) > 0.0 ) THEN
        liq_frac(l) =  lccwc1(lg, j) / (lccwc1(lg, j) + lccwc2(lg, j))
      ELSE
        liq_frac(l) = 0.0e+00
      END IF

      frac_cloud(l, i, ip_cloud_type_sw) = liq_frac(l)
      frac_cloud(l, i, ip_cloud_type_si) = (1.0e+00 - liq_frac(l))
      frac_cloud(l, i, ip_cloud_type_cw) = 0.0e+00
      frac_cloud(l, i, ip_cloud_type_ci) = 0.0e+00

    END DO      ! L (N_PROFILE)
  END DO      ! I (N_LAYER)

  ! Copy work_qcl to ccore_qcl_rad 
  IF (diag%l_ccore_qcl_rad) THEN
    DO i=1,nclds
      DO l=1, n_profile
        diag%ccore_qcl_rad(col_list(l),row_list(l),i) = work_qcl(l,i)
      END DO
    END DO
  END IF

  ! Convective core liquid water path
  IF (diag%l_ccore_qcl_rad_path) THEN
    DO i = 1, nclds
      DO l=1, n_profile
          diag%ccore_qcl_rad_path(col_list(l),row_list(l)) =               & 
               diag%ccore_qcl_rad_path(col_list(l),row_list(l)) +          &
               work_qcl(l,i) * d_mass(l,i) 
      END DO
    END DO
  END IF

  ! Copy work_qcf to ccore_qcf_rad 
  IF (diag%l_ccore_qcf_rad) THEN
    DO i=1,nclds
      DO l=1, n_profile
        diag%ccore_qcf_rad(col_list(l),row_list(l),i) = work_qcf(l,i)
      END DO
    END DO
  END IF

  ! Convective core ice water path
  IF (diag%l_ccore_qcf_rad_path) THEN
    DO i = 1, nclds
      DO l=1, n_profile
          diag%ccore_qcf_rad_path(col_list(l),row_list(l)) =               & 
               diag%ccore_qcf_rad_path(col_list(l),row_list(l)) +          &
               work_qcf(l,i) * d_mass(l,i) 
      END DO
    END DO
  END IF

  IF (lrad_ccrad .AND. blending_option > off .AND.                         &
      .NOT. l_add_cca_to_mcica) THEN

    !        Add on convective cloud water
    !        If large-scale cloud exists assume that is partially resolving
    !        the convection and so only add on that part of ccw that doesn't
    !        overlap with lca

    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        convective_cloud_layer(l) = cca(lg,n_layer+1-i)
        IF ( convective_cloud_layer(l) > w_cloud(l, i)+1.0e-04 ) THEN
          condensed_mix_rat_area(l, i, ip_clcmp_st_water)              &
            = condensed_mix_rat_area(l, i, ip_clcmp_st_water) +        &
              ccw(lg,n_layer+1-i)*(cca(lg,n_layer+1-i)-w_cloud(l, i))  &
                                 / cca(lg,n_layer+1-i)
          condensed_mix_rat_area(l, i, ip_clcmp_st_ice)                &
            = condensed_mix_rat_area(l, i, ip_clcmp_st_water)
          condensed_mix_rat_bulk(l, i, ip_clcmp_st_water)              &
            = condensed_mix_rat_bulk(l, i, ip_clcmp_st_water) +        &
              ccw(lg,n_layer+1-i)*(cca(lg,n_layer+1-i)-w_cloud(l, i))  &
                                 / cca(lg,n_layer+1-i)
          condensed_mix_rat_bulk(l, i, ip_clcmp_st_ice)                &
            = condensed_mix_rat_bulk(l, i, ip_clcmp_st_water)
        END IF
      END DO
    END DO

    !        Add on convective cloud fraction, maximally overlapped
    !        Use lrad_ccrad and l_local_cnv_partition code where there
    !        is no large-scale cloud
    !
    !           PARTITION BETWEEN ICE AND WATER USING THE RELATIONSHIPS
    !           GIVEN IN BOWER ET AL. (1996, Q.J. 122 p 1815-1844). ICE
    !           IS ALLOWED IN A LAYER WARMER THAN THE FREEZING POINT
    !           ONLY IF THE TOP OF THE CLOUD IS GLACIATED.

    DO l=1, n_profile
      IF ((n_layer+2-cct(i_gather(l))) <= n_layer) THEN
        l_glaciated_top(l) =                                    &
          (t(l, MAX(n_layer+2-cct(i_gather(l))                  &
                  , n_layer-nclds+1)) <  tm)
      END IF
    END DO

    DO i=n_layer+1-nclds, n_layer

      DO l=1, n_profile
        lg = i_gather(l)
        convective_cloud_layer(l) = cca(lg,n_layer+1-i)

        ! Keep large scale cloud partition between ice and liquid
        ! (ie as calculated above) where it exists and use standard
        ! convective cloud partition where there is no large-scale

        IF ( convective_cloud_layer(l) > w_cloud(l, i)+1.0e-04     &
                              .AND. w_cloud(l, i) < 1.0e-04  ) THEN
          ! Convective cloud exists (using same test as used above
          ! for CCW) but no large-scale so use convective partition

          ! NOTE: Ice is allowed above the freezing point ONLY if
          !       the TOP of the cloud is glaciated.
          liq_frac_conv(l) = MAX(  0.0e+00                         &
                                 , MIN(  1.0e+00                   &
                                       , 1.61e-02*(t(l, i)-tm)     &
                                         + 8.9e-01))
          ! Do not allow ice above 0 Celsius unless the top of
          ! the cloud is glaciated and force homogeneous
          ! nucleation at -40 Celsius
          IF ( (t(l, i) > tm) .AND.                                &
               (.NOT. l_glaciated_top(l)) ) THEN
            liq_frac_conv(l) = 1.0e+00
          ELSE IF (t(l, i) < tm-4.0e+01) THEN
            liq_frac_conv(l) = 0.0e+00
          END IF

          frac_cloud(l,i,ip_cloud_type_sw) = liq_frac_conv(l)

          frac_cloud(l, i, ip_cloud_type_si)                         &
                    = (1.0e+00 - liq_frac_conv(l))
        END IF  ! test on no large-scale cloud

        ! Update total cloud cover in gridbox
        w_cloud(l, i) = MAX( w_cloud(l, i),                          &
                             convective_cloud_layer(l) )


      END DO      ! L (N_PROFILE)
    END DO      ! I (N_LAYER)
    !--------------------------------------------------------------------
    ! Cloud fraction/Mixing Ratio consistency check
    !--------------------------------------------------------------------
    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        lg = i_gather(l)
        IF (frac_cloud(l, i, ip_cloud_type_sw) <= 1.0e-10)             &
                 condensed_mix_rat_area(l, i, ip_clcmp_st_water)  = 0.0
        IF (frac_cloud(l, i, ip_cloud_type_si) <= 1.0e-10)             &
                 condensed_mix_rat_area(l, i, ip_clcmp_st_ice)    = 0.0
      END DO      ! L (N_PROFILE)
    END DO      ! I (N_LAYER)

  END IF   ! test on lrad_ccrad and blending_option

END IF      ! I_CLOUD_REPRESENTATION


!     EFFECTIVE RADII OF WATER CLOUDS: A MICROPHYSICAL PARAMETRIZATION
!     IS AVAILABLE; OTHERWISE STANDARD VALUES ARE USED.

IF (l_microphysics) THEN

  !        STANDARD VALUES ARE USED FOR ICE CRYSTALS, BUT
  !        A PARAMETRIZATION PROVIDED BY UMIST AND MRF
  !        IS USED FOR DROPLETS.

  !        CALCULATE THE DENSITY OF AIR.
  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      density_air(l, i)=p(l, i)/(r*t(l, i))
    END DO
  END DO


  IF ((.NOT. lrad_ccrad) .AND.                                   &
      (i_cloud_representation /= ip_cloud_ice_water)) THEN
    ! Original code. With lrad_ccrad = T this has been moved into
    ! R2_RE_MRF_UMIST

    DO l=1, n_profile
      cc_depth(l) = 0.0e+00
    END DO

    DO l=1, n_profile
      lg = i_gather(l)

      ! This loop should be safe even when convective cloud is not
      ! present, since CCB should not exceed CCT.
      DO i=n_layer+2-cct(lg), MIN(n_layer+1-ccb(lg),n_layer)
        cc_depth(l) = cc_depth(l) + (d_mass(l, i)/density_air(l, i))
      END DO
    END DO
  END IF      ! lrad_ccrad


  DO l=1, n_profile
    lying_snow_g(l)=lying_snow(i_gather(l))
  END DO

  IF (diag%nc_diag_flag) THEN
    DO i=n_layer+1-nclds, n_layer
      DO l=1, n_profile
        IF (t(l, i)  >=  tm) THEN
          strat_liq_cloud_fraction(l, i)=w_cloud(l, i)       &
                      *frac_cloud(l, i, ip_cloud_type_sw)
          conv_liq_cloud_fraction(l, i)=w_cloud(l, i)        &
                      *frac_cloud(l, i, ip_cloud_type_cw)
        ELSE
          strat_liq_cloud_fraction(l, i)=0.0
          conv_liq_cloud_fraction(l, i)=0.0
        END IF
      END DO
    END DO
  END IF

  CALL r2_re_mrf_umist(n_profile, n_layer, nclds                 &
     , i_gather, col_list, row_list                              &
     , l_aerosol_ccn                                             &
     , l_biomass_ccn, l_ocff_ccn, l_nitrate_ccn                  &
     , sea_salt_film, sea_salt_jet                               &
     , l_seasalt_ccn, salt_dim1, salt_dim2, salt_dim3            &
     , l_use_biogenic, biogenic, biogenic_dim1, biogenic_dim2    &
     , accum_sulphate, diss_sulphate, aitken_sulphate            &
     , aged_bmass, cloud_bmass                                   &
     , aged_ocff, cloud_ocff                                     &
     , accum_nitrate, diss_nitrate                               &
     , l_easyaerosol_cdnc, easyaerosol_cdnc                      &
     , lying_snow_g, i_cloud_representation                      &
     , land_g, flandg_g, density_air                             &
     , condensed_mix_rat_bulk, cc_depth                          &
     , n_drop_pot, condensed_dim_char                            &
     , CDNC                                                      &
     , d_mass                                                    &
     , strat_liq_cloud_fraction                                  &
     , conv_liq_cloud_fraction                                   &
     , diag%nc_diag_flag                                         &
     , nc_diag_g                                                 &
     , nc_weight_g                                               &
     , ntot_diag_g                                               &
     , strat_lwc_diag_g                                          &
     , so4_ccn_diag_g                                            &
     , sulp_dim1, sulp_dim2                                      &
     , bmass_dim1, bmass_dim2                                    &
     , ocff_dim1, ocff_dim2                                      &
     , nitrate_dim1, nitrate_dim2                                &
     , npd_field, npd_profile, npd_layer, npd_aerosol_species    &
     , npd_cloud_component, id_ct                                &
     )

  ! Constrain the sizes of droplets to lie within the range of
  ! validity of the parametrization scheme.
  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l, i, ip_clcmp_st_water)              &
         =MAX(condensed_min_dim(ip_clcmp_st_water)             &
         , MIN(condensed_max_dim(ip_clcmp_st_water)            &
         , condensed_dim_char(l, i, ip_clcmp_st_water)))
      condensed_dim_char(l, i, ip_clcmp_cnv_water)             &
         =MAX(condensed_min_dim(ip_clcmp_cnv_water)            &
         , MIN(condensed_max_dim(ip_clcmp_cnv_water)           &
         , condensed_dim_char(l, i, ip_clcmp_cnv_water)))
    END DO
  END DO


  ! Set microphysical diagnostics. Weights for cloud calculated
  ! here are used solely for the microphysics and do not have
  ! an independent meaning.

  IF (diag%wgt_conv_flag) THEN
    IF (i_cloud_representation == ip_cloud_conv_strat) THEN
      DO i=1, nclds
        DO l=1, n_profile
          diag%wgt_conv(col_list(l), row_list(l), i)                           &
             =w_cloud(l, n_layer+1-i)                                          &
             *frac_cloud(l, n_layer+1-i, ip_cloud_type_conv)
        END DO
      END DO
    ELSE IF (i_cloud_representation == ip_cloud_csiw) THEN
      DO i=1, nclds
        DO l=1, n_profile
          diag%wgt_conv(col_list(l), row_list(l), i)                           &
             =w_cloud(l, n_layer+1-i)                                          &
             *frac_cloud(l, n_layer+1-i, ip_cloud_type_cw)
        END DO
      END DO
    ELSE IF (i_cloud_representation == ip_cloud_ice_water) THEN
      DO i=1, nclds
        DO l=1, n_profile
          diag%wgt_conv(col_list(l), row_list(l), i)=0.0e+00
        END DO
      END DO
    END IF
  END IF

  IF (diag%re_conv_flag) THEN
    DO i=1, nclds
      DO l=1, n_profile
        ! Effective radii are given in microns.
        diag%re_conv(col_list(l), row_list(l), i)                              &
           =condensed_dim_char(l, n_layer+1-i, ip_clcmp_cnv_water)             &
           *diag%wgt_conv(col_list(l), row_list(l), i)*1.0e+06
      END DO
    END DO
  END IF

  IF (diag%wgt_strat_flag) THEN
    IF (i_cloud_representation == ip_cloud_conv_strat) THEN
      DO i=1, nclds
        DO l=1, n_profile
          diag%wgt_strat(col_list(l), row_list(l), i)                          &
             =w_cloud(l, n_layer+1-i)                                          &
             *frac_cloud(l, n_layer+1-i, ip_cloud_type_strat)
        END DO
      END DO
    ELSE IF ((i_cloud_representation == ip_cloud_csiw) .OR.                    &
            (i_cloud_representation == ip_cloud_ice_water)) THEN
      DO i=1, nclds
        DO l=1, n_profile
          diag%wgt_strat(col_list(l), row_list(l), i)                          &
             =w_cloud(l, n_layer+1-i)                                          &
             *frac_cloud(l, n_layer+1-i, ip_cloud_type_sw)
        END DO
      END DO
    END IF
  END IF

  IF (diag%re_strat_flag) THEN
    DO i=1, nclds
      DO l=1, n_profile
        ! Effective radii are given in microns.
        diag%re_strat(col_list(l), row_list(l), i)                             &
           =condensed_dim_char(l, n_layer+1-i, ip_clcmp_st_water)              &
           *diag%wgt_strat(col_list(l), row_list(l), i)*1.0e+06
      END DO
    END DO
  END IF

  IF (diag%lwp_strat_flag) THEN
    DO i=1, nclds
      DO l=1, n_profile
        diag%lwp_strat(col_list(l), row_list(l), i)                            &
           =condensed_mix_rat_area(l, n_layer+1-i, ip_clcmp_st_water)          &
           *d_mass(l, n_layer+1-i)                                             &
           *diag%wgt_strat(col_list(l), row_list(l), i)
      END DO
    END DO
  END IF

  IF (diag%nc_diag_flag .AND. diag%nc_weight_flag) THEN
    DO l=1, n_profile
      diag%nc_diag(col_list(l), row_list(l)) = nc_diag_g(l)*nc_weight_g(l)
    END DO
  END IF

  IF (diag%nc_weight_flag) THEN
    DO l=1, n_profile
      diag%nc_weight(col_list(l), row_list(l)) = nc_weight_g(l)
    END DO
  END IF

  IF (diag%ntot_diag_flag) THEN
    DO i=1, nclds
      DO l=1, n_profile
        diag%ntot_diag(col_list(l), row_list(l), i)                            &
           =ntot_diag_g(l, n_layer+1-i)                                        &
           *diag%wgt_strat(col_list(l), row_list(l), i)
      END DO
    END DO
  END IF

  IF (diag%strat_lwc_diag_flag) THEN
    DO i=1, nclds
      DO l=1, n_profile
        diag%strat_lwc_diag(col_list(l), row_list(l), i)                       &
           =strat_lwc_diag_g(l, n_layer+1-i)                                   &
           *diag%wgt_strat(col_list(l), row_list(l), i)
      END DO
    END DO
  END IF

  ! Non-cloud diagnostics are "weighted" by the conditional sampling
  ! weight COND_SAMP_WGT, but as this is 1.0 if the SW radiation is
  ! active, and 0.0 if it is not, there is no need to actually
  ! multiply by it.

  IF (diag%cond_samp_wgt_flag) THEN
    DO i=1, nclds
      DO l=1, n_profile
        diag%cond_samp_wgt(col_list(l), row_list(l), i)=1.0
      END DO
    END DO
  END IF

  IF (diag%so4_ccn_diag_flag) THEN
    DO i=1, nclds
      DO l=1, n_profile
        diag%so4_ccn_diag(col_list(l), row_list(l), i)                         &
          = so4_ccn_diag_g(l, n_layer+1-i)
      END DO
    END DO
  END IF

ELSE

  ! All effective radii are set to standard values.

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l, i, ip_clcmp_st_water)=7.0e-6
      condensed_dim_char(l, i, ip_clcmp_cnv_water)=7.0e-6
    END DO
  END DO

END IF

!     SET THE CHARACTERISTIC DIMENSIONS OF ICE CRYSTALS:

!     ICE CRYSTALS IN STRATIFORM CLOUDS:

SELECT CASE (i_condensed_param(ip_clcmp_st_ice))

CASE (ip_slingo_schrecker_ice)

  !        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
  !        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l, i, ip_clcmp_st_ice)=30.0e-6
    END DO
  END DO

CASE (ip_ice_adt)

  !        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
  !        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
  !        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
  !        AT THE FREEZING LEVEL.

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l, i, ip_clcmp_st_ice)                &
         =MIN(7.198755e-04                                     &
         , EXP(5.522e-02*(t(l, i)-2.7965e+02))/9.702e+02)
    END DO
  END DO

CASE (ip_ice_agg_de, ip_ice_agg_de_sun)

  !      Aggregate parametrization based on effective dimension.
  !      In the initial form, the same approach is used for stratiform
  !      and convective cloud.

  !      The fit provided here is based on Stephan Havemann's fit of
  !      Dge with temperature, consistent with David Mitchell's treatment
  !      of the variation of the size distribution with temperature. The
  !      parametrization of the optical properties is based on De
  !      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
  !      (=(2*SQRT(3)/3)*volume/projected area), which explains the
  !      conversion factor. The fit to Dge is in two sections, because
  !      Mitchell's relationship predicts a cusp at 216.208 K. Limits
  !      of 8 and 124 microns are imposed on Dge: these are based on this
  !      relationship and should be reviewed if it is changed. Note also
  !      that the relationship given here is for polycrystals only.
  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      !          Preliminary calculation of Dge.
      IF (t(l, i) < t_switch) THEN
        condensed_dim_char(l, i, ip_clcmp_st_ice)                  &
          = a0_agg_cold*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_cold
      ELSE
        condensed_dim_char(l, i, ip_clcmp_st_ice)                  &
          = a0_agg_warm*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_warm
      END IF
      !          Limit and convert to De.
      condensed_dim_char(l, i, ip_clcmp_st_ice)                    &
        = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
          MIN(1.24e-04, MAX(8.0e-06,                               &
          condensed_dim_char(l, i, ip_clcmp_st_ice)))
    END DO
  END DO


CASE (ip_ice_sun_fu, ip_ice_chou_vis)

  !       Effective size of ice crystals from Sun & Rikus scheme (Sun,
  !       2001: Reply to comments by Greg M. McFarquhar on 'Parametrization
  !       of effective radius of cirrus clouds and its verification against
  !       observations'. Q. J. Royal. Meteor. Soc., 127, pp 267)

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile

      !           IWC is in g m-3
      iwc=condensed_mix_rat_area(l, i, ip_clcmp_st_ice)           &
        *density_air(l,i)*1.0e+03 + 1.0e-50

      !           1.e-6 for converting to metre for consistency
      condensed_dim_char(l, i, ip_clcmp_st_ice)                   &
        =1.0e-6*(EXP(LOG(a1_s)+a2_s*LOG(iwc))                      &
        +EXP(LOG(a3_s)+a4_s*LOG(iwc))*(t(l,i)-a5_s) )

      condensed_dim_char(l, i, ip_clcmp_st_ice)                   &
        = MERGE( condensed_dim_char(l, i, ip_clcmp_st_ice)        &
        * (a6_s+a7_s*(t(l,i)-273.15))                             &
        , condensed_dim_char(l, i, ip_clcmp_st_ice)               &
        , (ABS (alat(l)) >= a8_s )  )

    END DO
  END DO

CASE (ip_ice_baran)

  !        Diagnostic De=De(T) parameterisation developed by A. Baran.
  !        Fit to data based on Field et al., JAS, 2007, 10.1175/2007JAS2344.1.
  !        The weighting of the ensemble model follows that assumed in
  !        experiment 4 of Baran et al. (2014 and 2016). The weights applied
  !        at each size bin are 0.50,0.20,0.30, 0.0, 0.0, 0.0 (that is
  !        hexagons+rosettes account for 70% of the mass at each size bin).
  !        De is defined as 1.50*(mass of PSD/area of PSD).
  !        The size is limited to its value at the freezing level.

  condensed_min_dim(ip_clcmp_st_ice) = min_ice_baran
  condensed_max_dim(ip_clcmp_st_ice) = max_ice_baran
  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l,i,ip_clcmp_st_ice)=m_ice_baran*t(l, i)-n_ice_baran
    END DO
  END DO

END SELECT ! I_CONDENSED_PARAM(IP_CLCMP_ST_ICE)


!     ICE CRYSTALS IN CONVECTIVE CLOUDS:

SELECT CASE (i_condensed_param(ip_clcmp_cnv_ice))

CASE (ip_slingo_schrecker_ice)

  !        THIS PARAMETRIZATION IS BASED ON THE EFFECTIVE RADIUS
  !        AND A STANDARD VALUE OF 30-MICRONS IS ASSUMED.

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l, i, ip_clcmp_cnv_ice)=30.0e-6
    END DO
  END DO

CASE (ip_ice_adt)

  !        THIS PARAMETRIZATION IS BASED ON THE MEAN MAXIMUM
  !        DIMENSION OF THE CRYSTAL, DETERMINED AS A FUNCTION OF
  !        THE LOCAL TEMPERATURE. THE SIZE IS LIMITED TO ITS VALUE
  !        AT THE FREEZING LEVEL.

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l, i, ip_clcmp_cnv_ice)               &
         =MIN(7.198755e-04                                     &
         , EXP(5.522e-02*(t(l, i)-2.7965e+02))/9.702e+02)
    END DO
  END DO

CASE (ip_ice_agg_de, ip_ice_agg_de_sun)

  !      Aggregate parametrization based on effective dimension.
  !      In the initial form, the same approach is used for stratiform
  !      and convective cloud.

  !      The fit provided here is based on Stephan Havemann's fit of
  !      Dge with temperature, consistent with David Mitchell's treatment
  !      of the variation of the size distribution with temperature. The
  !      parametrization of the optical properties is based on De
  !      (=(3/2)volume/projected area), whereas Stephan's fit gives Dge
  !      (=(2*SQRT(3)/3)*volume/projected area), which explains the
  !      conversion factor. The fit to Dge is in two sections, because
  !      Mitchell's relationship predicts a cusp at 216.208 K. Limits
  !      of 8 and 124 microns are imposed on Dge: these are based on this
  !      relationship and should be reviewed if it is changed. Note also
  !      that the relationship given here is for polycrystals only.
  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      !          Preliminary calculation of Dge.
      IF (t(l, i) < t_switch) THEN
        condensed_dim_char(l, i, ip_clcmp_cnv_ice)                 &
          = a0_agg_cold*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_cold
      ELSE
        condensed_dim_char(l, i, ip_clcmp_cnv_ice)                 &
          = a0_agg_warm*EXP(s0_agg*(t(l, i)-t0_agg))+b0_agg_warm
      END IF
      !          Limit and convert to De.
      condensed_dim_char(l, i, ip_clcmp_cnv_ice)                   &
        = (3.0/2.0)*(3.0/(2.0*SQRT(3.0)))*                         &
          MIN(1.24e-04, MAX(8.0e-06,                               &
          condensed_dim_char(l, i, ip_clcmp_cnv_ice)))
    END DO
  END DO


CASE (ip_ice_sun_fu, ip_ice_chou_vis)

  !       Effective size of ice crystals from Sun & Rikus scheme (Sun,
  !       2001: Reply to comments by Greg M. McFarquhar on 'Parametrization
  !       of effective radius of cirrus clouds and its verification against
  !       observations'. Q. J. Royal. Meteor. Soc., 127, pp 267)

  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile

      !           IWC is in g m-3
      iwc=condensed_mix_rat_area(l, i, ip_clcmp_cnv_ice)          &
        *density_air(l,i)*1.0e+03 +1.0e-50

      !           1.e-6 for converting to metre for consistency
      condensed_dim_char(l, i, ip_clcmp_cnv_ice)                  &
        =1.0e-6*(EXP(LOG(a1_s)+a2_s*LOG(iwc))                      &
        +EXP(LOG(a3_s)+a4_s*LOG(iwc))*(t(l,i)-a5_s) )

      condensed_dim_char(l, i, ip_clcmp_cnv_ice)                  &
        = MERGE( condensed_dim_char(l, i, ip_clcmp_cnv_ice)       &
        * (a6_s+a7_s*(t(l,i)-273.15))                             &
        , condensed_dim_char(l, i, ip_clcmp_cnv_ice)              &
        , (ABS (alat(l)) >= a8_s )  )

    END DO
  END DO

CASE (ip_ice_baran)

  !        Diagnostic De=De(T) parameterisation developed by A. Baran.
  !        Fit to data based on Field et al., JAS, 2007, 10.1175/2007JAS2344.1.
  !        The weighting of the ensemble model follows that assumed in
  !        experiment 4 of Baran et al. (2014 and 2016). The weights applied
  !        at each size bin are 0.50,0.20,0.30, 0.0, 0.0, 0.0 (that is
  !        hexagons+rosettes account for 70% of the mass at each size bin).
  !        De is defined as 1.50*(mass of PSD/area of PSD).
  !        The size is limited to its value at the freezing level.

  condensed_min_dim(ip_clcmp_cnv_ice) = min_ice_baran
  condensed_max_dim(ip_clcmp_cnv_ice) = max_ice_baran
  DO i=n_layer+1-nclds, n_layer
    DO l=1, n_profile
      condensed_dim_char(l,i,ip_clcmp_cnv_ice)=m_ice_baran*t(l, i)-n_ice_baran
    END DO
  END DO

END SELECT ! I_CONDENSED_PARAM(IP_CLCMP_CNV_ICE)



!     CONSTRAIN THE SIZES OF ICE CRYSTALS TO LIE WITHIN THE RANGE
!     OF VALIDITY OF THE PARAMETRIZATION SCHEME.
DO i=n_layer+1-nclds, n_layer
  DO l=1, n_profile
    condensed_dim_char(l, i, ip_clcmp_st_ice)                   &
       =MAX(condensed_min_dim(ip_clcmp_st_ice)                  &
       , MIN(condensed_max_dim(ip_clcmp_st_ice)                 &
       , condensed_dim_char(l, i, ip_clcmp_st_ice)))
    condensed_dim_char(l, i, ip_clcmp_cnv_ice)                  &
       =MAX(condensed_min_dim(ip_clcmp_cnv_ice)                 &
       , MIN(condensed_max_dim(ip_clcmp_cnv_ice)                &
       , condensed_dim_char(l, i, ip_clcmp_cnv_ice)))
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE r2_set_cloud_field

END MODULE r2_set_cloud_field_mod
