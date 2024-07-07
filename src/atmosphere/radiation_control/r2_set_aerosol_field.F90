! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!- ---------------------------------------------------------------------
MODULE r2_set_aerosol_field_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'R2_SET_AEROSOL_FIELD_MOD'
CONTAINS

SUBROUTINE r2_set_aerosol_field(ierr                              &
   , n_profile, nlevs, n_layer, n_aerosol, n_aerosol_mr           &
   , type_aerosol, i_gather, row_list, col_list, l_extra_top      &
   , l_climat_aerosol, l_clim_aero_hgt, L_HadGEM1_Clim_Aero       &
   , bl_depth, t, n_levels_bl, l_murk_rad, aero_meso              &
   , l_dust, l_use_dust, dust_dim1, dust_dim2                     &
   , dust_1, dust_2, dust_3, dust_4, dust_5, dust_6               &
   , l_use_biogenic, biogenic_dim1, biogenic_dim2, biogenic       &
   , l_sulpc_so2, l_use_sulpc_direct                              &
   , sulp_dim1, sulp_dim2                                         &
   , accum_sulphate, aitken_sulphate                              &
   , l_use_seasalt_direct, salt_dim1, salt_dim2, salt_dim3        &
   , sea_salt_film, sea_salt_jet, p                               &
   , l_soot, l_use_soot_direct, soot_dim1, soot_dim2              &
   , fresh_soot, aged_soot                                        &
   , l_biomass, l_use_bmass_direct, bmass_dim1, bmass_dim2        &
   , fresh_bmass, aged_bmass                                      &
   , l_ocff, l_use_ocff_direct, ocff_dim1, ocff_dim2              &
   , fresh_ocff, aged_ocff                                        &
   , l_nitrate, l_use_nitrate_direct, nitrate_dim1, nitrate_dim2  &
   , accum_nitrate                                                &
   , n_arcl_species, n_arcl_compnts, i_arcl_compnts               &
   , l_use_arcl, arcl_dim1, arcl_dim2, arcl                       &
   , land, lying_snow, pstar                                      &
   , p_layer_boundaries                                           &
   , trindx, alat, previous_time                                  &
   , aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index  &
   , npd_field, npd_profile, npd_layer, npd_aerosol_species       &
   , npd_aerosol_mixratio, first_layer                            &
   )

!  Subroutine to set fields of aerosols.
!
! Purpose:
!   The mixing ratios of aerosols are transferred to the large array.
!
! Method:
!   Straightforward.
!
USE planet_constants_mod, ONLY: r, &
  l_planet_aerosol, planet_aerosol_mmr
USE conversions_mod, ONLY: pi
USE dust_parameters_mod, ONLY: l_twobin_dust
USE arcl_mod,           ONLY: npd_arcl_compnts, npd_arcl_species, &
                        ip_arcl_sulp, ip_arcl_dust, ip_arcl_sslt, &
                        ip_arcl_blck, ip_arcl_biom, ip_arcl_ocff, &
                        ip_arcl_dlta, &
                        ip_arcl_sulp_ac, ip_arcl_sulp_ak, &
                        ip_arcl_sulp_di, ip_arcl_dust_b1, &
                        ip_arcl_dust_b2, ip_arcl_dust_b3, &
                        ip_arcl_dust_b4, ip_arcl_dust_b5, &
                        ip_arcl_dust_b6, ip_arcl_sslt_fi, &
                        ip_arcl_sslt_jt, ip_arcl_blck_fr, &
                        ip_arcl_blck_ag, ip_arcl_biom_fr, &
                        ip_arcl_biom_ag, ip_arcl_biom_ic, &
                        ip_arcl_ocff_fr, ip_arcl_ocff_ag, &
                        ip_arcl_ocff_ic, ip_arcl_dlta_dl
USE rad_input_mod, ONLY: aeroscl_csk_clim
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE rad_pcf, ONLY: npd_aerosol_component, ip_water_soluble,       &
                   ip_dust_like, ip_oceanic, ip_soot, ip_ash,     &
                   ip_sulphuric, ip_accum_sulphate,               &
                   ip_aitken_sulphate, ip_fresh_soot,             &
                   ip_aged_soot, ip_seasalt_film,                 &
                   ip_seasalt_jet, ip_dust_div1, ip_dust_div2,    &
                   ip_dust_div3, ip_dust_div4, ip_dust_div5,      &
                   ip_dust_div6, ip_biomass_1, ip_biomass_2,      &
                   ip_biogenic, ip_ocff_fresh, ip_ocff_aged,      &
                   ip_delta, ip_nitrate, ip_twobindust_1,         &
                   ip_twobindust_2, ip_aersrc_cusack_ron,         &
                   ip_aersrc_cusack_roff,                         &
                   ip_aersrc_classic_ron,                         &
                   ip_aersrc_classic_roff,                        &
                   ip_aersrc_arcl_ron, ip_aersrc_arcl_roff,       &
                   i_err_fatal
USE umPrintMgr
USE errormessagelength_mod, ONLY: errormessagelength
USE r2_set_aero_clim_hadcm3_mod, ONLY: r2_set_aero_clim_hadcm3
IMPLICIT NONE
!
!     DUMMY ARGUMENTS.
!
!     SIZES OF ARRAYS:
INTEGER ::                                                        &
          !, INTENT(IN)
     npd_field                                                    &
!             FIELD SIZE IN CALLING PROGRAM
         , npd_profile                                                  &
!             SIZE OF ARRAY OF PROFILES
         , npd_layer                                                    &
!             MAXIMUM NUMBER OF LAYERS
         , npd_aerosol_species                                          &
!             MAXIMUM NUMBER OF AEROSOL SPECIES IN SPECTRAL INFORMATION
         , npd_aerosol_mixratio                                         &
!             MAXIMUM NUMBER OF AEROSOL SPECIES MIXING RATIO INFORMATION
         , first_layer
!             First layer for some variables
!             0 for flux_calc(3A/C), 1 for radiance_calc(3Z)
!
INTEGER ::                                                        &
          !, INTENT(OUT)
     ierr
!             ERROR FLAG
!
!     ACTUAL SIZES USED:
INTEGER ::                                                        &
          !, INTENT(IN)
     n_profile                                                    &
!             NUMBER OF PROFILES
         , nlevs                                                        &
!             Number of layers used outside the radiation scheme
         , n_layer                                                      &
!             Number of layers seen by radiation
         , n_levels_bl                                                  &
!             Number of layers occupied by boundary-layer aerosol
!             if L_CLIM_AERO_HGT is false.
         , n_aerosol                                                    &
!             NUMBER OF AEROSOLS IN SPECTRAL FILE
         , n_aerosol_mr                                                 &
!             NUMBER OF AEROSOLS IN AEROSOL_MIXING_RATIO ARRAY
         , type_aerosol(npd_aerosol_species)
!             ACTUAL TYPES OF AEROSOLS
!
!     GATHERING ARRAY:
INTEGER ::                                                        &
          !, INTENT(IN)
     i_gather(npd_field)                                          &
!             LIST OF POINTS TO GATHER
         , row_list(npd_profile)                                  & 
!             list of row indices of lit points  
         , col_list(npd_profile)
!             list of column indices of lit points
!
!     FLAG FOR THE CLIMATOLOGICAL AEROSOL DISTRIBUTION.
LOGICAL ::                                                        &
             !, INTENT(IN)
     l_climat_aerosol                                             &
!             FLAG FOR CLIMATOLOGICAL AEROSOL DISTRIBUTION
         , l_clim_aero_hgt                                              &
!             Flag to use the boundary layer depth in setting the
!             climatological aerosol
         , L_HadGEM1_Clim_Aero                                          &
!             Flag to use HadGEM1 setting for climatological aerosols
         , l_murk_rad                                                   &
!             FLAG FOR MESOSCALE MODEL AEROSOL
         , l_extra_top
!             Flag to include an extra top layer in the radiation scheme
!
! Declare mineral dust aerosol variables:
LOGICAL ::                                                        &
     l_dust                                                       &
!             Mineral dust is available for radn. effects or diagnositcs
         , l_use_dust
!             Use direct radiative effect of mineral dust
INTEGER :: dust_dim1,dust_dim2
!        Dimensions for mineral dust arrays (P_FIELD,P_LEVELS or 1,1)
REAL :: dust_1(dust_dim1, dust_dim2)                                 &
                                  !Mass mixing ratio of div1 dust
   , dust_2(dust_dim1, dust_dim2)                                 &
                                  !Mass mixing ratio of div2 dust
   , dust_3(dust_dim1, dust_dim2)                                 &
                                  !Mass mixing ratio of div3 dust
   , dust_4(dust_dim1, dust_dim2)                                 &
                                  !Mass mixing ratio of div4 dust
   , dust_5(dust_dim1, dust_dim2)                                 &
                                  !Mass mixing ratio of div5 dust
   , dust_6(dust_dim1, dust_dim2) !Mass mixing ratio of div6 dust
!     VARIABLES FOR THE BIOGENIC AEROSOL
LOGICAL :: l_use_biogenic
INTEGER :: biogenic_dim1, biogenic_dim2
REAL :: biogenic(biogenic_dim1, biogenic_dim2)

!     VARIABLES FOR THE SULPHUR CYCLE:
LOGICAL ::                                                        &
             !, INTENT(IN)
     l_sulpc_so2                                                  &
!             SULPHUR CYCLE AVAILABLE FOR RADN. EFFECTS OR DIAGNOSTICS
         , l_use_sulpc_direct
!             FLAG TO USE SULPHUR CYCLE FOR DIRECT EFFECT
INTEGER ::                                                        &
             !, INTENT(IN)
     sulp_dim1,sulp_dim2
!             DIMENSIONS FOR _SULPHATE ARRAYS, (P_FIELD,P_LEVELS or 1,1)
REAL ::                                                           &
          !, INTENT(IN)
     accum_sulphate(sulp_dim1, sulp_dim2)                         &
!             MASS MIXING RATIOS OF ACCUMULATION MODE AEROSOL
         , aitken_sulphate(sulp_dim1, sulp_dim2)
!             MASS MIXING RATIOS OF AITKEN MODE AEROSOL
!
! Declare soot variables:
LOGICAL ::                                                        &
     l_soot                                                       &
!             SOOT AEROSOL AVAILABLE FOR RADN. EFFECTS OR DIAGNOSTICS
         , l_use_soot_direct
!             USE DIRECT RAD. EFFECT OF SOOT AEROSOL
INTEGER :: soot_dim1,soot_dim2
          !DIMENSIONS FOR SOOT ARRAYS, (P_FIELD,P_LEVELS or 1,1)
REAL :: fresh_soot(soot_dim1, soot_dim2)                             &
                                           ! MMR OF FRESH SOOT
   , aged_soot(soot_dim1, soot_dim2)       ! MMR OF AGED SOOT
!
! Declare biomass smoke aerosol variables:
LOGICAL ::                                                        &
     l_biomass                                                    &
!             Biomass smoke available for radn. effects or diagnostics
         , l_use_bmass_direct
!              Use direct radiative effect of biomass smoke
INTEGER :: bmass_dim1,bmass_dim2
!              Dimensions for biomass arrays (P_FIELD,P_LEVELS or 1,1)
REAL :: fresh_bmass(bmass_dim1, bmass_dim2)                          &
                                              ! MMR OF FRESH SMOKE
   , aged_bmass(bmass_dim1, bmass_dim2)       ! MMR OF AGED SMOKE
!
! Declare fossil-fuel organic carbon aerosol variables:
LOGICAL ::                                                        &
     l_ocff                                                       &
!              OCFF aerosol available for radn. effects or diagnostics
         , l_use_ocff_direct
!              Use direct radiative effect of OCFF aerosol
INTEGER :: ocff_dim1, ocff_dim2
!              Dimensions for OCFF arrays (P_FIELD,P_LEVELS or 1,1)
REAL :: fresh_ocff(ocff_dim1, ocff_dim2)                             &
                                              ! MMR OF FRESH OCFF
   , aged_ocff(ocff_dim1, ocff_dim2)          ! MMR OF AGED OCFF
!
LOGICAL ::                                                        &
     l_use_seasalt_direct
!             Flag for direct effect of interactive sea-salt aerosol
!
INTEGER ::                                                        &
     salt_dim1, salt_dim2, salt_dim3
!             Array sizes of sea-salt aerosols
!                (either P_FIELD,P_LEVELS or 1,1)
!
REAL ::                                                           &
     sea_salt_film(salt_dim1, salt_dim2, salt_dim3)               &
!             On input, number concentration (m-3) of film-mode
!             sea-salt aerosols; converted to mass mixing ratio.
         , sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)
!             On input, number concentration (m-3) of jet-mode
!             sea-salt aerosols; converted to mass mixing ratio.
! Declare nitrate aerosol variables:
LOGICAL ::                                                        &
     l_nitrate                                                    &
!              Nitrate aerosol available for radn. effects or diagnostcs
         , l_use_nitrate_direct
!              Use direct radiative effect of nitrate aerosols
INTEGER :: nitrate_dim1, nitrate_dim2
!              Dimensions for nitrate arrays (P_FIELD,P_LEVELS or 1,1)
REAL :: accum_nitrate(nitrate_dim1, nitrate_dim2)
!              Mass-mixing ratio of accumulation nitrate
!
!
! Aerosol climatology for NWP:
!

      ! Number of requested species within the climatology
INTEGER :: n_arcl_species

! Corresponding number of requested components
INTEGER :: n_arcl_compnts

! Model switches for each species
LOGICAL :: l_use_arcl(npd_arcl_species)

! Array index of components
INTEGER :: i_arcl_compnts(npd_arcl_compnts)

! Array dimensions
INTEGER :: arcl_dim1, arcl_dim2

! Mass-mixing ratios
REAL ::                                                           &
    arcl(arcl_dim1, arcl_dim2, n_arcl_compnts)

REAL ::                                                           &
     aero_meso(npd_field, nlevs)
!             MIXING RATIO OF 'URBAN' AEROSOL OF MESOSCALE MODEL
!     GENERAL ATMOSPHERIC PROPERTIES:
INTEGER ::                                                        &
          !, INTENT(IN)
     trindx(npd_field)                                            &
!             LAYER BOUNDARY OF TROPOPAUSE
      ,    previous_time(7)
!             Time
REAL ::                                                           &
          !, INTENT(IN)
     pstar(npd_field)                                             &
!             SURFACE PRESSURES
      ,    p_layer_boundaries(npd_field,0:nlevs)                        &
!             PRESSURE AT BOUNDARIES OF LAYERS
      ,    alat(npd_profile)
!             Latitudes in degrees
!
!     SURFACE FIELDS
LOGICAL ::                                                        &
          !, INTENT(IN)
     land(npd_field)
!             LAND SEA MASK
REAL ::                                                           &
          !, INTENT(IN)
     lying_snow(npd_field)                                        &
!             DEPTH OF LYING SNOW
         , bl_depth(npd_field)                                          &
!             Depth of the boundary layer
         , t(npd_profile, first_layer:npd_layer)                        &
!             Temperatures of atmospheric layers
         , p(npd_profile, first_layer:npd_layer)
!             Pressures at layer centres
!
REAL ::                                                           &
          !, INTENT(OUT)
     aerosol_mix_ratio(npd_profile, first_layer:npd_layer         &
        , npd_aerosol_mixratio)
!             MIXING RATIOS OF AEROSOLS
INTEGER ::                                                        &
             !, INTENT(OUT)
     aerosol_mr_type_index(npd_aerosol_mixratio)                  &
!             Index relating aerosol_mix_ratio aerosols to aerosols in
!             the spectral information
         , aerosol_mr_source(npd_aerosol_mixratio)
!             Source of the aerosol data (compared to IP_AERSRC_... )
!
!
!     LOCAL VARIABLES:
INTEGER ::                                                        &
     i                                                            &
!             LOOP VARIABLE
         , j                                                            &
!             LOOP VARIABLE
         , j_plus                                                       &
!             ADDITION TO LOOP VARIABLE, FOR INDEXING AEROSOL_MIX_RATIO
         , j_mr                                                         &
!             LOOP VARIABLE FOR LOOPING OVER J+J_PLUS
         , l                                                            &
!             LOOP VARIABLES
         , lg, lh                                                           &
!             INDEX FOR GATHERING
         , bltop                                                        &
!             INDEX OF UPPER BOUNDARY OF PLANETARY BOUNDARY LAYER
         , i_aerosol                                                    &
!             ACTUAL TYPE OF AEROSOL BEING CONSIDERED
         , i_top_copy
!             Topmost radiative layer to be set directly from the
!             profile supplied
!
!
!     ARRAYS FOR THE CLIMATOLOGICAL AEROSOL MODEL
INTEGER :: ii
!       Variable used in initialization of arrays.
LOGICAL :: l_in_climat(npd_aerosol_component)

!       Flags to indicate which aerosols are included in the
!       climatology: this may be used to enable various components
!       to be replaced by fully prognostic schemes.
INTEGER :: i_clim_pointer(npd_aerosol_component) =                &
    (/ 1, 2, 3, 4, 0, 5,                                          &
    (0, ii=7, npd_aerosol_component) /)

!       Pointers to the indices of the original climatological
!       aerosol model.
REAL ::                                                           &
     aerosol_mix_ratio_clim(npd_profile, first_layer:npd_layer, 5)
!             MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
!
REAL ::                                                           &
     mode_radius_ss_film                                          &
!            Mode radius of film-mode sea-salt aerosol (m)
         , mode_radius_ss_jet                                           &
!            Mode radius of jet-mode sea-salt aerosol (m)
         , sigma_ss_film                                                &
!            Geometric standard deviation of film-mode sea-salt aerosol
         , sigma_ss_jet                                                 &
!            Geometric standard deviation of jet-mode sea-salt aerosol
         , density_sea_salt
!            Bulk density of sodium chloride, taken as representative
!            of sea-salt aerosol (kg m-3); taken from Pruppacher & Klett
!            2nd edition (1997), p. 244.
!
INTEGER :: im
REAL :: Meso_frac(npd_aerosol_component) =                        &
  (/ 0.61, 0.17, 0.0, 0.22, (0.0, im=5,npd_aerosol_component) /)
REAL :: inv17

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*),PARAMETER   :: RoutineName = 'R2_SET_AEROSOL_FIELD'
CHARACTER(LEN=errormessagelength) :: cmessage

PARAMETER(                                                        &
     mode_radius_ss_film=0.1e-06                                  &
   , mode_radius_ss_jet=1.0e-06                                   &
   , sigma_ss_film=1.9                                            &
   , sigma_ss_jet=2.0                                             &
   , density_sea_salt=2165.0                                      &
   )
!
!

!
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!
!     INITIALISE SOME VALUES:
j_plus=0
!
!     If an extra layer is used in the radiation scheme, any
!     aerosol profiles supplied will be copied into the profiles
!     seen by the radiation scheme starting with the second layer
!     down, rather than the first.
i_top_copy=1
IF (l_extra_top) i_top_copy=2


l_in_climat(1:4) = .TRUE.
l_in_climat(5)   = .FALSE.
l_in_climat(6)   = .TRUE.
l_in_climat(7:npd_aerosol_component) = .FALSE.

!  Use HadGEM1 settings to remove boundary layer climatological aerosols.
IF (L_HadGEM1_Clim_Aero) THEN
  DO ii=1,n_aerosol
    IF ((type_aerosol(ii) == ip_water_soluble) .OR.                &
        (type_aerosol(ii) == ip_dust_like) .OR.                   &
        (type_aerosol(ii) == ip_oceanic) .OR.                      &
        (type_aerosol(ii) == ip_soot)) THEN
      l_in_climat(ii) = .FALSE.
    END IF
  END DO
  !  Check to see if MURK is on in radiation
  IF (l_murk_rad) THEN
    WRITE(cmessage, '(A)') '*** ERROR: MURK AEROSOL SWITCHED'      &
      // 'ON IN RADIATION BUT THIS IS CURRENTLY INCONSISTENT WITH' &
      // 'ONLY HAVING THE STRATOSPHERIC CUSACK CLIMATLOGY - THE'   &
      // 'MURK SPLITS ITS MIXING RATIO AMONGST THE BOUNDARY LAYER' &
      // 'CUSACK AEROSOLS'
    ierr=i_err_fatal
    CALL ereport(RoutineName,ierr,cmessage)
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF
END IF
!
!  Use climatological soot if climatological aerosols are on and not
!  using interactive soot.
l_in_climat(ip_soot) = l_in_climat(ip_soot)                      &
                    .AND. (.NOT. l_use_soot_direct)
!
!  Use climatological "oceanic" aerosol if climatological aerosols have
!  been selected and not using sea-salt parametrization.
l_in_climat(ip_oceanic) = l_in_climat(ip_oceanic)                 &
                          .AND. (.NOT. l_use_seasalt_direct)
!
!  The climatological water-soluble aerosol should not be used
!  if the direct effects of sulphate aerosols are included.
!  (Note that this differs the situation applying in earlier
!  versions of the model, such as HadAM3).
l_in_climat(ip_water_soluble) = l_in_climat(ip_water_soluble)     &
                     .AND. (.NOT. l_use_sulpc_direct)
!
!
!
IF (l_climat_aerosol) THEN
  !
  !        SET THE MIXING RATIOS OF THE CLIMATOLOGICAL AEROSOLS
  !        USED IN THE CLIMATOLOGY OF HADCM3. A SEPARATE SUBROUTINE
  !        IS USED TO ENSURE BIT-REPRODUCIBLE RESULTS BY USING
  !        EARLIER CODE. THIS COULD BE ALTERED IF A NEW CLIMATOLOGY WERE
  !        USED.
  !
  CALL r2_set_aero_clim_hadcm3(n_profile, nlevs, n_layer         &
     , i_gather, l_extra_top                                     &
     , l_clim_aero_hgt, bl_depth, t, n_levels_bl                 &
     , land, lying_snow, pstar, p_layer_boundaries, trindx       &
     , alat, previous_time                                       &
     , aerosol_mix_ratio_clim                                    &
     , npd_field, npd_profile, npd_layer, first_layer            &
     )
  !
END IF
!
!
!     THE AEROSOLS REQUIRED BY FOR THE CALCULATION SHOULD HAVE BEEN
!     SELECTED WHEN THE SPECTRAL FILE WAS READ IN. EACH TYPE SHOULD
!     BE SET APPROPRIATELY.
!
DO j=1, n_aerosol
  !
  i_aerosol=type_aerosol(j)

  IF (l_planet_aerosol) THEN
    aerosol_mr_type_index(j+j_plus)=j
    aerosol_mr_source(j+j_plus) = ip_aersrc_cusack_ron
    aerosol_mix_ratio(1:n_profile, 1:n_layer , j+j_plus) = &
      planet_aerosol_mmr

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE IF (l_climat_aerosol .AND. l_in_climat(i_aerosol)) THEN
    aerosol_mr_type_index(j+j_plus)=j
    aerosol_mr_source(j+j_plus) = ip_aersrc_cusack_ron
    !
    !           Here mxing ratios in all layers are set because the
    !           possible presence of an extra layer has been allowed
    !           for in the subroutine.
    DO i=1, n_layer
      DO l=1, n_profile
        aerosol_mix_ratio(l, i, j+j_plus)                     &
           =aerosol_mix_ratio_clim(l, i                       &
           , i_clim_pointer(i_aerosol))                       &
           * aeroscl_csk_clim(i_clim_pointer(i_aerosol))
      END DO
    END DO
    !
    ! For the other aerosol species, they may be supplied from dedicated
    ! schemes (e.g. climate model) or climatologies (e.g. NWP model).
    ! The former are indicated by L_USE_<SPECIES>, the latter by
    ! L_USE_ARCL(IP_ARCL_<SPECIES>). Should both logical be set to true
    ! for a given aerosol species (i.e. an aerosol is provided by a
    ! dedicated scheme AND a climatology), the climatology wins and will
    ! be seen by radiation.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     Mineral dust aerosols:
    !
    !        Six-bin dust, bin 1, clim or prognostic
  ELSE IF (i_aerosol == ip_dust_div1 .AND.                       &
           (l_use_arcl(ip_arcl_dust) .OR. ( l_dust .AND.           &
            (.NOT. l_twobin_dust) ) ) ) THEN
    IF (l_use_arcl(ip_arcl_dust)) THEN
       ! RECORD THE J FOR THE AEROSOL MIXING RATIO AS THIS IS THE
       ! INDEX THAT RELATES IT TO SPECTRAL INFORMATION
      aerosol_mr_type_index(j+j_plus)=j
      ! AEROSOL CLIM ALWAYS SEEN BY RADIATION
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      ! NOW SET THE AEROSOL DATA VALUES
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =arcl(lg, n_layer+1-i,                             &
                   i_arcl_compnts(ip_arcl_dust_b1))
        END DO
      END DO
      IF ( l_dust .AND. .NOT. l_twobin_dust ) THEN
         ! MAKE ROOM FOR THE 6 BIN PROGNOSTIC DUST ASWELL
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_dust .AND. (.NOT. l_twobin_dust) ) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_dust) ) THEN
         ! PROG. NOT USED TO EFFECT RAD FLUXES IF CLIM IS ON
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_dust) THEN
         ! DUST DIRECT EFFECT IS ON, SO USED IN RADIATION:
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =dust_1(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
    !        Two-bin dust, bin 1. prognostic only, 6 bin clim used in radn
  ELSE IF (i_aerosol == ip_twobindust_1 .AND.                    &
           l_dust .AND. l_twobin_dust ) THEN
    aerosol_mr_type_index(j+j_plus)=j
    IF ( l_use_arcl(ip_arcl_dust) ) THEN
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
    ELSE IF (l_use_dust) THEN
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
    ELSE
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
    END IF
    DO i=i_top_copy, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        aerosol_mix_ratio(l, i, j+j_plus)                     &
           =dust_1(lg, n_layer+1-i)
      END DO
    END DO

    !
    !        Six-bin dust, bin 2, clim or prognostic
  ELSE IF (i_aerosol == ip_dust_div2 .AND.                       &
           (l_use_arcl(ip_arcl_dust) .OR. ( l_dust .AND.           &
            (.NOT. l_twobin_dust) ) ) ) THEN
    IF (l_use_arcl(ip_arcl_dust)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =arcl(lg, n_layer+1-i,                             &
                   i_arcl_compnts(ip_arcl_dust_b2))
        END DO
      END DO
      IF ( l_dust .AND. .NOT. l_twobin_dust ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_dust .AND. (.NOT. l_twobin_dust) ) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_dust) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_dust) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =dust_2(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
    !        Two-bin dust, bin 2. prognostic only, 6 bin clim used in radn
  ELSE IF (i_aerosol == ip_twobindust_2 .AND.                    &
           l_dust .AND. l_twobin_dust ) THEN
    aerosol_mr_type_index(j+j_plus)=j
    IF ( l_use_arcl(ip_arcl_dust) ) THEN
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
    ELSE IF (l_use_dust) THEN
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
    ELSE
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
    END IF
    DO i=i_top_copy, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        aerosol_mix_ratio(l, i, j+j_plus)                     &
           =dust_2(lg, n_layer+1-i)
      END DO
    END DO
    !
    !        Six-bin dust, bin 3, clim or prognostic
  ELSE IF (i_aerosol == ip_dust_div3 .AND.                       &
           (l_use_arcl(ip_arcl_dust) .OR. ( l_dust .AND.           &
            (.NOT. l_twobin_dust) ) ) ) THEN
    IF (l_use_arcl(ip_arcl_dust)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =arcl(lg, n_layer+1-i,                             &
                   i_arcl_compnts(ip_arcl_dust_b3))
        END DO
      END DO
      IF ( l_dust .AND. .NOT. l_twobin_dust ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_dust .AND. (.NOT. l_twobin_dust) ) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_dust) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_dust) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =dust_3(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
    !        Six-bin dust, bin 4, clim or prognostic
  ELSE IF (i_aerosol == ip_dust_div4 .AND.                       &
           (l_use_arcl(ip_arcl_dust) .OR. ( l_dust .AND.           &
            (.NOT. l_twobin_dust) ) ) ) THEN
    IF (l_use_arcl(ip_arcl_dust)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =arcl(lg, n_layer+1-i,                             &
                   i_arcl_compnts(ip_arcl_dust_b4))
        END DO
      END DO
      IF ( l_dust .AND. .NOT. l_twobin_dust ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_dust .AND. (.NOT. l_twobin_dust) ) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_dust) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_dust) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =dust_4(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
    !        Six-bin dust, bin 5, clim or prognostic
  ELSE IF (i_aerosol == ip_dust_div5 .AND.                       &
           (l_use_arcl(ip_arcl_dust) .OR. ( l_dust .AND.           &
            (.NOT. l_twobin_dust) ) ) ) THEN
    IF (l_use_arcl(ip_arcl_dust)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =arcl(lg, n_layer+1-i,                             &
                   i_arcl_compnts(ip_arcl_dust_b5))
        END DO
      END DO
      IF ( l_dust .AND. .NOT. l_twobin_dust ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_dust .AND. (.NOT. l_twobin_dust) ) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_dust) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_dust) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =dust_5(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
    !        Six-bin dust, bin 6, clim or prognostic
  ELSE IF (i_aerosol == ip_dust_div6 .AND.                       &
           (l_use_arcl(ip_arcl_dust) .OR. ( l_dust .AND.           &
            (.NOT. l_twobin_dust) ) ) ) THEN
    IF (l_use_arcl(ip_arcl_dust)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =arcl(lg, n_layer+1-i,                             &
                   i_arcl_compnts(ip_arcl_dust_b6))
        END DO
      END DO
      IF ( l_dust .AND. .NOT. l_twobin_dust ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_dust .AND. (.NOT. l_twobin_dust) ) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_dust) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_dust) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                     &
             =dust_6(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Aerosols related to the sulphur cycle (note that dissolved
    !        sulphate does not contribute to the direct effect):
  ELSE IF (i_aerosol == ip_accum_sulphate .AND.                  &
           (l_use_arcl(ip_arcl_sulp) .OR. l_sulpc_so2)) THEN
    IF (l_use_arcl(ip_arcl_sulp)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =arcl(lg, n_layer+1-i,                             &
                  i_arcl_compnts(ip_arcl_sulp_ac))
        END DO
      END DO
      IF ( l_sulpc_so2 ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_sulpc_so2) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_sulp) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_sulpc_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =accum_sulphate(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
  ELSE IF (i_aerosol == ip_aitken_sulphate .AND.                 &
           (l_use_arcl(ip_arcl_sulp) .OR. l_sulpc_so2)) THEN

    IF (l_use_arcl(ip_arcl_sulp)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =arcl(lg, n_layer+1-i,                             &
                  i_arcl_compnts(ip_arcl_sulp_ak))
        END DO
      END DO
      IF ( l_sulpc_so2 ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_sulpc_so2) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_sulp) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_sulpc_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =aitken_sulphate(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Black Carbon/Soot aerosols
  ELSE IF (i_aerosol == ip_fresh_soot .AND.                      &
           (l_use_arcl(ip_arcl_blck) .OR. l_soot)) THEN
    IF (l_use_arcl(ip_arcl_blck)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =arcl(lg, n_layer+1-i,                             &
                  i_arcl_compnts(ip_arcl_blck_fr))
        END DO
      END DO
      IF ( l_soot ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_soot) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_blck) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_soot_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =fresh_soot(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
  ELSE IF (i_aerosol == ip_aged_soot .AND.                       &
           (l_use_arcl(ip_arcl_blck) .OR. l_soot)) THEN
    IF (l_use_arcl(ip_arcl_blck)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =arcl(lg, n_layer+1-i,                             &
                  i_arcl_compnts(ip_arcl_blck_ag))
        END DO
      END DO
      IF ( l_soot ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_soot) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_blck) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_soot_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =aged_soot(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Biomass burning aerosols
  ELSE IF (i_aerosol == ip_biomass_1 .AND.                       &
           (l_use_arcl(ip_arcl_biom) .OR. l_biomass)) THEN

    IF (l_use_arcl(ip_arcl_biom)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =arcl(lg, n_layer+1-i,                             &
                  i_arcl_compnts(ip_arcl_biom_fr))
        END DO
      END DO
      IF ( l_biomass ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_biomass) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_biom) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_bmass_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =fresh_bmass(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
  ELSE IF (i_aerosol == ip_biomass_2 .AND.                       &
           (l_use_arcl(ip_arcl_biom) .OR. l_biomass)) THEN

    IF (l_use_arcl(ip_arcl_biom)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =arcl(lg, n_layer+1-i,                             &
                  i_arcl_compnts(ip_arcl_biom_ag))
        END DO
      END DO
      IF ( l_biomass ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_biomass) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_biom) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_bmass_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)                    &
            =aged_bmass(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Sea salt aerosols
  ELSE IF (i_aerosol == ip_seasalt_film .AND.                    &
         (l_use_arcl(ip_arcl_sslt) .OR. l_use_seasalt_direct)) THEN

    IF (l_use_arcl(ip_arcl_sslt)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                   &
                 (arcl(lg, n_layer+1-i,                        &
                  i_arcl_compnts(ip_arcl_sslt_fi))*4.0*pi      &
                *density_sea_salt*(mode_radius_ss_film**3)   &
                *EXP(4.5*(LOG(sigma_ss_film))**2))           &
                /(3.0*p(l, i)/(r*t(l, i)))
        END DO
      END DO
      IF ( l_use_seasalt_direct ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_use_seasalt_direct) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_sslt) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=col_list(l)
          lh=row_list(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                   &
                 (sea_salt_film(lg, lh, n_layer+1-i)*4.0*pi    &
                *density_sea_salt*(mode_radius_ss_film**3)     &
                *EXP(4.5*(LOG(sigma_ss_film))**2))             &
                /(3.0*p(l, i)/(r*t(l, i)))
        END DO
      END DO
    END IF
    !
  ELSE IF (i_aerosol == ip_seasalt_jet .AND.                     &
           (l_use_arcl(ip_arcl_sslt) .OR. l_use_seasalt_direct)) THEN

    IF (l_use_arcl(ip_arcl_sslt)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                   &
                 (arcl(lg, n_layer+1-i,                        &
                  i_arcl_compnts(ip_arcl_sslt_jt))*4.0*pi      &
                *density_sea_salt*(mode_radius_ss_jet**3)    &
                *EXP(4.5*(LOG(sigma_ss_jet))**2))            &
                /(3.0*p(l, i)/(r*t(l, i)))
        END DO
      END DO
      IF ( l_use_seasalt_direct ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_use_seasalt_direct) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_sslt) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=col_list(l)
          lh=row_list(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                   &
                 (sea_salt_jet(lg, lh, n_layer+1-i)*4.0*pi     & 
                *density_sea_salt*(mode_radius_ss_jet**3)      &
                *EXP(4.5*(LOG(sigma_ss_jet))**2))              &
                /(3.0*p(l, i)/(r*t(l, i)))
        END DO
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Organic Carbon from Fossil Fuel aerosol:
  ELSE IF (i_aerosol == ip_ocff_fresh .AND.                      &
           (l_use_arcl(ip_arcl_ocff) .OR. l_ocff)) THEN

    IF (l_use_arcl(ip_arcl_ocff)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                    &
                 arcl(lg, n_layer+1-i,                          &
                   i_arcl_compnts(ip_arcl_ocff_fr))
        END DO
      END DO
      IF ( l_ocff ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_ocff) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_ocff) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_ocff_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                     &
                 fresh_ocff(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !
  ELSE IF (i_aerosol == ip_ocff_aged .AND.                       &
           (l_use_arcl(ip_arcl_ocff) .OR. l_ocff)) THEN

    IF (l_use_arcl(ip_arcl_ocff)) THEN
      aerosol_mr_type_index(j+j_plus)=j
      aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                    &
                 arcl(lg, n_layer+1-i,                          &
                   i_arcl_compnts(ip_arcl_ocff_ag))
        END DO
      END DO
      IF ( l_ocff ) THEN
        j_plus=j_plus+1
      END IF
    END IF
    IF (l_ocff) THEN
      aerosol_mr_type_index(j+j_plus)=j
      IF ( l_use_arcl(ip_arcl_ocff) ) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      ELSE IF (l_use_ocff_direct) THEN
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
      ELSE
        aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
      END IF
      DO i=i_top_copy, n_layer
        DO l=1, n_profile
          lg=i_gather(l)
          aerosol_mix_ratio(l, i, j+j_plus)=                     &
                  aged_ocff(lg, n_layer+1-i)
        END DO
      END DO
    END IF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Nitrate aerosol: only available as prognostic
  ELSE IF ((i_aerosol == ip_nitrate) .AND. l_nitrate) THEN
    aerosol_mr_type_index(j+j_plus)=j
    IF (l_use_nitrate_direct) THEN
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_ron
    ELSE
      aerosol_mr_source(j+j_plus) = ip_aersrc_classic_roff
    END IF
    DO i=i_top_copy, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        aerosol_mix_ratio(l, i, j+j_plus)=                    &
                accum_nitrate(lg, n_layer+1-i)
      END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Biogenic aerosols, only a climatology:
  ELSE IF ((i_aerosol == ip_biogenic) .AND. l_use_biogenic) THEN
    aerosol_mr_type_index(j+j_plus)=j
    aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
    DO i=i_top_copy, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        aerosol_mix_ratio(l, i, j+j_plus)=                      &
                 biogenic(lg, n_layer+1-i)
      END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        'Delta' aerosol: only ever a climatology
  ELSE IF ( (i_aerosol == ip_delta) .AND.                        &
             l_use_arcl(ip_arcl_dlta) ) THEN
    aerosol_mr_type_index(j+j_plus)=j
    aerosol_mr_source(j+j_plus) = ip_aersrc_arcl_ron
    DO i=i_top_copy, n_layer
      DO l=1, n_profile
        lg=i_gather(l)
        aerosol_mix_ratio(l, i, j+j_plus)=                       &
          arcl(lg, n_layer+1-i, i_arcl_compnts(ip_arcl_dlta_dl))
      END DO
    END DO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE
    !           The options to the radiation code do not require this
    !           aerosol to be considered: its mixing ratio is not set.
    aerosol_mr_type_index(j+j_plus)=j
    aerosol_mr_source(j+j_plus) = ip_aersrc_cusack_roff
  END IF
  !
  !        If using an extra top layer extrapolate the mixing ratio
  !        from the adjacent layer, except in the case of Cusack
  !        climatological aerosols which are specifically set.
  IF (l_extra_top .AND.                                          &
     aerosol_mr_source(j+j_plus) /= ip_aersrc_cusack_ron .AND.   &
     aerosol_mr_source(j+j_plus) /= ip_aersrc_cusack_roff ) THEN
    DO l=1, n_profile
      aerosol_mix_ratio(l, 1, j+j_plus)                          &
        =aerosol_mix_ratio(l, 2, j+j_plus)
    END DO
  END IF
  !
END DO

! Check the number of aerosols actually set matches the number
! expected and report diagnostic info if the numbers mismatch
IF (n_aerosol_mr /= n_aerosol+j_plus) THEN
  ! writing details of the error to stdout:
  CALL umPrint( '*** ERROR: '                                            &
    // 'AN UNEXPECTED NUMBER OF AEROSOL SPECIES HAS BEEN SET.***'        &
    // '(details below)',src='r2_set_Aerosol_field')
  WRITE(umMessage, '(A,I3)') 'NUMBER OF AEROSOLS IN SPEC FILE:',         &
      n_aerosol
  CALL umPrint(umMessage,src='r2_set_Aerosol_field')
  WRITE(umMessage, '(A,I3)') 'NUMBER OF AEROSOLS EXPECTED (CLIMS'        &
      // ' AND PROGNSOTICS):' ,  n_aerosol_mr
  CALL umPrint(umMessage,src='r2_set_Aerosol_field')
  WRITE(umMessage, '(A,I3)') 'NUMBER OF AEROSOLS ACTUALLY SET:',         &
      n_aerosol+j_plus
  CALL umPrint(umMessage,src='r2_set_Aerosol_field')
  CALL umPrint(' OUTPUT AEROSOL ARRAYS FROM R2_SET_AEROSOL_FIELD: ',     &
      src='r2_set_Aerosol_field')
  CALL umPrint( ' (SEE rad_pcf.F90 TO INTERPRET THE'                     &
      // ' AEROSOL TYPES AND SOURCES)',src='r2_set_Aerosol_field')
  CALL umPrint( ' J_MR, AEROSOL_MR_TYPE_INDEX, '                         &
    // 'TYPE_AEROSOL(AEROSOL_MR_TYPE_INDEX), AEROSOL_MR_SOURCE',         &
    src='r2_set_Aerosol_field')
  DO j_mr=1, n_aerosol_mr
    IF (aerosol_mr_type_index(j_mr) >= 1 .AND.                   &
        aerosol_mr_type_index(j_mr) <= n_aerosol) THEN
      WRITE(umMessage,'(I3,I3,I3,I3)') j_mr,                            &
        aerosol_mr_type_index(j_mr),                                    &
        type_aerosol(aerosol_mr_type_index(j_mr)),                      &
        aerosol_mr_source(j_mr)
      CALL umPrint(umMessage,src='r2_set_Aerosol_field')
    ELSE
      WRITE(umMessage,'(I3,I3,A3,I3)') j_mr,                            &
        aerosol_mr_type_index(j_mr),' - ', aerosol_mr_source(j_mr)
      CALL umPrint(umMessage,src='r2_set_Aerosol_field')
    END IF
  END DO
  ! now writing out the error message itself
  WRITE(cmessage, '(A)') '*** ERROR: '                            &
    // 'AN UNEXPECTED NUMBER OF AEROSOL SPECIES HAS BEEN SET.*** '&
    // ' (see stdout for details)'
  ierr=i_err_fatal
  CALL ereport(RoutineName,ierr,cmessage)
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

!
! Add the MURK aerosol onto the Cusack climatology (via meso_frac),
! but only in the boundary layer levels
IF (l_murk_rad) THEN
  inv17=1.0/1.7
  DO j_mr=1, n_aerosol_mr
    DO i=(nlevs-n_levels_bl+1),nlevs
      DO l=1, n_profile
        lg=i_gather(l)
        ! Note that aerosol mass mixing ratios are in units of 1E9 kg/kg
        aerosol_mix_ratio(l,i,j_mr) = aerosol_mix_ratio(l,i,j_mr)   &
          + aero_meso(lg,nlevs+1-i)*1e-9*                           &
            meso_frac(aerosol_mr_type_index(j_mr))*inv17
      END DO
    END DO
  END DO
END IF
!
!

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE r2_set_aerosol_field
END MODULE r2_set_aerosol_field_mod
