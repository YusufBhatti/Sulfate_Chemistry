! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE compute_all_aod_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'COMPUTE_ALL_AOD_MOD'
CONTAINS

SUBROUTINE compute_all_aod(diag, row_list, col_list,                    &
       n_aerosol,      n_aerosol_mr,                                    &
       nd_aerosol_species, nd_aerosol_mixratio,                         &
       n_aod_wavel,    nd_aod_wavel,    aod_wavel,                      &
       n_profile,      nd_profile,                                      &
       n_layer,        first_layer,     nd_layer,                       &
       n_humidities,   nd_humidities,                                   &
       type_aerosol,   i_aod_type,                                      &
       aerosol_mix_ratio, d_mass, air_density,                          &
       aerosol_mr_source, aerosol_mr_type_index,                        &
       i_aerosol_parametrization, aod_absorption, aod_scattering,       &
       i_humidity_pointer, mean_rel_humidity,                           &
       humidities,     delta_humidity,   l_use_arcl)

USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE dust_parameters_mod, ONLY: l_twobin_dust
USE arcl_mod, ONLY: npd_arcl_species, ip_arcl_dust
USE def_diag, ONLY: StrDiag
USE classic_3d_diags_mod, ONLY: classic_3d_diags
USE aodtype_mod, ONLY: ip_type_allaod, ip_type_sulphate, ip_type_dust,   &
                       ip_type_seasalt, ip_type_soot, ip_type_biomass,   &
                       ip_type_biogenic, ip_type_ocff, ip_type_delta,    &
                       ip_type_nitrate, ip_type_twobdust
USE compute_aod_mod, ONLY: compute_aod
IMPLICIT NONE

! A wrapper subroutine to calculate the aerosol optical depth diagnostic
! a number of different times for different aerosols
!
! Description:
!   Acts as a wrapper subroutine so that compute_aod can be called multiple
!   times for different aerosol types.
!
! Method:
!   Calls the compute_aod subroutine multiple times, for different aerosol
!   types.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.
!
! Declarations:
!   These are of the form:-
!     INTEGER, INTENT(IN) :: ExampleVariable      !Description of variable
!
! Global variables (#include statements etc):

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='COMPUTE_ALL_AOD'

!     General properties of input data
INTEGER, INTENT(IN) :: n_aod_wavel  ! Number of wavelengths for aod calcualtions
INTEGER, INTENT(IN) :: n_aerosol    ! Number of aerosols in spectral information
INTEGER, INTENT(IN) :: n_aerosol_mr ! Number of aerosols in aerosol_mix_ratio array
INTEGER, INTENT(IN) :: n_profile    ! Number of Atmospheric profiles
INTEGER, INTENT(IN) :: n_layer      ! Number of Atmospheric layers
INTEGER, INTENT(IN) :: first_layer  ! First atmosphere layer
INTEGER, INTENT(IN) :: n_humidities ! Number of humidities used in computing aod
!
!     Size allocated for...
INTEGER, INTENT(IN) :: nd_profile   ! ... atmospheric profiles
INTEGER, INTENT(IN) :: nd_layer     ! ... atmospheric layers
INTEGER, INTENT(IN) :: nd_aerosol_species  ! ... aerosol species in spectral info
INTEGER, INTENT(IN) :: nd_aerosol_mixratio ! ... aerosol species in mixing array
INTEGER, INTENT(IN) :: nd_aod_wavel        ! ... the AOD wavelengths
INTEGER, INTENT(IN) :: nd_humidities       ! ... humidities for computing AOD

!     Input data about aerosol properties
INTEGER, INTENT(IN) :: type_aerosol(nd_aerosol_species)
                                        ! Aerosol type in spectral information
INTEGER, INTENT(IN) :: i_aod_type(nd_aerosol_species)
                                        ! Map between aerosol type and aod type
INTEGER, INTENT(IN) :: i_aerosol_parametrization(nd_aerosol_species)
                                            ! Parametrization flags for aerosol
INTEGER, INTENT(IN) :: aerosol_mr_type_index(nd_aerosol_mixratio)
                        ! Index relating aerosol_mix_ratio aerosols to aerosols in
                        ! the spectral information
INTEGER, INTENT(IN) :: aerosol_mr_source(nd_aerosol_mixratio)
                        ! Scheme/source of the aerosol data, to determine use in
                        ! changing radiative fluxes and use in diagnostics


REAL, INTENT(IN) :: aod_absorption(nd_humidities, nd_aerosol_species,  &
                       nd_aod_wavel)    ! Aerosol absorption array for AOD
REAL, INTENT(IN) :: aod_scattering(nd_humidities, nd_aerosol_species,  &
                       nd_aod_wavel)    ! Aerosol scattering array for AOD
REAL, INTENT(IN) :: humidities(nd_humidities, nd_aerosol_species)
                                        ! Humidities for aerpsol species
LOGICAL, INTENT(IN) :: l_use_arcl(npd_arcl_species)
                                        ! Logical for aerosol clim. species

!     Input data itself:
INTEGER, INTENT(IN):: i_humidity_pointer(nd_profile, nd_layer)
                                        ! Pointer to look-up table for aerosols
REAL, INTENT(IN) ::  aerosol_mix_ratio(nd_profile, nd_layer,         &
                       nd_aerosol_mixratio) ! Mixing ratios of aerosols
REAL, INTENT(IN) ::  d_mass(nd_profile, nd_layer) ! Mass thickness of each layer
REAL, INTENT(IN) ::  air_density(nd_profile, nd_layer) ! Air density (kg/m3)
REAL, INTENT(IN) ::  mean_rel_humidity(nd_profile, nd_layer)
                                            ! Mean relative humidity of layers
REAL, INTENT(IN) ::  delta_humidity         ! Increment in look-up table for hum.

REAL, INTENT(IN) ::  aod_wavel(nd_aod_wavel) ! Wavelengths for AOD calculations

! Diagnostics:
TYPE (StrDiag), INTENT(INOUT) :: diag

INTEGER, INTENT(IN) :: row_list(nd_profile)
!                          list of row indices of lit points
INTEGER, INTENT(IN) :: col_list(nd_profile)
!                          list of column indices of lit points

! AOD data to be calculated:
REAL :: aod_sulphate(nd_profile, nd_aod_wavel)
REAL :: aod_dust(nd_profile, nd_aod_wavel)
REAL :: aod_seasalt(nd_profile, nd_aod_wavel)
REAL :: aod_soot(nd_profile, nd_aod_wavel)
REAL :: aod_biomass(nd_profile, nd_aod_wavel)
REAL :: aod_biogenic(nd_profile, nd_aod_wavel)
REAL :: aod_ocff(nd_profile, nd_aod_wavel)
REAL :: aod_delta(nd_profile, nd_aod_wavel)
REAL :: aod_nitrate(nd_profile, nd_aod_wavel)
REAL :: aod_total_radn(nd_profile, nd_aod_wavel)
REAL :: angst_total_radn(nd_profile, nd_aod_wavel)
REAL :: aod_prog_sulphate(nd_profile, nd_aod_wavel)
REAL :: aod_prog_dust(nd_profile, nd_aod_wavel)
REAL :: aod_prog_seasalt(nd_profile, nd_aod_wavel)
REAL :: aod_prog_soot(nd_profile, nd_aod_wavel)
REAL :: aod_prog_biomass(nd_profile, nd_aod_wavel)
REAL :: aod_prog_ocff(nd_profile, nd_aod_wavel)
REAL :: aod_prog_nitrate(nd_profile, nd_aod_wavel)

! AAOD data to be calculated:
REAL :: aaod_sulphate(nd_profile, nd_aod_wavel)
REAL :: aaod_dust(nd_profile, nd_aod_wavel)
REAL :: aaod_seasalt(nd_profile, nd_aod_wavel)
REAL :: aaod_soot(nd_profile, nd_aod_wavel)
REAL :: aaod_biomass(nd_profile, nd_aod_wavel)
REAL :: aaod_biogenic(nd_profile, nd_aod_wavel)
REAL :: aaod_ocff(nd_profile, nd_aod_wavel)
REAL :: aaod_nitrate(nd_profile, nd_aod_wavel)
REAL :: aaod_dummy(nd_profile, nd_aod_wavel)

! 3D extinction diagnostic
REAL ::                                                            &
  clas_aerosol_ext(nd_profile, nd_layer)

! 3D absorption diagnostic
REAL ::                                                            &
  clas_aerosol_abs(nd_profile, nd_layer)

! 3D scattering diagnostic
REAL ::                                                            &
  clas_aerosol_sca(nd_profile, nd_layer)

!
LOGICAL :: l_radn_on
!  logical determining if the aod to be calculated is only from
!  radiatively active aerosol mixing ratio data

LOGICAL :: l_prognostic
!  logical determining if the aod to be calculated is only from
!  prognostic aerosol mixing ratio data, as opposed to climatologies

!  loop indices
INTEGER :: i, j, k, l, iwv_ext, iwv_abs, iwv_sca

! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!     For each aerosol type, compute the optical depth if
!     it was requested.
!
IF (diag%l_aod_sulphate .OR. diag%l_aaod_sulphate) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_sulphate, l_radn_on,l_prognostic, &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_sulphate, aaod_sulphate)
  IF (diag%l_aod_sulphate) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_sulphate(col_list(l), row_list(l), i) = aod_sulphate(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_sulphate) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_sulphate(col_list(l), row_list(l), i) = aaod_sulphate(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_dust .OR. diag%l_aaod_dust) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  IF (l_twobin_dust .AND. .NOT. l_use_arcl(ip_arcl_dust) ) THEN
     ! two-bin prognostic dust in radiation
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_twobdust,l_radn_on,l_prognostic,&
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_dust, aaod_dust)
  ELSE  ! not two-bin dust
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_dust, l_radn_on,l_prognostic,   &
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_dust, aaod_dust)
  END IF   ! ends not two-bin dust
  IF (diag%l_aod_dust) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_dust(col_list(l), row_list(l), i) = aod_dust(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_dust) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_dust(col_list(l), row_list(l), i) = aaod_dust(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_seasalt .OR. diag%l_aaod_seasalt) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_seasalt, l_radn_on,l_prognostic,  &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_seasalt, aaod_seasalt)
  IF (diag%l_aod_seasalt) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_seasalt(col_list(l), row_list(l), i) = aod_seasalt(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_seasalt) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_seasalt(col_list(l), row_list(l), i) = aaod_seasalt(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_soot .OR. diag%l_aaod_soot) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_soot, l_radn_on,l_prognostic,     &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_soot, aaod_soot)
  IF (diag%l_aod_soot) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_soot(col_list(l), row_list(l), i) = aod_soot(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_soot) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_soot(col_list(l), row_list(l), i) = aaod_soot(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_biomass .OR. diag%l_aaod_biomass) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_biomass, l_radn_on,l_prognostic,  &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_biomass, aaod_biomass)
  IF (diag%l_aod_biomass) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_biomass(col_list(l), row_list(l), i) = aod_biomass(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_biomass) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_biomass(col_list(l), row_list(l), i) = aaod_biomass(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_biogenic .OR. diag%l_aaod_biogenic) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_biogenic,l_radn_on,l_prognostic,  &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_biogenic, aaod_biogenic)
  IF (diag%l_aod_biogenic) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_biogenic(col_list(l), row_list(l), i) = aod_biogenic(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_biogenic) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_biogenic(col_list(l), row_list(l), i) = aaod_biogenic(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_ocff .OR. diag%l_aaod_ocff) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_ocff, l_radn_on,l_prognostic,     &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_ocff, aaod_ocff)
  IF (diag%l_aod_ocff) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_ocff(col_list(l), row_list(l), i) = aod_ocff(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_ocff) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_ocff(col_list(l), row_list(l), i) = aaod_ocff(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_delta) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_delta, l_radn_on,l_prognostic,    &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_delta, aaod_dummy)
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_delta(col_list(l), row_list(l), i) = aod_delta(l, i)
    END DO
  END DO
END IF
IF (diag%l_aod_nitrate .OR. diag%l_aaod_nitrate) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_nitrate,l_radn_on,l_prognostic,   &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_nitrate, aaod_nitrate)
  IF (diag%l_aod_nitrate) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_nitrate(col_list(l), row_list(l), i) = aod_nitrate(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_aaod_nitrate) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aaod_nitrate(col_list(l), row_list(l), i) = aaod_nitrate(l, i)
      END DO
    END DO
  END IF
END IF
IF (diag%l_aod_total_radn .OR. diag%l_angst_total_radn) THEN
  l_radn_on = .TRUE.
  l_prognostic = .FALSE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_allaod, l_radn_on,l_prognostic,   &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_total_radn, aaod_dummy)
  IF (diag%l_aod_total_radn) THEN
    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%aod_total_radn(col_list(l), row_list(l), i) = aod_total_radn(l, i)
      END DO
    END DO
  END IF
  IF (diag%l_angst_total_radn) THEN
    ! calculate the angstrom coefficient for different aod wavelengths
    ! angst=alog(aot_w1/aot_w2)/alog(w2/w1)

    IF (n_aod_wavel > 1 ) THEN
      ! for now, set the aod wavelengths hardcoded here:

      ! first wavelength, non-centred gradient:
      k = 1
      DO l = 1, n_profile
        angst_total_radn(l, k)=LOG(aod_total_radn(l,k)/aod_total_radn(l,k+1))/ &
                               LOG(aod_wavel(k+1)/aod_wavel(k))
      END DO
      ! centre wavelengths, centred gradients:
      DO k=2, n_aod_wavel-1
        DO l = 1, n_profile
          angst_total_radn(l,k)=                                               &
            LOG(aod_total_radn(l,k-1)/aod_total_radn(l,k+1))/                  &
            LOG(aod_wavel(k+1)/aod_wavel(k-1))
        END DO
      END DO
      ! final wavelength, non-centred gradient:
      k = n_aod_wavel
      DO l = 1, n_profile
        angst_total_radn(l, k)=LOG(aod_total_radn(l,k-1)/aod_total_radn(l,k))/ &
                               LOG(aod_wavel(k)/aod_wavel(k-1))
      END DO

    ELSE
      DO k = 1, n_aod_wavel
        DO l = 1, n_profile
          angst_total_radn(l, k) = 0.0
        END DO
      END DO
    END IF

    DO i=1, n_aod_wavel
      DO l=1, n_profile
        diag%angst_total_radn(col_list(l), row_list(l), i)                     &
          = angst_total_radn(l, i)
      END DO
    END DO

  END IF
END IF
IF (diag%l_aod_prog_sulphate) THEN
  l_radn_on = .FALSE.
  l_prognostic = .TRUE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_sulphate, l_radn_on,l_prognostic, &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_sulphate, aaod_dummy)
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_prog_sulphate(col_list(l), row_list(l), i)                      &
        = aod_prog_sulphate(l, i)
    END DO
  END DO
END IF
IF (diag%l_aod_prog_dust) THEN
  l_radn_on = .FALSE.
  l_prognostic = .TRUE.
  IF (l_twobin_dust) THEN
     ! two-bin prognostic dust
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_twobdust,l_radn_on,l_prognostic,&
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_prog_dust, aaod_dummy)
  ELSE  ! not two-bin dust
    CALL compute_aod(                                                   &
      n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio, &
      n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                 &
      n_layer, first_layer, nd_layer,                                   &
      n_humidities, nd_humidities,                                      &
      type_aerosol, i_aod_type, ip_type_dust, l_radn_on,l_prognostic,   &
      aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,      &
      d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,&
      i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,&
      aod_prog_dust, aaod_dummy)
  END IF   ! ends not two-bin dust
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_prog_dust(col_list(l), row_list(l), i) = aod_prog_dust(l, i)
    END DO
  END DO
END IF
IF (diag%l_aod_prog_seasalt) THEN
  l_radn_on = .FALSE.
  l_prognostic = .TRUE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_seasalt, l_radn_on,l_prognostic,  &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_seasalt, aaod_dummy)
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_prog_seasalt(col_list(l), row_list(l), i)                       &
        = aod_prog_seasalt(l, i)
    END DO
  END DO
END IF
IF (diag%l_aod_prog_soot) THEN
  l_radn_on = .FALSE.
  l_prognostic = .TRUE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_soot, l_radn_on,l_prognostic,     &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_soot, aaod_dummy)
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_prog_soot(col_list(l), row_list(l), i) = aod_prog_soot(l, i)
    END DO
  END DO
END IF
IF (diag%l_aod_prog_biomass) THEN
  l_radn_on = .FALSE.
  l_prognostic = .TRUE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_biomass, l_radn_on,l_prognostic,  &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_biomass, aaod_dummy)
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_prog_biomass(col_list(l), row_list(l), i)                       &
        = aod_prog_biomass(l, i)
    END DO
  END DO
END IF
IF (diag%l_aod_prog_ocff) THEN
  l_radn_on = .FALSE.
  l_prognostic = .TRUE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_ocff, l_radn_on,l_prognostic,     &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_ocff, aaod_dummy)
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_prog_ocff(col_list(l), row_list(l), i)                          &
        = aod_prog_ocff(l, i)
    END DO
  END DO
END IF
IF (diag%l_aod_prog_nitrate) THEN
  l_radn_on = .FALSE.
  l_prognostic = .TRUE.
  CALL compute_aod(                                                     &
    n_aerosol, n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,   &
    n_aod_wavel, nd_aod_wavel, n_profile, nd_profile,                   &
    n_layer, first_layer, nd_layer,                                     &
    n_humidities, nd_humidities,                                        &
    type_aerosol, i_aod_type, ip_type_nitrate, l_radn_on,l_prognostic,  &
    aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,        &
    d_mass, i_aerosol_parametrization, aod_absorption, aod_scattering,  &
    i_humidity_pointer, mean_rel_humidity, humidities, delta_humidity,  &
    aod_prog_nitrate, aaod_dummy)
  DO i=1, n_aod_wavel
    DO l=1, n_profile
      diag%aod_prog_nitrate(col_list(l), row_list(l), i)                       &
        = aod_prog_nitrate(l, i)
    END DO
  END DO
END IF

! Counters for array indicies
iwv_ext = 0
iwv_abs = 0
iwv_sca = 0

DO i = 1, n_aod_wavel

  ! If requested then calculate 3D aerosol radiation diagnostics
  ! separately for each wavelength from classic_3D_diags

  IF (diag%l_clas_aerosol_ext(i)  .OR.                                    &
      diag%l_clas_aerosol_abs(i)  .OR.                                    &
      diag%l_clas_aerosol_sca(i)) THEN

    CALL classic_3D_diags (                                               &
           n_aerosol_mr, nd_aerosol_species, nd_aerosol_mixratio,         &
           nd_aod_wavel, n_profile, nd_profile,                           &
           n_layer, first_layer, nd_layer,                                &
           nd_humidities, type_aerosol,                                   &
           ! Index of wavelength to consider 
           i,                                                             &
           aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,   &
           air_density, i_aerosol_parametrization,                        &
           aod_absorption, aod_scattering,                                &
           i_humidity_pointer, mean_rel_humidity,                         &
           humidities, delta_humidity,                                    &
           clas_aerosol_ext, clas_aerosol_abs,                            &
           clas_aerosol_sca                                               &
         )

    ! Put requested 3D aerosol optical property diagnostics in diag structure

    IF (diag%l_clas_aerosol_ext(i)) THEN
      iwv_ext = iwv_ext + 1
      DO k=1, n_layer
        DO l=1, n_profile
          diag%clas_aerosol_ext(col_list(l), row_list(l), k, iwv_ext)      &
             = clas_aerosol_ext(l, n_layer+1-k)
        END DO
      END DO
    END IF
    IF (diag%l_clas_aerosol_abs(i)) THEN
      iwv_abs = iwv_abs + 1
      DO k=1, n_layer
        DO l=1, n_profile
          diag%clas_aerosol_abs(col_list(l), row_list(l), k, iwv_abs)     &
             = clas_aerosol_abs(l, n_layer+1-k)
        END DO
      END DO
    END IF
    IF (diag%l_clas_aerosol_sca(i)) THEN
      iwv_sca = iwv_sca + 1
      DO k=1, n_layer
        DO l=1, n_profile
          diag%clas_aerosol_sca(col_list(l), row_list(l), k, iwv_sca)      &
             = clas_aerosol_sca(l, n_layer+1-k)
        END DO
      END DO
    END IF
  END IF ! If requesting 3D diags for this wavelength 
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE compute_all_aod
END MODULE compute_all_aod_mod
