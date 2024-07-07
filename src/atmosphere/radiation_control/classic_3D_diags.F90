! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Calculates the 3D field of aerosol extinction coefficient, absorption, 
!  scattering, and gsca for CLASSIC aerosols.
!
MODULE classic_3D_diags_mod           

IMPLICIT NONE

! Description:
!   Outputs extinction, absorption and scattering as 3D fields
!   in units (m-1) for CLASSIC 3D aerosol-radiation diagnostics 
!
! Method:
!   For each radiatively active CLASSIC aerosol component the
!   specific scattering and absorption (units of m2/kg) originating 
!   from the spectral file are multiplied by mass mixing ratio and 
!   air density to calculate scattering, absorption and extinction
!   in units (m-1).
! 
! Output:  clas_aerosol_ext, clas_aerosol_abs
!          clas_aerosol_abs
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
! Code description: 
!   Language: Fortran 95. 
!   This code is written to UMDP3 standards. 

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CLASSIC_3D_DIAGS_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE classic_3D_diags (                                           &
! actual and fixed array dimensions
         n_aerosol_mr, npd_aerosol, npd_aerosol_mr,                     &
         nd_aod_wavel, n_profile, npd_profile,                          &
         n_layer, first_layer, npd_layer,                               &
         npd_humidities,                                                &
! variables with intent in
         type_aerosol, i_wavel,                                         &
         aerosol_mix_ratio, aerosol_mr_source, aerosol_mr_type_index,   &
         air_density, i_aerosol_parametrization,                        &
         aod_absorption, aod_scattering,                                &
         i_humidity_pointer, mean_rh,                                   &
         humidities, delta_humidity,                                    &
! variable with intent out
         clas_aerosol_ext, clas_aerosol_abs,                            &
         clas_aerosol_sca                                               &
       )

USE rad_pcf
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments

! Arguments with intent in

INTEGER, INTENT(IN) :: n_aerosol_mr
!  number of aerosol components in mixing ration information
INTEGER, INTENT(IN) :: npd_aerosol
!  number of aerosol components in the spectral information
INTEGER, INTENT(IN) :: npd_aerosol_mr
!  number of aerosol components in aerosol_mix_ratio
INTEGER, INTENT(IN) :: nd_aod_wavel        ! ... the AOD wavelengths
!  number of grid-boxes
INTEGER, INTENT(IN) :: n_profile
INTEGER, INTENT(IN) :: npd_profile
!  number of vertical layers
INTEGER, INTENT(IN) :: n_layer
INTEGER, INTENT(IN) :: first_layer
INTEGER, INTENT(IN) :: npd_layer
!  number of humidities used in moist aerosol parameterisation
INTEGER, INTENT(IN) :: npd_humidities
!  array giving the type of each aerosol component (see AERCMP3A)
INTEGER, INTENT(IN) :: type_aerosol(npd_aerosol)
! Index of wavelength to consider 
INTEGER, INTENT(IN) :: i_wavel
!  aerosol component mass mixing ratio
REAL, INTENT(IN) :: aerosol_mix_ratio(npd_profile, first_layer:npd_layer,      &
                       npd_aerosol_mr)
! Index relating aerosol_mix_ratio aerosols to aerosols in
! the spectral information
INTEGER, INTENT(IN) :: aerosol_mr_type_index(npd_aerosol_mr)
! Scheme/source of the aerosol data, to determine use in
! changing radiative fluxes and use in diagnostics
INTEGER, INTENT(IN) :: aerosol_mr_source(npd_aerosol_mr)
!  Density of the air (kg/m3)
REAL, INTENT(IN) :: air_density(npd_profile, npd_layer)
!  aerosol parameterisation (dry or moist)
INTEGER, INTENT(IN) :: i_aerosol_parametrization(npd_aerosol)
!  aerosol specific coefficients for absorption and scattering
!  (monochromatic)
REAL, INTENT(IN) :: aod_absorption(npd_humidities, npd_aerosol, nd_aod_wavel)
REAL, INTENT(IN) :: aod_scattering(npd_humidities, npd_aerosol, nd_aod_wavel)
!  mean relative humidities, and indices to the look-up tables
!  it may be the grid-box mean or clear-sky mean relative humidity,
!  depending on calculations made in FLUX_CALC
INTEGER, INTENT(IN) :: i_humidity_pointer(npd_profile, npd_layer)
REAL, INTENT(IN) :: mean_rh(npd_profile, npd_layer)
REAL, INTENT(IN) :: humidities(npd_humidities, npd_aerosol)
REAL, INTENT(IN) :: delta_humidity

! Arguments with intent out

!  computed aerosol extinction (m-1) summed over all aerosol species
REAL, INTENT(OUT) :: clas_aerosol_ext(npd_profile, npd_layer)
!  computed aerosol absorption (m-1) summed over all aerosol species
REAL, INTENT(OUT) :: clas_aerosol_abs(npd_profile, npd_layer)
!  computed aerosol scattering (m-1) summed over all aerosol species
REAL, INTENT(OUT) :: clas_aerosol_sca(npd_profile, npd_layer)

! Local variables

! Loop indices
INTEGER :: i, j, j_mr, k, l

! aerosol extinction (absorption + scattering)
REAL :: extinction
REAL :: absorption
REAL :: scattering

!  variables needed for the interpolation on humidities
INTEGER :: i_pointer
REAL :: weight_upper
REAL :: weight_lower

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CLASSIC_3D_DIAGS'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialize output arrays
DO l = 1, n_profile
  DO i = 1, n_layer
    clas_aerosol_ext(l, i) = 0.0
    clas_aerosol_abs(l, i) = 0.0
    clas_aerosol_sca(l, i) = 0.0
  END DO
END DO

! Begin large loop over CLASSIC aerosol components
DO j_mr = 1, n_aerosol_mr
  j=aerosol_mr_type_index(j_mr)

  ! DO NOT CONSIDER CLIMATOLOGICAL AEROSOLS (TYPE_AEROSOL IS
  ! SMALLER THAN 10)
  IF (type_aerosol(j) >= 10) THEN

    ! Only want aerosols that are radiatively active
    ! and do not care if they are prognostic or not
    IF (aerosol_mr_source(j_mr) == ip_aersrc_cusack_ron .OR.      &
      aerosol_mr_source(j_mr) == ip_aersrc_classic_ron  .OR.      &
      aerosol_mr_source(j_mr) == ip_aersrc_arcl_ron) THEN

      IF (i_aerosol_parametrization(j)  ==                        &
        ip_aerosol_param_dry) THEN

        ! Non-hygroscopic aerosol

        extinction = aod_absorption(1, j, i_wavel) +              &
                     aod_scattering(1, j, i_wavel)
        absorption = aod_absorption(1, j, i_wavel)
        scattering = aod_scattering(1, j, i_wavel)

        DO i = 1, n_layer
          DO l = 1, n_profile

            clas_aerosol_ext(l, i) = clas_aerosol_ext(l, i) +     &
              aerosol_mix_ratio(l, i, j_mr) *                     &
              air_density(l, i) * extinction

            clas_aerosol_abs(l, i) = clas_aerosol_abs(l, i) +     &
              aerosol_mix_ratio(l, i, j_mr) *                     &
              air_density(l, i) * absorption

            clas_aerosol_sca(l, i) = clas_aerosol_sca(l, i) +     &
              aerosol_mix_ratio(l, i, j_mr) *                     &
              air_density(l, i) * scattering

          END DO ! L
        END DO ! I

      ELSE IF (i_aerosol_parametrization(j)  ==                   &
            ip_aerosol_param_moist) THEN

      ! Hygroscopic aerosol
      ! interpolation on the mean relative humidity

        DO i = 1, n_layer
          DO l = 1, n_profile
            i_pointer = i_humidity_pointer(l, i)
            weight_upper = ( mean_rh(l, i)                        &
                   - humidities(i_pointer, j))                    &
                   / delta_humidity
            weight_lower = 1.00e+00 - weight_upper

            extinction =                                          &
              (aod_absorption(i_pointer, j, i_wavel) +            &
               aod_scattering(i_pointer, j, i_wavel))             &
              * weight_lower + weight_upper *                     &
              (aod_absorption(i_pointer+1,j,i_wavel) +            &
               aod_scattering(i_pointer+1,j,i_wavel))

            absorption =                                          &
              (aod_absorption(i_pointer, j, i_wavel)              &
              * weight_lower) + (weight_upper *                   &
              aod_absorption(i_pointer+1, j, i_wavel))

            scattering =                                          &
              (aod_scattering(i_pointer, j, i_wavel)              &
              * weight_lower) + (weight_upper *                   &
              aod_scattering(i_pointer+1, j, i_wavel))

            clas_aerosol_ext(l, i) = clas_aerosol_ext(l, i) +  &
              aerosol_mix_ratio(l, i, j_mr) *                        &
              air_density(l, i) * extinction

            clas_aerosol_abs(l, i) = clas_aerosol_abs(l, i) +  &
              aerosol_mix_ratio(l, i, j_mr) *                        &
              air_density(l, i) * absorption

            clas_aerosol_sca(l, i) = clas_aerosol_sca(l, i) +  &
              aerosol_mix_ratio(l, i, j_mr) *                        &
              air_density(l, i) * scattering

          END DO ! L 
        END DO ! I 
      END IF ! Hygroscopic / dry aerosol
    END IF ! Radiatively active aerosol sources
  END IF ! CLASSIC prognostic aerosol types
END DO ! J_MR

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE classic_3D_diags

END MODULE classic_3D_diags_mod
