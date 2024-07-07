! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to set the input atmospheric profiles for the core radiation code
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_atm_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_ATM_MOD'

CONTAINS

SUBROUTINE set_atm(                                                            &

! Structures for the core radiation code interface
  control, dimen, spectrum, atm,                                               &

! Grid
  n_profile, n_layer, nozone, nd_field, list, true_latitude, true_longitude,   &

! Thermodynamic fields
  p_layer_boundaries, p_layer_centres, t_layer_boundaries, t_layer_centres,    &
  p_extra_layer, t_extra_layer, d_mass, density,                               &
  r_layer_centres, r_layer_boundaries,                                         &

! Gas mass mixing ratios
  h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio,                   &
  cfc11_mix_ratio, cfc12_mix_ratio, o2_mix_ratio,                              &
  cfc113_mix_ratio, cfc114_mix_ratio, hcfc22_mix_ratio,                        &
  hfc125_mix_ratio, hfc134a_mix_ratio,                                         &

! 3D CO2
  co2_dim1, co2_dim2, co2_3d, l_co2_3d,                                        &

! Chemical greenhouse gas fields
  ngrgas, grgas_field,                                                         &

! Fields to calculate satellite viewing directions
  sindec, seconds_since_midnight)


USE rad_pcf
USE def_control,  ONLY: StrCtrl
USE def_dimen,    ONLY: StrDim
USE def_spectrum, ONLY: StrSpecData
USE def_atm,      ONLY: StrAtm, allocate_atm
USE yomhook,      ONLY: lhook, dr_hook
USE parkind1,     ONLY: jprb, jpim
USE ereport_mod,  ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE planet_constants_mod, ONLY: planet_radius, g, r, c_virtual
USE conversions_mod, ONLY: pi, rsec_per_day
USE nlsizes_namelist_mod, ONLY: model_levels
USE rad_input_mod, ONLY: l_extra_top

USE r2_set_gas_mix_ratio_mod, ONLY: r2_set_gas_mix_ratio
IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(IN)    :: control

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Atmospheric properties:
TYPE(StrAtm),       INTENT(INOUT) :: atm

INTEGER, INTENT(IN) :: n_profile
!   Number of atmospheric profiles for radiation calculations
INTEGER, INTENT(IN) :: n_layer
!   Number of atmospheric layers for radiation calculations
INTEGER, INTENT(IN) :: nozone
!   Number of levels with ozone
INTEGER, INTENT(IN) :: nd_field
!   Field size in calling program
INTEGER, INTENT(IN) :: list(nd_field)
!   List of points where radiation is to be calculated
REAL, INTENT(IN) :: true_latitude(nd_field)
!   Latitudes of grid-points
REAL, INTENT(IN) :: true_longitude(nd_field)
!   Longitudes of grid-points

REAL, INTENT(IN) :: p_layer_boundaries(nd_field, 0:model_levels)
!   Pressure at layer boundaries
REAL, INTENT(IN) :: p_layer_centres(nd_field, 0:model_levels)
!   Pressure at layer centres
REAL, INTENT(IN) :: t_layer_boundaries(nd_field, 0:model_levels)
!   Temperature at layer boundaries
REAL, INTENT(IN) :: t_layer_centres(nd_field, model_levels)
!   Temperature at layer centres
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

! Mixing ratios supplied:
REAL, INTENT(IN) :: h2o(nd_field, model_levels)
!   Mass mixing ratio of water vapour
REAL, INTENT(IN) :: co2
!   Mass mixing ratio of carbon dioxide
REAL, INTENT(IN) :: o3(nd_field, nozone)
!   Mass mixing ratio of ozone
REAL, INTENT(IN) :: n2o_mix_ratio
!   Mass mixing ratio of nitrous oxide
REAL, INTENT(IN) :: ch4_mix_ratio
!   Mass mixing ratio of methane
REAL, INTENT(IN) :: so2_mix_ratio
!   Mass mixing ratio of sulphur dioxide
REAL, INTENT(IN) :: cfc11_mix_ratio
!   Mass mixing ratio of CFC11
REAL, INTENT(IN) :: cfc12_mix_ratio
!   Mass mixing ratio of CFC12
REAL, INTENT(IN) :: o2_mix_ratio
!   Mass mixing ratio of O2
REAL, INTENT(IN) :: cfc113_mix_ratio
!   Mass mixing ratio of CFC113
REAL, INTENT(IN) :: cfc114_mix_ratio
!   Mass mixing ratio of CFC114
REAL, INTENT(IN) :: hcfc22_mix_ratio
!   Mass mixing ratio of HCFC22
REAL, INTENT(IN) :: hfc125_mix_ratio
!   Mass mixing ratio of HFC125
REAL, INTENT(IN) :: hfc134a_mix_ratio
!   Mass mixing ratio of HFC134a

! 3D CO2
INTEGER, INTENT(IN) :: co2_dim1, co2_dim2
!   Dimensions of CO2_3D field
REAL, INTENT(IN) :: co2_3d(co2_dim1, co2_dim2)
!   3D mass mixing ratio of CO2 (full field)
LOGICAL, INTENT(IN) :: l_co2_3d
!   Controls use of 3D co2 field

! Chemical greenhouse gas fields
INTEGER, INTENT(IN) :: ngrgas
REAL, INTENT(IN) :: grgas_field(nd_field, model_levels, ngrgas)

REAL, INTENT(IN) :: sindec
!   Sine of solar declination
REAL, INTENT(IN) :: seconds_since_midnight
!   Seconds since midnight



! Local variables.
INTEGER :: i, l, lg, i_top_copy
!   Loop variables

REAL :: r_sat(3)
!   Cartesian coordinates of the satellite
REAL :: r_obs(3)
!   Cartesian coordinates of observed point
REAL :: r_diff(3)
!   Cartesian coordinates of separation of the points
REAL :: r_sun(3)
!   Cartesian coordinates of the sun
REAL :: r_hor(3)
!   Projection of r_diff on local horizontal plane
REAL :: r_hor_s(3)
!   Projection of r_sun on local horizontal plane
REAL :: mag_rdiff
!   Magnitude of the separation
REAL :: recip_mag_robs
!   Reciprocal of the magnitude of r_obs
REAL :: mag_rhor
!   Magnitude or r_hor
REAL :: mag_rhor_s
!   Magnitude of r_hor_s
REAL :: solar_zen
!   Solar zenith angle
REAL :: h_angle
!   Hour Angle

INTEGER                      :: ierr = i_normal
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SET_ATM'
CHARACTER (LEN=errormessagelength) :: cmessage

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL allocate_atm(atm, dimen, spectrum)

! Setup atmosphere grid
atm%n_profile = n_profile
atm%n_layer   = n_layer
DO l=1, n_profile
  atm%lat(l) = true_latitude(list(l))
  atm%lon(l) = true_longitude(list(l))
END DO

! Set the pressures and temperatures
IF (l_extra_top) THEN
  ! Here, an extra level is inserted at the top of the atmosphere by replacing
  ! the top layer used in the rest of the model with two layers. The absolute
  ! top pressure level is set to zero.
  DO l=1, n_profile
    lg=list(l)
    atm%p_level(l, 0) = 0.0
    atm%t_level(l, 0) = 0.0
    atm%p_level(l, 1) = p_layer_boundaries(lg, model_levels)
    atm%t_level(l, 1) = t_layer_boundaries(lg, model_levels)
    atm%p(l, 1) = p_extra_layer(lg)
    atm%t(l, 1) = t_extra_layer(lg)
  END DO
  ! The second radiative layer will be the first to have properties
  ! set by copying input fields.
  i_top_copy=2
ELSE
  DO l=1, n_profile
    lg=list(l)
    atm%p_level(l, 0) = p_layer_boundaries(lg, model_levels)
    atm%t_level(l, 0) = t_layer_boundaries(lg, model_levels)
  END DO
  i_top_copy=1
END IF
DO i=i_top_copy, n_layer
  DO l=1, n_profile
    lg=list(l)
    atm%p_level(l, i) = p_layer_boundaries(lg, n_layer-i)
    atm%t_level(l, i) = t_layer_boundaries(lg, n_layer-i)
    atm%p(l, i)       = p_layer_centres(lg, n_layer+1-i)
    atm%t(l, i)       = t_layer_centres(lg, n_layer+1-i)
  END DO
END DO

! Set the layer masses (per square metre), densities, and heights
DO l=1, n_profile
  lg=list(l)
  atm%r_level(l, 0) = r_layer_boundaries(lg, n_layer)
END DO
DO i=1, n_layer
  DO l=1, n_profile
    lg=list(l)
    atm%mass(l, i)    = d_mass(lg, n_layer+1-i)
    atm%density(l, i) = density(lg, n_layer+1-i)
    atm%r_level(l, i) = r_layer_boundaries(lg, n_layer-i)
    atm%r_layer(l, i) = r_layer_centres(lg, n_layer+1-i)
  END DO
END DO

! Set the mixing ratios of gases.
CALL r2_set_gas_mix_ratio(control, spectrum, atm, &
! Grid
  n_profile, n_layer, nozone, nd_field, list, &
! Gas mass mixing ratios
  h2o, co2, o3, n2o_mix_ratio, ch4_mix_ratio, so2_mix_ratio, &
  cfc11_mix_ratio, cfc12_mix_ratio, o2_mix_ratio, &
  cfc113_mix_ratio, cfc114_mix_ratio, hcfc22_mix_ratio, &
  hfc125_mix_ratio, hfc134a_mix_ratio, &
! 3D CO2
  co2_dim1, co2_dim2, co2_3d, l_co2_3d, &
! Chemical greenhouse gas fields
  ngrgas, grgas_field)


IF (control%i_angular_integration == ip_spherical_harmonic) THEN

  IF (control%i_sph_mode == ip_sph_mode_flux) THEN
    atm%n_viewing_level = n_layer+1
    DO i=1, n_layer + 1
      atm%viewing_level(i) = REAL(i-1)
    END DO
    atm%n_direction = 1
  ELSE
    ! One viewing level at the top of the atmosphere.
    atm%n_viewing_level    = 1
    atm%viewing_level(1)   = 0.0
    atm%n_direction        = 1

    SELECT CASE(control%l_geostationary)
    CASE (.TRUE.)
      ! Define the Cartesian position of the satellite.
      r_sat(1) = ( planet_radius + control%sat_hgt ) * &
                   COS(control%sat_lat) * COS(control%sat_lon)
      r_sat(2) = ( planet_radius + control%sat_hgt ) * &
                   COS(control%sat_lat) * SIN(control%sat_lon)
      r_sat(3) = ( planet_radius + control%sat_hgt ) * &
                   SIN(control%sat_lat)
      DO l = 1, atm%n_profile
        ! Define point in satelite footprint the satelite is looking at:
        r_obs(1) = planet_radius * &
                     COS(atm%lat(l)) * COS(atm%lon(l))
        r_obs(2) = planet_radius * &
                     COS(atm%lat(l)) * SIN(atm%lon(l))
        r_obs(3) = planet_radius * &
                     SIN(atm%lat(l))

        r_diff(1) = r_sat(1) - r_obs(1)
        r_diff(2) = r_sat(2) - r_obs(2)
        r_diff(3) = r_sat(3) - r_obs(3)
        mag_rdiff = SQRT( r_diff(1) * r_diff(1) + &
                          r_diff(2) * r_diff(2) + &
                          r_diff(3) * r_diff(3) )
        recip_mag_robs  = 1.0 / SQRT( r_obs(1) * r_obs(1) + &
                                      r_obs(2) * r_obs(2) + &
                                      r_obs(3) * r_obs(3) )
        atm%direction(l, 1, 1) = ( r_obs(1) * r_diff(1) + &
                                    r_obs(2) * r_diff(2) + &
                                    r_obs(3) * r_diff(3) ) / &
                                  ( planet_radius * mag_rdiff )

        ! Calculate Azimuthal angle
        IF (control%isolir == ip_solar) THEN
          ! Define the normalised vector pointing to the sun
          h_angle = 2.0*pi*(0.5-seconds_since_midnight/rsec_per_day)
          r_sun(1)=COS(sindec)*COS(h_angle)
          r_sun(2)=COS(sindec)*SIN(h_angle)
          r_sun(3)=SIN(sindec)

          ! Projection of r_diff on local horizontal plane
          r_hor(1)=r_diff(1)-atm%direction(l,1,1) &
             *r_obs(1)*mag_rdiff*recip_mag_robs
          r_hor(2)=r_diff(2)-atm%direction(l,1,1) &
             *r_obs(2)*mag_rdiff*recip_mag_robs
          r_hor(3)=r_diff(3)-atm%direction(l,1,1) &
             *r_obs(3)*mag_rdiff*recip_mag_robs

          ! Projection of unity vector r_sun on horizontal plane
          solar_zen=(r_obs(1)*r_sun(1)+r_obs(2)*r_sun(2)  &
                   +r_obs(3)*r_sun(3))*recip_mag_robs

          ! 1.0 refers to magnitude of unit vector pointing to the sun.
          r_hor_s(1)=r_sun(1)-solar_zen*r_obs(1)*recip_mag_robs
          r_hor_s(2)=r_sun(2)-solar_zen*r_obs(2)*recip_mag_robs
          r_hor_s(3)=r_sun(3)-solar_zen*r_obs(3)*recip_mag_robs

          ! Calculate the azimuth angle and store in atm%direction(l,1,2)
          mag_rhor  = SQRT( r_hor(1) * r_hor(1) + &
                            r_hor(2) * r_hor(2) + &
                            r_hor(3) * r_hor(3) )
          mag_rhor_s  = SQRT( r_hor_s(1) * r_hor_s(1) + &
                              r_hor_s(2) * r_hor_s(2) + &
                              r_hor_s(3) * r_hor_s(3) )
          atm%direction(l,1,2)= - (r_hor(1)*r_hor_s(1) &
                                  + r_hor(2)*r_hor_s(2) &
                                  + r_hor(3)*r_hor_s(3)) &
                                  /mag_rhor/mag_rhor_s
          atm%direction(l,1,2)=ACOS(atm%direction(l,1,2))
        ELSE
          ! For infra-red simulations we set the azimuth to 0.0.
          atm%direction(l, 1, 2) = 0.0
        END IF
      END DO
    CASE (.FALSE.)
      cmessage = 'Only geostationary satellites are available so far'
      ierr=i_err_fatal
      CALL ereport(RoutineName, ierr, cmessage)
    END SELECT

  END IF

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_atm

END MODULE set_atm_mod
