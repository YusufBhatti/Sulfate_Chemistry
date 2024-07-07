! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set the boundary fields (surface and top-of-atmosphere).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_bound_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SET_BOUND_MOD'
CONTAINS

SUBROUTINE set_bound(control, dimen, spectrum, bound, diag,                    &
    nd_field, n_points, n_layer, list, row_list, col_list,                     &
    land, flandg, ice_fraction,                                                &
    emis_land, t_rad_surf, t_rad_land, t_rad_sice, t_rad_sea,                  &
    open_sea_albedo, sea_ice_albedo, land_albedo,                              &
    coszin, lit, solar_constant, scs, cos_zen_sph, day_frac_sph)

USE rad_pcf
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_bound,    ONLY: StrBound, allocate_bound
USE def_diag,     ONLY: StrDiag
USE planet_constants_mod, ONLY: &
  l_planet_grey_surface, planet_emissivity, planet_albedo
USE jules_sea_seaice_mod, ONLY: l_ctile, emis_sea, emis_sice
USE solinc_data,  ONLY: orog_corr

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE


! Control options:
TYPE(StrCtrl),      INTENT(IN)  :: control

! Dimensions:
TYPE(StrDim),       INTENT(IN)  :: dimen

! Spectral data:
TYPE (StrSpecData), INTENT(IN)  :: spectrum

! Boundary properties:
TYPE(StrBound),     INTENT(OUT) :: bound

! Diagnostics:
TYPE (StrDiag),   INTENT(INOUT) :: diag

! Dimensions of arrays:
INTEGER, INTENT(IN) :: nd_field
!   Allocated size of fields of data

! Actual sizes used:
INTEGER, INTENT(IN) :: n_points
!   Number of atmospheric points
INTEGER, INTENT(IN) :: n_layer
!   Number of atmospheric layers for radiation calculations
INTEGER, INTENT(IN) :: list(nd_field)
!   List of points
INTEGER, INTENT(IN) :: row_list(nd_field)
!   List of row indices
INTEGER, INTENT(IN) :: col_list(nd_field)
!   List of column indices

LOGICAL, INTENT(IN) :: land(nd_field)
!   Land flag
REAL, INTENT(IN) :: flandg(nd_field)
!   Fraction of land in a grid-box
REAL, INTENT(IN) :: ice_fraction(nd_field)
!   Fraction of sea ice

REAL, INTENT(IN) :: emis_land(nd_field)
!   Mean land emissivity in a gridbox
REAL, INTENT(IN) :: t_rad_surf(nd_field)
!   Effective radiative temperature over whole grid-box
REAL, INTENT(IN) :: t_rad_land(nd_field)
!   Effective radiative temperature of land
REAL, INTENT(IN) :: t_rad_sice(nd_field)
!   Effective radiative temperature of sea-ice
REAL, INTENT(IN) :: t_rad_sea(nd_field)
!   Effective radiative temperature of sea

REAL, INTENT(IN) :: open_sea_albedo(nd_field, 2, spectrum%dim%nd_band)
!   Diffuse albedo field
REAL, INTENT(IN) :: sea_ice_albedo(nd_field, 4)
!   Sea-ice albedos
REAL, INTENT(IN) :: land_albedo(nd_field, 4)
!   Land surface albedo fields

REAL, INTENT(IN) :: coszin(nd_field)
!   Cosine of zenith angle
REAL, INTENT(IN) :: solar_constant
!   Total solar irradiance at 1 AU
REAL, INTENT(IN) :: scs
!   Scaling of solar incident field
REAL, INTENT(IN) :: lit(nd_field)
!   Fraction of time point is lit
REAL, INTENT(IN) :: cos_zen_sph(nd_field, 0:n_layer+1)
!   Cosines of zenith angle for each layer
REAL, INTENT(IN) :: day_frac_sph(nd_field, 0:n_layer+1)
!   Fraction of time point is lit within each layer

! Local variables.
INTEGER :: i, l, ll
!   Loop variables

! Arrays related to tiling of the surface
INTEGER :: list_tile_outer(dimen%nd_point_tile)
!   List of points with surface tiling in the full list
INTEGER :: index_tile(dimen%nd_tile_type)
!   The indexing number of tiles of the given type


CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SET_BOUND'
CHARACTER (LEN=errormessagelength) :: cmessage
INTEGER :: ierr

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate structure for the core radiation code interface
CALL allocate_bound(bound, dimen, spectrum)

! Set the radiative characteristics of the surface.
SELECT CASE (control%isolir)
CASE (ip_solar)
  ! Define weightings for the basis functions of the surface
  ! BRDFs: in effect these are surafce albedos. Each grid-box
  ! may contain land, sea or sea-ice, which are treated by tiles
  ! within the radiation scheme (which are used more generally
  ! then simply when the coastal tiling scheme is on).

  ! Note: Without coastal tiling the land albedo contains the
  ! albedo of the solid surface. If coastal tiling is enabled,
  ! the land albedo refers only to the land surface and there
  ! is a separate sea-ice albedo.

  ! The land_albedo array allows for a spectral dependence,
  ! split between the VIS and NIR, as well as direct/diffuse.

  IF (l_planet_grey_surface) THEN
    bound%rho_alb(1:n_points, ip_surf_alb_diff, 1:spectrum%basic%n_band)       &
      = planet_albedo
    bound%rho_alb(1:n_points, ip_surf_alb_dir,  1:spectrum%basic%n_band)       &
      = planet_albedo
  ELSE
    DO i=1, spectrum%basic%n_band
      DO l=1, n_points
        ! Oceanic surface.
        IF (flandg(list(l)) < 1.0) THEN
          bound%rho_alb(l, ip_surf_alb_diff, i)                                &
            =(spectrum%solar%weight_blue(i)*sea_ice_albedo(list(l), 2)         &
            +(1.0-spectrum%solar%weight_blue(i))*sea_ice_albedo(list(l), 4))   &
            *ice_fraction(list(l))                                             &
            + MAX(0.0, MIN(1.0, open_sea_albedo(list(l),2,i) ) ) *             &
                                      (1.0-ice_fraction(list(l)))
          bound%rho_alb(l, ip_surf_alb_dir, i)                                 &
            =(spectrum%solar%weight_blue(i)*sea_ice_albedo(list(l), 1)         &
            +(1.0-spectrum%solar%weight_blue(i))*sea_ice_albedo(list(l), 3))   &
            *ice_fraction(list(l))                                             &
            + MAX(0.0, MIN(1.0, open_sea_albedo(list(l),1,i) ) ) *             &
                                      (1.0-ice_fraction(list(l)))
        ELSE
          bound%rho_alb(l, ip_surf_alb_diff, i) = 0.0
          bound%rho_alb(l, ip_surf_alb_dir, i) = 0.0
        END IF
  
        ! Add contributions from the land.
        IF (flandg(list(l)) > 0.0) THEN
          IF (l_ctile) THEN
            bound%rho_alb(l, ip_surf_alb_diff, i)                              &
              =(1.0-flandg(list(l)))*bound%rho_alb(l, ip_surf_alb_diff, i)     &
              +flandg(list(l)) * (spectrum%solar%weight_blue(i)                &
              *land_albedo(list(l), 2)                                         &
              +(1.0-spectrum%solar%weight_blue(i))*land_albedo(list(l), 4))
            bound%rho_alb(l, ip_surf_alb_dir, i)                               &
              =(1.0-flandg(list(l)))*bound%rho_alb(l, ip_surf_alb_dir, i)      &
              +flandg(list(l)) * (spectrum%solar%weight_blue(i)                &
              *land_albedo(list(l), 1)                                         &
              +(1.0-spectrum%solar%weight_blue(i))*land_albedo(list(l), 3))
          ELSE
            bound%rho_alb(l, ip_surf_alb_diff, i)                              &
              =spectrum%solar%weight_blue(i)*land_albedo(list(l), 2)           &
              +(1.0-spectrum%solar%weight_blue(i))*land_albedo(list(l), 4)
            bound%rho_alb(l, ip_surf_alb_dir, i)                               &
              =spectrum%solar%weight_blue(i)*land_albedo(list(l), 1)           &
              +(1.0-spectrum%solar%weight_blue(i))*land_albedo(list(l), 3)
          END IF
        END IF
      END DO
    END DO
  END IF

  ! Check albedo is in physical limits:
  DO i=1, spectrum%basic%n_band
    DO l=1, n_points
      bound%rho_alb(l, ip_surf_alb_dir, i) = MAX(0.0, MIN(1.0,                 &
                             bound%rho_alb(l, ip_surf_alb_dir, i) ) )
      bound%rho_alb(l, ip_surf_alb_diff, i) = MAX(0.0, MIN(1.0,                &
                             bound%rho_alb(l, ip_surf_alb_diff, i) ) )
    END DO
  END DO

CASE (ip_infra_red)
  ! Surface temperature
  DO l=1, n_points
    bound%t_ground(l) = t_rad_surf(list(l))
  END DO

  ! Zero the irrelevant direct albedo.
  bound%rho_alb(1:n_points, ip_surf_alb_dir, 1:spectrum%basic%n_band) = 0.0

  ! Set the diffuse albedo
  IF (l_planet_grey_surface) THEN
    bound%rho_alb(1:n_points, ip_surf_alb_diff, 1:spectrum%basic%n_band)       &
      = 1.0 - planet_emissivity
  ELSE
    DO i=1,spectrum%basic%n_band
      DO l=1,n_points
        bound%rho_alb(l, ip_surf_alb_diff, i)                                  &
          = flandg(list(l)) * (1.0 - emis_land(list(l))) +                     &
            (1.0 - flandg(list(l))) * (                                        &
            (1.0 - ice_fraction(list(l))) * (1.0 - emis_sea) +                 &
            ice_fraction(list(l)) * (1.0 - emis_sice)                          &
            )
      END DO
    END DO
  END IF
END SELECT

! Set the surface basis functions for a Lambertian surface.
bound%n_brdf_basis_fnc=1
! By defining F_{1,0,0,0} to be 4, rho_alb becomes equal to the
! diffuse albedo.
bound%f_brdf(1, 0, 0, 0)=4.0
IF (control%ls_brdf_trunc /= 0) THEN
  cmessage = 'The order of surface truncation is too high.'
  ierr=i_err_fatal
  CALL ereport(RoutineName, ierr, cmessage)
END IF


IF (control%l_tile) THEN

  ! Set up the surface tiling variables. There are multiple levels
  ! of indexing. Over all points in the domain only those in the
  ! array list require radiative calculations and of these points
  ! only those in the array list_file require tiling, this array
  ! being indexed over points where radiative calculations are to
  ! be done. list_tile_outer is indexed over the tiled points and
  ! gives the index in the whole domain.

  IF (l_ctile) THEN

    ! With coastal tiling we can have land, open sea or sea ice
    ! in the grid-box.

    bound%n_tile=3
    index_tile(ip_ocean_tile)=1
    index_tile(ip_seaice_tile)=2
    index_tile(ip_land_tile)=3
    bound%n_point_tile=0
    DO ll=1, n_points
      l=list(ll)
      IF ( (flandg(l) < 1.0) .AND.                                             &
           ( (flandg(l) > 0.0) .OR.                                            &
             ( (ice_fraction(l) > 0.0) .AND.                                   &
               (ice_fraction(l) < 1.0) ) ) ) THEN
        bound%n_point_tile=bound%n_point_tile+1
        bound%list_tile(bound%n_point_tile)=ll
        list_tile_outer(bound%n_point_tile)=l
        ! Assign tiled fractions consistent with the indices above.
        bound%frac_tile(bound%n_point_tile, 1)                                 &
          =(1.0-flandg(l))*(1.0-ice_fraction(l))
        bound%frac_tile(bound%n_point_tile, 2)                                 &
          =(1.0-flandg(l))*ice_fraction(l)
        bound%frac_tile(bound%n_point_tile, 3)                                 &
          =flandg(l)
      END IF
    END DO

  ELSE

    ! Without coastal tiling we have only open sea or sea ice
    ! forming the coastal tiling.

    bound%n_tile=2
    index_tile(ip_ocean_tile)=1
    index_tile(ip_seaice_tile)=2
    bound%n_point_tile=0
    DO ll=1, n_points
      l=list(ll)
      IF ( (.NOT. land(l)) .AND.                                               &
           (ice_fraction(l) >  0.0) .AND.                                      &
           (ice_fraction(l) <  1.0) ) THEN
        bound%n_point_tile=bound%n_point_tile+1
        bound%list_tile(bound%n_point_tile)=ll
        list_tile_outer(bound%n_point_tile)=l
        bound%frac_tile(bound%n_point_tile, 1)                                 &
          =(1.0-ice_fraction(l))
        bound%frac_tile(bound%n_point_tile, 2)                                 &
          =ice_fraction(l)
      END IF
    END DO

  END IF


  ! Now assign the tiled surface properties at points where tiling is active.
  SELECT CASE (control%isolir)
  CASE (ip_solar)
    DO i=1, spectrum%basic%n_band
    
      ! The oceanic surface.
      bound%rho_alb_tile(1:bound%n_point_tile                                  &
        , ip_surf_alb_dir, ip_ocean_tile, i)                                   &
        =open_sea_albedo(list_tile_outer(1:bound%n_point_tile),1,i)
      bound%rho_alb_tile(1:bound%n_point_tile                                  &
        , ip_surf_alb_diff, ip_ocean_tile, i)                                  &
        =open_sea_albedo(list_tile_outer(1:bound%n_point_tile),2,i)
    
      IF (l_ctile) THEN
    
        ! With coastal tiling there is a real distinction between
        ! land and seaice.
    
        ! Seaice
        bound%rho_alb_tile(1:bound%n_point_tile                                &
          , ip_surf_alb_dir, ip_seaice_tile, i)                                &
          =spectrum%solar%weight_blue(i)                                       &
          *sea_ice_albedo(list_tile_outer(1:bound%n_point_tile), 1)            &
          +(1.0-spectrum%solar%weight_blue(i))                                 &
          *sea_ice_albedo(list_tile_outer(1:bound%n_point_tile), 3)
        bound%rho_alb_tile(1:bound%n_point_tile                                &
          , ip_surf_alb_diff, ip_seaice_tile, i)                               &
          =spectrum%solar%weight_blue(i)                                       &
          *sea_ice_albedo(list_tile_outer(1:bound%n_point_tile), 2)            &
          +(1.0-spectrum%solar%weight_blue(i))                                 &
          *sea_ice_albedo(list_tile_outer(1:bound%n_point_tile), 4)
    
        ! Land
        bound%rho_alb_tile(1:bound%n_point_tile                                &
          , ip_surf_alb_dir, ip_land_tile, i)                                  &
          =spectrum%solar%weight_blue(i)                                       &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 1)               &
          +(1.0-spectrum%solar%weight_blue(i))                                 &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 3)
        bound%rho_alb_tile(1:bound%n_point_tile                                &
          , ip_surf_alb_diff, ip_land_tile, i)                                 &
          =spectrum%solar%weight_blue(i)                                       &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 2)               &
          +(1.0-spectrum%solar%weight_blue(i))                                 &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 4)
      ELSE
    
        ! The land albedo fields contain the values for sea-ice.
        bound%rho_alb_tile(1:bound%n_point_tile                                &
          , ip_surf_alb_dir, ip_seaice_tile, i)                                &
          =spectrum%solar%weight_blue(i)                                       &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 1)               &
          +(1.0-spectrum%solar%weight_blue(i))                                 &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 3)
        bound%rho_alb_tile(1:bound%n_point_tile                                &
          , ip_surf_alb_diff, ip_seaice_tile, i)                               &
          =spectrum%solar%weight_blue(i)                                       &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 2)               &
          +(1.0-spectrum%solar%weight_blue(i))                                 &
          *land_albedo(list_tile_outer(1:bound%n_point_tile), 4)
    
      END IF
    
    END DO
    
    ! Check albedo is in physical limits:
    DO i=1, spectrum%basic%n_band
      DO l=1, bound%n_point_tile
        DO ll =1, bound%n_tile
          bound%rho_alb_tile(l, ip_surf_alb_dir, ll, i) = MAX(0.0,MIN(1.0,     &
                        bound%rho_alb_tile(l, ip_surf_alb_dir, ll, i) ) )
          bound%rho_alb_tile(l, ip_surf_alb_diff, ll, i) =MAX(0.0,MIN(1.0,     &
                        bound%rho_alb_tile(l, ip_surf_alb_diff, ll, i) ) )
        END DO
      END DO
    END DO

  CASE (ip_infra_red)
    ! Open sea and sea ice always need to be set with
    ! radiative tiling, but specific land fields are required only
    ! if coastal tiling is on.
    
    ! Zero the irrelevant direct albedos
    bound%rho_alb_tile(1:bound%n_point_tile, ip_surf_alb_dir                   &
      , 1:bound%n_tile, 1:spectrum%basic%n_band) = 0.0
    
    ! The oceanic surface.
    bound%rho_alb_tile(1:bound%n_point_tile, ip_surf_alb_diff                  &
      , ip_ocean_tile, 1:spectrum%basic%n_band) = (1.0 - emis_sea)
    
    ! Sea-ice.
    bound%rho_alb_tile(1:bound%n_point_tile, ip_surf_alb_diff                  &
      , ip_seaice_tile, 1:spectrum%basic%n_band) = (1.0 - emis_sice)
    
    ! Land points. The test on l_ctile is required only
    ! because without coastal tiling land points will be wholly land
    ! and do not require radiative tiling.
    IF ( l_ctile ) THEN
      DO i = 1, bound%n_point_tile
        bound%rho_alb_tile(i, ip_surf_alb_diff, ip_land_tile,                  &
          1:spectrum%basic%n_band) = 1.0 - emis_land(list_tile_outer(i))
      END DO
    END IF
    
    ! Tiled temperatures are required to handle emission from the
    ! surface. Ensure that the indexing of the second subscript is
    ! consistent with the assignment of index_tile above.
    bound%t_tile(1:bound%n_point_tile, 1)                                      &
      = t_rad_sea(list_tile_outer(1:bound%n_point_tile))
    IF (l_ctile) THEN
      bound%t_tile(1:bound%n_point_tile, 2)                                    &
        = t_rad_sice(list_tile_outer(1:bound%n_point_tile))
      bound%t_tile(1:bound%n_point_tile, 3)                                    &
        = t_rad_land(list_tile_outer(1:bound%n_point_tile))
    ELSE
      ! Without coastal tiling, the solid part of the grid-box
      ! can only be sea ice.
      bound%t_tile(1:bound%n_point_tile, 2)                                    &
        = t_rad_sice(list_tile_outer(1:bound%n_point_tile))
    END IF
  END SELECT
END IF


! Set the incident solar flux.
IF (control%isolir == ip_solar .OR. control%l_solar_tail_flux) THEN
  IF (control%l_spherical_solar) THEN
    DO i=0, n_layer+1
      DO l=1, n_points
        bound%cos_zen(l, i) = cos_zen_sph(list(l),n_layer+1-i)
        bound%lit(l, i) = day_frac_sph(list(l),n_layer+1-i)
      END DO
    END DO
  END IF
  DO l=1, n_points
    IF (control%l_spherical_solar) THEN
      bound%solar_irrad(l)=scs*solar_constant
    ELSE IF (coszin(list(l)) > 0.0) THEN
      bound%solar_irrad(l)=scs*solar_constant*lit(list(l))
      bound%zen_0(l)=1.0/coszin(list(l))
    ELSE
      bound%solar_irrad(l)=0.0
      bound%zen_0(l)=1.0
    END IF
    ! Gather the orography correction factor into lit points.
    IF (control%l_orog) THEN
      bound%orog_corr(l)=orog_corr(col_list(l),row_list(l))
    END IF
  END DO
  IF (control%i_angular_integration == ip_spherical_harmonic) THEN
    DO l=1, n_points
      bound%zen_0(l)=coszin(list(l))
    END DO
  END IF
END IF


! Write out diagnostics for the boundary fields
IF (diag%l_direct_albedo) THEN
  DO i=1, spectrum%basic%n_band
    DO l=1, n_points
      diag%direct_albedo(col_list(l), row_list(l), i)                          &
            = bound%rho_alb(l, ip_surf_alb_dir, i)
    END DO
  END DO
END IF
IF (diag%l_diffuse_albedo) THEN
  DO i=1, spectrum%basic%n_band
    DO l=1, n_points
      diag%diffuse_albedo(col_list(l), row_list(l), i)                         &
            = bound%rho_alb(l, ip_surf_alb_diff, i)
    END DO
  END DO
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE set_bound
END MODULE set_bound_mod
