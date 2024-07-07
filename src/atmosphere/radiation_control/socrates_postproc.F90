! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Post process the data after the core radiation calculations.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------

MODULE socrates_postproc_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'SOCRATES_POSTPROC_MOD'
CONTAINS

SUBROUTINE socrates_postproc(                                                  &
!                   Thermodynamic variables
    layer_heat_capacity,                                                       &
!                   Options for COSP
    l_cosp,                                                                    &
!                   Surface fields
    flandg, ice_fraction, land,                                                &
!                   Solar fields and grid-dependent arrays
    coszin, lit, day_frac_sph, list,                                           &
    obs_solid_angle, trans_solid_angle, dir_flux_to_trans,                     &
    trindx,                                                                    &
!                   Spectrum
    spectrum,                                                                  &
!                   Algorithmic options
    control, pts,                                                              &
!                   diagnostics
    sw_diag, row_list, col_list,                                               &
!                   Physical dimensions
    dimen, nlit, n_points, n_layer, nclds,                                     &
    row_length, rows,                                                          &
    nd_field, nd_field_flux_diag,                                              &
!                   Radiance core data
    atm, cld, aer, bound, radout, radout_clean, radout_forc,                   &
!                   Output
    surf_down_sw, flux_below_690nm_surf,                                       &
    netsw, top_absorption, swsea, swout,                                       &
!                   COSP input arguments
    cosp_gbx, cosp_sgx)

USE rad_pcf
USE conversions_mod, ONLY: pi
USE def_spectrum, ONLY: StrSpecData
USE def_dimen,    ONLY: StrDim
USE def_control,  ONLY: StrCtrl
USE def_atm,      ONLY: StrAtm,     deallocate_atm
USE def_cld,      ONLY: StrCld,     deallocate_cld,                            &
                                    deallocate_cld_prsc,                       &
                                    deallocate_cld_mcica
USE def_aer,      ONLY: StrAer,     deallocate_aer,                            &
                                    deallocate_aer_prsc
USE def_bound,    ONLY: StrBound,   deallocate_bound
USE def_out,      ONLY: StrOut,     deallocate_out
USE def_diag,     ONLY: StrDiag
USE solinc_data,  ONLY: f_orog
USE mcica_mod
USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE cosp_types_mod, ONLY: cosp_gridbox, cosp_subgrid
USE ereport_mod, ONLY: ereport
USE arcl_mod,    ONLY: npd_arcl_species, npd_arcl_compnts
USE rad_input_mod, ONLY: l_extra_top
USE nlsizes_namelist_mod, ONLY: model_levels
USE r2_calc_total_cloud_cover_mod, ONLY: r2_calc_total_cloud_cover
USE set_diag_mod, ONLY: set_diag

IMPLICIT NONE


! Dummy arguments

! Dimensions of arrays:
INTEGER, INTENT(IN) :: row_length
!                          length of rows on each domain
INTEGER, INTENT(IN) :: rows
!                          number of rows in the domain
INTEGER ::                                                                     &
          !, intent(in)
      nd_field                                                                 &
!       Field size in calling program
    , nd_field_flux_diag                                                       
!       Field size for flux diagnostics

! Actual sizes used:
INTEGER ::                                                                     &
          !, intent(in)
    n_points                                                                   &
!       Number of points to be diagnosed including unlit points
    , n_layer                                                                  &
!       number of layers seen in the radiation scheme
    , nclds                                                                    
!       Number of cloudy levels

REAL, INTENT(IN) :: layer_heat_capacity(nd_field, model_levels)
!   Heat capacity of layer

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

! Diagnostics:
TYPE (StrDiag), INTENT(INOUT) :: sw_diag

! Output fields from core radiation code:
TYPE(StrOut), INTENT(INOUT) :: radout

! Output fields from diagnostic calls:
TYPE(StrOut), INTENT(INOUT) :: radout_clean
TYPE(StrOut), INTENT(INOUT) :: radout_forc


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
    , lit(nd_field)                                                            &
!       Fraction of time point is lit
    , day_frac_sph(nd_field, 0:n_layer+1)
!       Fraction of time point is lit within each layer


! Flag for COSP
LOGICAL,INTENT(IN) :: l_cosp


! Properties of the surface:
LOGICAL, INTENT(IN) :: land(nd_field)
!                         Land mask

REAL ::                                                                        &
          !, intent(in)
    ice_fraction(nd_field)                                                     &
!         fraction of sea ice in sea portion of grid box

    , flandg(nd_field)                                                         
!         land fraction in grid box

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
!       Temperature increments
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


INTEGER, INTENT(IN) :: row_list(nd_field)
!                          list of row indices of lit points
INTEGER, INTENT(IN) :: col_list(nd_field)
!                          list of column indices of lit points

! Structure with COSP inputs
TYPE(cosp_gridbox),INTENT(INOUT) :: cosp_gbx
TYPE(cosp_subgrid),INTENT(INOUT) :: cosp_sgx


! Local variables.

! Array related to tiling of the surface
INTEGER ::                                                                     &
      index_tile(dimen%nd_tile_type)
!       The indexing number of tiles of the given type

INTEGER :: i, j, k, l, ll, lll, ic
!       Loop variables


REAL :: frac_solid
!   Total solid fraction in a grid-box

REAL :: total_cloud_cover_g(dimen%nd_profile)
!   Cloud fraction at gathered points

! Fluxes:
REAL ::                                                                        &
    flux_net(dimen%nd_flux_profile, 0: dimen%nd_layer, dimen%nd_channel),      &
!       Net flux
    flux_net_clear(dimen%nd_flux_profile, 0: dimen%nd_layer, dimen%nd_channel)
!       Clear-sky net total flux

! Small real number used in COSP diagnostics calculations
REAL, PARAMETER :: epsreal = EPSILON(1.0)

! Temporary variable used for COSP diagnostics when McICA is used
REAL :: cosp_temp

! Temporary field associated with the orography correction.
REAL ::                                                                        &
     swout_temp(nd_field)

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'SOCRATES_POSTPROC'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Output diagnostics common to SW and LW
CALL set_diag(nlit, n_layer, list, col_list, row_list, model_levels, &
  nd_field, obs_solid_angle, trans_solid_angle, dir_flux_to_trans, &
  control, atm, spectrum, radout, radout_clean, radout_forc, sw_diag)


! Prepare the output arrays:

! Processing depends on whether the code has been invoked to
! calculate radiances or fluxes.
IF ( (control%i_angular_integration == ip_two_stream) .OR.                     &
     ( (control%i_angular_integration == ip_spherical_harmonic) .AND.          &
        (control%i_sph_mode == ip_sph_mode_flux) ) ) THEN

  ! Convert downward fluxes to net fluxes.
  DO i=0, n_layer
    DO l=1, nlit
      flux_net(l, i, 1)=radout%flux_down(l, i, 1)-radout%flux_up(l, i, 1)
    END DO
  END DO
  IF (control%l_clear) THEN
    DO i=0, n_layer
      DO l=1, nlit
        flux_net_clear(l, i, 1)                                                &
          =radout%flux_down_clear(l, i, 1)-radout%flux_up_clear(l, i, 1)
      END DO
    END DO
  END IF


  ! Scatter the net downward flux at each level into SWOUT.
  DO i=1, model_levels+1
    DO l=1, nlit
      swout(list(l), i)=flux_net(l, n_layer+1-i, 1)
    END DO
  END DO
  DO l=1, nlit
    swout(list(l), model_levels+2)=swout(list(l), model_levels+1)
  END DO


  ! Net shortwave radiation absorbed by the planet
  ! (i. e. earth and atmosphere together):
  IF (control%l_spherical_solar) THEN
    DO l=1, nlit
      netsw(list(l)) = swout(list(l), model_levels+1)                          &
        + radout%flux_direct_sph(l, n_layer+1, 1)
    END DO
    DO i=model_levels, 1, -1
      DO l=1, nlit
        netsw(list(l)) = netsw(list(l))                                        &
          + radout%flux_direct_div(l, n_layer+1-i, 1)
      END DO
    END DO
  ELSE
    DO l=1, nlit
      netsw(list(l)) = swout(list(l), model_levels+1)
    END DO
  END IF

  IF (l_extra_top) THEN
    ! calculate the radiation absorbed in the extra layer
    ! above the top of the rest of the model.
    IF (control%l_spherical_solar) THEN
      DO l=1, nlit
        top_absorption(list(l))=flux_net(l, 0, 1)                              &
          - flux_net(l, n_layer-model_levels, 1)                               &
          + radout%flux_direct_div(l, 1, 1)
      END DO
    ELSE
      DO l=1, nlit
        top_absorption(list(l))=flux_net(l, 0, 1)                              &
          -flux_net(l, n_layer-model_levels, 1)
      END DO
    END IF
  END IF


  ! Extra direct SW flux reaching the surface due to the
  ! orography correction:
  IF (control%l_orog .AND. .NOT. control%l_spherical_solar) THEN
    DO l=1, nlit
      f_orog(col_list(l),row_list(l)) = radout%flux_direct(l, n_layer,1)       &
        * (bound%orog_corr(l) - 1.0)/bound%orog_corr(l)
    END DO
  END IF


  ! Assignment of diagnostics

  ! Note: purely diagnostic quantities allocated dynamically in
  ! RAD_CTL are zeroed there and need to be filled only at lit
  ! points.


  ! Outgoing solar radiation at TOA:
  IF (sw_diag%l_solar_out_toa) THEN
    DO l=1, nlit
      sw_diag%solar_out_toa(col_list(l), row_list(l))                          &
        = radout%flux_up(l, 0, 1)
    END DO
  END IF


  ! Clear-sky outgoing solar radiation at TOA:
  IF (sw_diag%l_solar_out_clear) THEN
    DO l=1, nlit
      sw_diag%solar_out_clear(col_list(l), row_list(l))                        &
        = radout%flux_up_clear(l, 0, 1)
    END DO
  END IF


  ! Total cloud cover:
  SELECT CASE (control%i_cloud)
  CASE (ip_cloud_mix_max, ip_cloud_mix_random, ip_cloud_triple,                &
        ip_cloud_part_corr, ip_cloud_part_corr_cnv)
    DO l=1, nlit
      total_cloud_cover_g(l)=radout%tot_cloud_cover(l)
    END DO
  CASE (ip_cloud_mcica)
    DO l=1,nlit
      i=((row_list(l)-1)*row_length)+col_list(l)
      total_cloud_cover_g(l)=REAL(ncldy(i))/REAL(tot_subcol_gen)
    END DO
  CASE DEFAULT
    CALL r2_calc_total_cloud_cover(nlit, nclds, nclds,                         &
      control%i_cloud, cld%w_cloud, total_cloud_cover_g,                       &
      dimen%nd_profile, dimen%nd_layer)
  END SELECT

  ! Total clear area
  IF (sw_diag%l_total_clear_area) THEN
    DO l=1, nlit
      sw_diag%total_clear_area(col_list(l), row_list(l))                       &
        =1.0-total_cloud_cover_g(l)
    END DO
  END IF

  ! Clear-sky flux at TOA weighted by total clear area
  IF (sw_diag%l_toa_clear_weighted) THEN
    DO l=1, nlit
      sw_diag%toa_clear_weighted(col_list(l), row_list(l))                     &
        =radout%flux_up_clear(l, 0, 1) * (1.0-total_cloud_cover_g(l))
    END DO
  END IF


  ! Surface flux below 690nm.
  IF (control%l_tile) THEN
    DO ll=1, nlit
      l=list(ll)
      IF (control%l_spherical_solar) THEN
        flux_below_690nm_surf(l) = radout%flux_direct_blue_surf(ll)            &
          + radout%flux_down_blue_surf(ll)-radout%flux_up_blue_surf(ll)
      ELSE
        flux_below_690nm_surf(l)                                               &
          = radout%flux_down_blue_surf(ll)-radout%flux_up_blue_surf(ll)
      END IF
      IF ( (flandg(l) < TINY(flandg)) .AND.                                    &
           (ice_fraction(l) < TINY(ice_fraction)) ) THEN
        ! This point is open sea with no sea ice.
        IF (sw_diag%l_FlxSeaBelow690nmSurf) THEN
          sw_diag%FlxSeaBelow690nmSurf(col_list(ll),row_list(ll))              &
            = flux_below_690nm_surf(l)
        END IF
      ELSE
        IF (sw_diag%l_FlxSolBelow690nmSurf) THEN
          sw_diag%FlxSolBelow690nmSurf(col_list(ll),row_list(ll))              &
            = flux_below_690nm_surf(l)
        END IF
      END IF
    END DO

    index_tile(ip_ocean_tile) = 1

    ! Tiled points will have both land and sea. Note that the
    ! channel index of flux_up_tile is hard-wired to 1 because
    ! we don't envisage calling the code in other cases.
    DO lll=1, bound%n_point_tile
      ll=bound%list_tile(lll)
      l=list(ll)
      IF (control%l_spherical_solar) THEN
        flux_below_690nm_surf(l)                                               &
          =(1.0-ice_fraction(l))*(radout%flux_down_blue_surf(ll)               &
          +radout%flux_direct_blue_surf(ll)                                    &
          -radout%flux_up_blue_tile(lll, index_tile(ip_ocean_tile), 1))
      ELSE
        flux_below_690nm_surf(l)                                               &
          =(1.0-ice_fraction(l))*(radout%flux_down_blue_surf(ll)               &
          -radout%flux_up_blue_tile(lll, index_tile(ip_ocean_tile), 1))
      END IF
      IF (sw_diag%l_FlxSeaBelow690nmSurf) THEN
        sw_diag%FlxSeaBelow690nmSurf(col_list(ll), row_list(ll))               &
          =flux_below_690nm_surf(l)
      END IF
      IF (sw_diag%l_FlxSolBelow690nmSurf) THEN
        frac_solid=flandg(l)+(1.0-flandg(l))*ice_fraction(l)
        IF (frac_solid > 0.0) THEN
          sw_diag%FlxSolBelow690nmSurf(col_list(ll),row_list(ll))              &
            =(sw_diag%FlxSolBelow690nmSurf(col_list(ll),row_list(ll))          &
            -(1.0-flandg(l))*flux_below_690nm_surf(l))/frac_solid
        END IF
      END IF
    END DO
  ELSE
    IF (control%l_spherical_solar) THEN
      DO ll=1, nlit
        l=list(ll)
        IF (land(l)) THEN
          flux_below_690nm_surf(l) = radout%flux_direct_blue_surf(ll)          &
            + radout%flux_down_blue_surf(ll) - radout%flux_up_blue_surf(ll)
        ELSE
          flux_below_690nm_surf(l) = (radout%flux_direct_blue_surf(ll)         &
            + radout%flux_down_blue_surf(ll) - radout%flux_up_blue_surf(ll) )  &
            * ( 1.0-ice_fraction(l) )
        END IF
      END DO
    ELSE
      DO ll=1, nlit
        l=list(ll)
        IF (land(l)) THEN
          flux_below_690nm_surf(l)                                             &
            =radout%flux_down_blue_surf(ll)-radout%flux_up_blue_surf(ll)
        ELSE
          flux_below_690nm_surf(l)                                             &
            =(radout%flux_down_blue_surf(ll)-radout%flux_up_blue_surf(ll))     &
            *(1.0-ice_fraction(l))
        END IF
      END DO
    END IF
  END IF


  ! Orography correction to direct SW flux:
  IF (sw_diag%l_orog_corr .AND. control%l_orog) THEN
    DO l=1, nlit
      sw_diag%orog_corr(col_list(l), row_list(l)) = bound%orog_corr(l)
    END DO
  END IF

  ! Components of downward flux at the surface:
  IF (control%l_spherical_solar) THEN
    surf_down_sw(list(1:nlit), 1)                                              &
      = radout%flux_direct_blue_surf(1:nlit)
    surf_down_sw(list(1:nlit), 2)                                              &
      = radout%flux_down_blue_surf(1:nlit)
    surf_down_sw(list(1:nlit), 3)                                              &
      = radout%flux_direct_sph(1:nlit, n_layer+1, 1)                           &
      - radout%flux_direct_blue_surf(1:nlit)
    surf_down_sw(list(1:nlit), 4)                                              &
      = radout%flux_down(1:nlit, n_layer, 1)                                   &
      - radout%flux_down_blue_surf(1:nlit)
  ELSE
    surf_down_sw(list(1:nlit), 1)                                              &
      =radout%flux_direct_blue_surf(1:nlit)
    surf_down_sw(list(1:nlit), 2)                                              &
      =radout%flux_down_blue_surf(1:nlit)-radout%flux_direct_blue_surf(1:nlit)
    surf_down_sw(list(1:nlit), 3)                                              &
      =radout%flux_direct(1:nlit, n_layer, 1)                                  &
      -radout%flux_direct_blue_surf(1:nlit)
    surf_down_sw(list(1:nlit), 4)                                              &
      =flux_net(1:nlit, n_layer, 1)+radout%flux_up(1:nlit, n_layer, 1)         &
      -radout%flux_down_blue_surf(1:nlit)                                      &
      -radout%flux_direct(1:nlit, n_layer, 1)                                  &
      +radout%flux_direct_blue_surf(1:nlit)
  END IF

  ! Downward flux at the surface:
  IF (sw_diag%l_surface_down_flux) THEN
    IF (control%l_spherical_solar) THEN
      DO l=1, nlit
        sw_diag%surface_down_flux(col_list(l), row_list(l))                    &
          = radout%flux_down(l, n_layer, 1)                                    &
          + radout%flux_direct_sph(l, n_layer+1, 1)
      END DO
    ELSE
      DO l=1, nlit
        sw_diag%surface_down_flux(col_list(l), row_list(l))                    &
          =flux_net(l, n_layer, 1)+radout%flux_up(l, n_layer, 1)
      END DO
    END IF
  END IF


  ! Clear-sky downward flux at the surface:
  IF (sw_diag%l_surf_down_clr) THEN
    IF (control%l_spherical_solar) THEN
      DO l=1, nlit
        sw_diag%surf_down_clr(col_list(l), row_list(l))                        &
          = radout%flux_down_clear(l, n_layer, 1)                              &
          + radout%flux_direct_clear_sph(l, n_layer+1, 1)
      END DO
    ELSE
      DO l=1, nlit
        sw_diag%surf_down_clr(col_list(l), row_list(l))                        &
          =flux_net_clear(l, n_layer, 1)                                       &
          +radout%flux_up_clear(l, n_layer, 1)
      END DO
    END IF
  END IF


  ! Clear-sky upward flux at the surface:
  IF (sw_diag%l_surf_up_clr) THEN
    DO l=1, nlit
      sw_diag%surf_up_clr(col_list(l), row_list(l))                            &
        =radout%flux_up_clear(l, n_layer, 1)
    END DO
  END IF


  ! Net flux at the tropopause:
  IF (sw_diag%l_net_flux_trop) THEN
    DO l=1, nlit
      sw_diag%net_flux_trop(col_list(l), row_list(l))                          &
        =flux_net(l, n_layer+1-trindx(list(l)), 1)
    END DO
  END IF


  ! Upward flux at the tropopause:
  IF (sw_diag%l_up_flux_trop) THEN
    DO l=1, nlit
      sw_diag%up_flux_trop(col_list(l), row_list(l))                           &
        =radout%flux_up(l, n_layer+1-trindx(list(l)), 1)
    END DO
  END IF

  ! Direct and diffuse downward flux
  IF (sw_diag%l_flux_direct) THEN
    DO i=1, model_levels+1
      DO l=1, nlit
        sw_diag%flux_direct(col_list(l), row_list(l),i)                        &
          = radout%flux_direct(l,n_layer+1-i,1)
      END DO
    END DO
  END IF
  IF (sw_diag%l_flux_diffuse) THEN
    DO i=1, model_levels+1
      DO l=1, nlit
        sw_diag%flux_diffuse(col_list(l), row_list(l),i)                       &
          = flux_net(l,n_layer+1-i,1) + radout%flux_up(l,n_layer+1-i,1)        &
          - radout%flux_direct(l,n_layer+1-i,1)
      END DO
    END DO
  END IF

  ! UV-Fluxes
  IF (sw_diag%l_uvflux_direct) THEN
    DO i=1, model_levels+1
      DO l=1, nlit
        sw_diag%uvflux_direct(col_list(l), row_list(l),i) = 0.0
      END DO
    END DO
    DO k=1, spectrum%basic%n_band
      DO i=1, model_levels+1
        DO l=1, nlit
          sw_diag%uvflux_direct(col_list(l), row_list(l),i)                    &
            = sw_diag%uvflux_direct(col_list(l), row_list(l),i)                &
            + control%weight_diag(k)*radout%flux_direct_band(l,n_layer+1-i,k)
        END DO
      END DO
    END DO
  END IF
  IF (sw_diag%l_uvflux_up) THEN
    DO i=1,model_levels+1
      DO l=1, nlit
        sw_diag%uvflux_up(col_list(l), row_list(l),i) = 0.0
      END DO
    END DO
    DO k=1, spectrum%basic%n_band
      DO i=1,model_levels+1
        DO l=1, nlit
          sw_diag%uvflux_up(col_list(l), row_list(l),i)                        &
            = sw_diag%uvflux_up(col_list(l), row_list(l),i)                    &
            + control%weight_diag(k)*radout%flux_up_band(l,n_layer+1-i,k)
        END DO
      END DO
    END DO
  END IF
  IF (sw_diag%l_uvflux_down) THEN
    DO i=1,model_levels+1
      DO l=1, nlit
        sw_diag%uvflux_down(col_list(l), row_list(l),i) = 0.0
      END DO
    END DO
    DO k=1, spectrum%basic%n_band
      DO i=1,model_levels+1
        DO l=1, nlit
          sw_diag%uvflux_down(col_list(l), row_list(l),i)                      &
            = sw_diag%uvflux_down(col_list(l), row_list(l),i)                  &
            + control%weight_diag(k)*radout%flux_down_band(l,n_layer+1-i,k)
        END DO
      END DO
    END DO
  END IF
  IF (sw_diag%l_surf_uv) THEN
    DO l=1, nlit
      sw_diag%surf_uv(col_list(l), row_list(l)) = 0.0
    END DO
    IF (control%l_spherical_solar) THEN
      DO k=1, spectrum%basic%n_band
        DO l=1, nlit
          sw_diag%surf_uv(col_list(l), row_list(l))                            &
            = sw_diag%surf_uv(col_list(l), row_list(l))                        &
            + control%weight_diag(k)*(radout%flux_down_band(l,n_layer,k)       &
            + radout%flux_direct_sph_band(l,n_layer+1,k))
        END DO
      END DO
    ELSE
      DO k=1, spectrum%basic%n_band
        DO l=1, nlit
          sw_diag%surf_uv(col_list(l), row_list(l))                            &
            = sw_diag%surf_uv(col_list(l), row_list(l))                        &
            + control%weight_diag(k)*radout%flux_down_band(l,n_layer,k)
        END DO
      END DO
    END IF
  END IF
  IF (sw_diag%l_surf_uv_clr) THEN
    DO l=1, nlit
      sw_diag%surf_uv_clr(col_list(l), row_list(l)) = 0.0
    END DO
    IF (control%l_spherical_solar) THEN
      DO k=1, spectrum%basic%n_band
        DO l=1, nlit
          sw_diag%surf_uv_clr(col_list(l), row_list(l))                        &
            = sw_diag%surf_uv_clr(col_list(l), row_list(l))                    &
            + control%weight_diag(k)*(radout%flux_down_clear_band(l,n_layer,k) &
            + radout%flux_direct_clear_sph_band(l,n_layer+1,k))
        END DO
      END DO
    ELSE
      DO k=1, spectrum%basic%n_band
        DO l=1, nlit
          sw_diag%surf_uv_clr(col_list(l), row_list(l))                        &
            = sw_diag%surf_uv_clr(col_list(l), row_list(l))                    &
            + control%weight_diag(k)*radout%flux_down_clear_band(l,n_layer,k)
        END DO
      END DO
    END IF
  END IF

  ! Cloud extinction diagnostics
  IF (sw_diag%l_cloud_extinction) THEN
    DO i=1, nclds
      DO l=1, nlit
        sw_diag%cloud_extinction                                               &
           (col_list(l), row_list(l), i)                                       &
           =radout%cloud_extinction(l, n_layer+1-i)
        sw_diag%cloud_weight_extinction                                        &
          (col_list(l), row_list(l), i)                                        &
          =radout%cloud_weight_extinction(l, n_layer+1-i)
      END DO
    END DO
  END IF

  IF (sw_diag%l_ls_cloud_extinction) THEN
    DO i=1, nclds
      DO l=1, nlit
        sw_diag%ls_cloud_extinction                                            &
           (col_list(l), row_list(l), i)                                       &
           =radout%ls_cloud_extinction(l, n_layer+1-i)
        sw_diag%ls_cloud_weight_extinction                                     &
           (col_list(l), row_list(l), i)                                       &
           =radout%ls_cloud_weight_extinction(l, n_layer+1-i)
      END DO
    END DO
  END IF

  IF (sw_diag%l_cnv_cloud_extinction) THEN
    DO i=1, nclds
      DO l=1, nlit
        sw_diag%cnv_cloud_extinction                                           &
           (col_list(l), row_list(l), i)                                       &
           =radout%cnv_cloud_extinction(l, n_layer+1-i)
        sw_diag%cnv_cloud_weight_extinction                                    &
           (col_list(l), row_list(l), i)                                       &
           =radout%cnv_cloud_weight_extinction(l, n_layer+1-i)
      END DO
    END DO
  END IF


  ! COSP arguments
  IF (l_cosp) THEN
    IF (control%i_cloud == IP_cloud_mcica .AND. &
      cosp_sgx%Ncolumns == tot_subcol_gen) THEN
      DO i=1, nclds
        DO l=1, nlit
          IF (radout%ls_cloud_weight_extinction(l,n_layer+1-i) > epsreal) THEN
            j = (row_list(l)-1)*row_length + col_list(l)
            ! Sub-grid cloud optical depth for COSP assumes
            ! same generated sub-columns for liquid and ice.
            cosp_temp = atm%mass(l,n_layer+1-i)*                               &
                radout%ls_cloud_extinction(l, n_layer+1-i) /                   &
                radout%ls_cloud_weight_extinction(l, n_layer+1-i)
            DO lll=1, ncldy(j)
              cosp_sgx%dtau(j,lll,i) =                                         &
                clw_sub_full(j,n_layer+1-i,lll)*cosp_temp
            END DO
          END IF
        END DO
      END DO
    ELSE
      DO i=1, nclds
        DO l=1, nlit
          j = (row_list(l)-1)*row_length + col_list(l)
          ! Large-scale cloud optical depth
          IF (radout%ls_cloud_weight_extinction(l,n_layer+1-i) > epsreal) THEN
            cosp_gbx%dtau_s(j,i) = atm%mass(l,n_layer+1-i)                     &
               * radout%ls_cloud_extinction(l, n_layer+1-i)                    &
               / radout%ls_cloud_weight_extinction(l, n_layer+1-i)
          END IF
          ! Convective cloud optical depth
          IF (radout%cnv_cloud_weight_extinction(l,n_layer+1-i) > epsreal) THEN
            cosp_gbx%dtau_c(j,i) = atm%mass(l,n_layer+1-i)                     &
               * radout%cnv_cloud_extinction(l, n_layer+1-i)                   &
               / radout%cnv_cloud_weight_extinction(l, n_layer+1-i)
          END IF
        END DO
      END DO
    END IF
  END IF





  ! Final processing of output fields

  IF (control%l_spherical_solar) THEN

    ! Convert the fluxes to increments.
    DO i=model_levels, 1, -1
      DO l=1, nlit
        swout(list(l),i+1) = ( swout(list(l),i+1) - swout(list(l),i)           &
          + radout%flux_direct_div(l, n_layer+1-i, 1) )                        &
          * pts / layer_heat_capacity(list(l), i)
      END DO
    END DO
    DO l=1, nlit
      swout(list(l), 1) = swout(list(l), 1)                                    &
        + radout%flux_direct_sph(l, n_layer+1, 1)
    END DO
    IF (sw_diag%l_clear_hr) THEN
      DO i=model_levels, 1, -1
        DO l=1, nlit
          sw_diag%clear_hr(col_list(l), row_list(l), i) =                      &
            ( flux_net_clear(l, n_layer-i,   1)                                &
            - flux_net_clear(l, n_layer+1-i, 1)                                &
            + radout%flux_direct_clear_div(l, n_layer+1-i, 1) )                &
            / layer_heat_capacity(list(l), i)
        END DO
      END DO
    END IF

  ELSE

    IF (control%l_orog) THEN
      DO l = 1, nlit
        swout_temp(list(l)) = swout(list(l),1)
        swout(list(l),1)    = swout(list(l),1) - f_orog(col_list(l),row_list(l))
      END DO
    END IF
    ! Convert the fluxes to increments.
    DO i=model_levels, 1, -1
      DO l=1, nlit
        swout(list(l),i+1)=(swout(list(l),i+1)-swout(list(l),i))               &
          *pts/layer_heat_capacity(list(l), i)
      END DO
      IF (sw_diag%l_clear_hr) THEN
        DO l=1, nlit
          sw_diag%clear_hr(col_list(l), row_list(l), i)                        &
            =(flux_net_clear(l, n_layer-i, 1)                                  &
            -flux_net_clear(l, n_layer+1-i, 1))                                &
            /layer_heat_capacity(list(l), i)
        END DO
      END IF
    END DO
    IF (control%l_orog) THEN
      DO l = 1, nlit
        swout(list(l),1)=swout_temp(list(l))
      END DO
    END IF

  END IF



  ! Separate the contributions over open sea and sea ice.
  ! Fluxes returned from the radiation code itself are not
  ! weighted by the fraction of the tile, but here are converted
  ! to grid-box mean values. This can only be done if the surface
  ! has been tiled.

  IF (control%l_tile) THEN

    ! The variable flandg is set even if coastal tiling is not
    ! used, so fairly generic code can be written. Note that
    ! swsea is initialized in the calling routine, so its
    ! setting at land points is implicit.

    DO ll=1, nlit
      l=list(ll)
      IF ( (flandg(l) < TINY(flandg)) .AND.                                    &
           (ice_fraction(l) < TINY(ice_fraction) )                             &
         ) THEN
        ! This point is open sea with no sea ice.
        swsea(l)=swout(l, 1)
        IF (.NOT. control%l_spherical_solar) swout(l, 1)=0.0
      END IF
    END DO

    index_tile(ip_ocean_tile) = 1

    ! Tiled points will have both land and sea. Note that the
    ! channel index of flux_up_tile is hard-wired to 1 because
    ! we don't envisage calling the code in other cases.
    DO lll=1, bound%n_point_tile
      ll=bound%list_tile(lll)
      l=list(ll)
      swsea(l)=(1.0-ice_fraction(l))*(swout(l, 1)                              &
        +radout%flux_up(ll, n_layer, 1)                                        &
        -radout%flux_up_tile(lll, index_tile(ip_ocean_tile), 1))
      IF (.NOT. control%l_spherical_solar) THEN
        swout(l, 1)=swout(l, 1)-(1.0-flandg(l))*swsea(l)
      END IF
    END DO

    ! The remaining points are entirely land points and swout
    ! need not be altered.

  ELSE

    ! Without radiative tiling we must assume that fluxes are
    ! uniform across the grid-box.
    IF (control%l_spherical_solar) THEN
      WHERE (flandg(list(1:nlit)) < 1.0-TINY(flandg))
        swsea(list(1:nlit))=(1.0-ice_fraction(list(1:nlit)))                   &
          *swout(list(1:nlit), 1)
      END WHERE
    ELSE
      WHERE (flandg(list(1:nlit)) < 1.0-TINY(flandg))
        swsea(list(1:nlit))=(1.0-ice_fraction(list(1:nlit)))                   &
          *swout(list(1:nlit), 1)
        swout(list(1:nlit), 1)=swout(list(1:nlit), 1)                          &
          -(1.0-flandg(list(1:nlit)))*swsea(list(1:nlit))
      END WHERE
    END IF

  END IF


  IF (.NOT. control%l_spherical_solar) THEN
    ! Divide by cosine of solar zenith angle to provide values for
    ! upper routines. This applies only to SWOUT.
    DO i=1, model_levels+2
      DO l=1, n_points
        swout(l, i)=swout(l, i)/(coszin(l)*lit(l)+TINY(1.0))
      END DO
    END DO
    
    ! Process surface fluxes
    DO i=1, 4
      surf_down_sw(1:n_points, i) = surf_down_sw(1:n_points, i)                &
        / (coszin(1:n_points) * lit(1:n_points) + TINY(coszin) )
    END DO
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


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE socrates_postproc
END MODULE socrates_postproc_mod
