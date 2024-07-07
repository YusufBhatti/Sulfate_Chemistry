! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!-----------------------------------------------------------------------

MODULE fill_missing_data_sw_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FILL_MISSING_DATA_SW_MOD'
CONTAINS

SUBROUTINE fill_missing_data_sw(                                  &
  off_x, off_y, row_length, rows, model_levels, ntiles,           &
  salt_dim1, salt_dim2, salt_dim3,                                &
  cloud_levels,first_row,last_row,                                &
  first_data_interp, es_space_interp,                             &
  l_complete_north, l_complete_south, l_complete_deg,             &
  l_flux_below_690nm_surf, n_channel, j_sw, l_extra_top,          &
  sw_incs, netsw, swsea,flux_below_690nm_surf,                    &
  top_absorption, surf_down_sw, sea_salt_film, sea_salt_jet)

!  Subroutine to fill in missing data

! Method:

!   Radiative fluxes may not have been calculated at all
!   points: we now fill in as required. This part of the
!   code was originally located in RAD_CTL2 (v6.1 and below)
!   but has been move into a subroutine in order to make
!   RAD_CTL2 more readable.


USE spec_sw_lw, ONLY: sw_spectrum
USE solinc_data, ONLY:                                            &
  f_orog, l_orog
USE sw_diag_mod, ONLY:                                            &
  sw_diag

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE rad3d_inp_mod, ONLY: rad3d_inp
IMPLICIT NONE


! VARIABLES WITH INTENT IN

INTEGER ::                                                        &
  off_x                                                           &
, off_y                                                           &
, row_length                                                      &
, rows                                                            &
, model_levels                                                    &
, ntiles                                                          &
, cloud_levels                                                    &
, salt_dim1                                                       &
, salt_dim2                                                       &
, salt_dim3

INTEGER ::                                                        &
  first_row                                                       &
, last_row                                                        &
, first_data_interp                                               &
, n_channel                                                       &
, j_sw

REAL ::                                                           &
   es_space_interp(4, row_length, rows)

LOGICAL ::                                                        &
  l_complete_north                                                &
, l_complete_south                                                &
, l_complete_deg                                                  &
, l_flux_below_690nm_surf                                         &
, l_extra_top

! VARIABLES WITH INTENT IN/OUT

REAL ::                                                           &
  sw_incs(row_length, rows, 0:model_levels+1)                     &
, netsw(row_length, rows)                                         &
, swsea(row_length, rows)                                         &
, surf_down_sw(row_length,rows,4)                                 &
, top_absorption(row_length, rows)                                &
, flux_below_690nm_surf(row_length,rows)                          &
, sea_salt_film(salt_dim1, salt_dim2, salt_dim3)                  &
, sea_salt_jet(salt_dim1, salt_dim2, salt_dim3)

INTEGER :: iwv
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_MISSING_DATA_SW'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!           Primary Fields:

CALL rad3d_inp(                                                 &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     model_levels+2,                                            &
     sw_incs                                                    &
     )
CALL rad3d_inp(                                                 &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     netsw                                                      &
     )
CALL rad3d_inp(                                                 &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     swsea                                                      &
     )

! The next field is purely diagnostic in some versions of
! the model, but is always required with MOSES.

IF ( l_flux_below_690nm_surf ) THEN

  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     flux_below_690nm_surf                                      &
     )
END IF

IF (l_extra_top) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     top_absorption                                             &
     )
END IF

CALL rad3d_inp(                                                 &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     4,                                                         &
     surf_down_sw                                               &
     )

IF (l_orog) THEN
  CALL rad3d_inp(                                               &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, 1,                    &
      f_orog                                                    &
      )
END IF

! Complete the diagnostic fields as required.

IF ( sw_diag(j_sw)%l_rad_mask ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%rad_mask )
END IF
IF ( sw_diag(j_sw)%l_flux_up ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_up )
END IF
IF ( sw_diag(j_sw)%l_flux_down ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_down )
END IF
IF ( sw_diag(j_sw)%l_flux_up_clear ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_up_clear )
END IF
IF ( sw_diag(j_sw)%l_flux_down_clear ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_down_clear )
END IF
IF ( sw_diag(j_sw)%l_solar_out_toa ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%solar_out_toa                                &
     )
END IF
IF ( sw_diag(j_sw)%l_solar_out_clear ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%solar_out_clear                              &
     )
END IF
IF ( sw_diag(j_sw)%l_surface_down_flux ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%surface_down_flux                            &
     )
END IF
IF ( sw_diag(j_sw)%l_surf_down_clr ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%surf_down_clr                                &
     )
END IF
IF ( sw_diag(j_sw)%l_surf_up_clr ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%surf_up_clr                                  &
     )
END IF
IF ( sw_diag(j_sw)%l_clear_hr ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     model_levels,                                              &
     sw_diag(j_sw)%clear_hr                                     &
     )
END IF
IF ( sw_diag(j_sw)%l_net_flux_trop ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%net_flux_trop                                &
     )
END IF
IF ( sw_diag(j_sw)%l_up_flux_trop ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%up_flux_trop                                 &
     )
END IF

!           Microphysical diagnostics

IF ( sw_diag(j_sw)%re_conv_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%re_conv                                      &
     )
END IF
IF ( sw_diag(j_sw)%re_strat_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%re_strat                                     &
     )
END IF
IF ( sw_diag(j_sw)%wgt_conv_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%wgt_conv                                     &
     )
END IF
IF ( sw_diag(j_sw)%wgt_strat_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%wgt_strat                                    &
     )
END IF
IF ( sw_diag(j_sw)%lwp_strat_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%lwp_strat                                    &
     )
END IF
IF ( sw_diag(j_sw)%weighted_re_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%weighted_re                                  &
     )
END IF
IF ( sw_diag(j_sw)%sum_weight_re_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%sum_weight_re                                &
     )
END IF
IF ( sw_diag(j_sw)%wgtd_warm_re_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%weighted_warm_re                             &
     )
END IF
IF ( sw_diag(j_sw)%sum_wgt_warm_re_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%sum_weight_warm_re                           &
     )
END IF
IF ( sw_diag(j_sw)%cdnc_ct_diag_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%cdnc_ct_diag                                 &
     )
END IF
IF ( sw_diag(j_sw)%cdnc_ct_weight_flag) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%cdnc_ct_weight                               &
     )
END IF
IF ( sw_diag(j_sw)%nc_diag_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%nc_diag                                      &
     )
END IF
IF ( sw_diag(j_sw)%nc_weight_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%nc_weight                                    &
     )
END IF
IF ( sw_diag(j_sw)%ntot_diag_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%ntot_diag                                    &
     )
END IF
IF ( sw_diag(j_sw)%strat_lwc_diag_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%strat_lwc_diag                               &
     )
END IF
IF ( sw_diag(j_sw)%so4_ccn_diag_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%so4_ccn_diag                                 &
     )
END IF
IF ( sw_diag(j_sw)%cond_samp_wgt_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     cloud_levels,                                              &
     sw_diag(j_sw)%cond_samp_wgt                                &
     )
END IF
IF ( sw_diag(j_sw)%seasalt_film_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     salt_dim3,                                                 &
     sea_salt_film                                              &
     )
END IF
IF ( sw_diag(j_sw)%seasalt_jet_flag ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp,                        &
     salt_dim3,                                                 &
     sea_salt_jet                                               &
     )
END IF

! Diagnostics for MOSES

IF ( sw_diag(j_sw)%l_FlxSolBelow690nmSurf ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%FlxSolBelow690nmSurf                         &
     )
END IF
IF ( sw_diag(j_sw)%l_FlxSeaBelow690nmSurf ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%FlxSeaBelow690nmSurf                         &
     )
END IF

! Diagnostic for Radiance

IF ( sw_diag(j_sw)%l_toa_radiance ) THEN
  CALL rad3d_inp(                                           &
    l_complete_north, l_complete_south, l_complete_deg,     &
    row_length, rows, off_x, off_y, first_row, last_row,    &
    first_data_interp, es_space_interp, n_channel,          &
    sw_diag(j_sw)%toa_radiance                              &
    )
END IF

! Diagnostic for Fluxes

IF ( sw_diag(j_sw)%l_uvflux_direct ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1,                                           &
     sw_diag(j_sw)%uvflux_direct                               &
    )
END IF
IF ( sw_diag(j_sw)%l_uvflux_up ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1,                                           &
     sw_diag(j_sw)%uvflux_up                                   &
    )
END IF
IF ( sw_diag(j_sw)%l_uvflux_down ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1,                                           &
     sw_diag(j_sw)%uvflux_down                                 &
    )
END IF
IF ( sw_diag(j_sw)%l_surf_uv ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%surf_uv                                      &
     )
END IF
IF ( sw_diag(j_sw)%l_surf_uv_clr ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%surf_uv_clr                                  &
     )
END IF
IF ( sw_diag(j_sw)%l_flux_direct ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1,                                           &
     sw_diag(j_sw)%flux_direct                                 &
    )
END IF
IF ( sw_diag(j_sw)%l_flux_diffuse ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1,                                           &
     sw_diag(j_sw)%flux_diffuse                                &
    )
END IF

! Diagnostics for orography correction

IF ( sw_diag(j_sw)%l_orog_corr ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, 1,                     &
     sw_diag(j_sw)%orog_corr                                    &
     )
END IF

! Extinction diagnostics:

IF ( sw_diag(j_sw)%l_cloud_extinction ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, cloud_levels,          &
     sw_diag(j_sw)%cloud_extinction                             &
     )
END IF
IF ( sw_diag(j_sw)%l_cloud_weight_extinction ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, cloud_levels,          &
     sw_diag(j_sw)%cloud_weight_extinction                      &
     )
END IF
IF ( sw_diag(j_sw)%l_ls_cloud_extinction ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, cloud_levels,          &
     sw_diag(j_sw)%ls_cloud_extinction                          &
     )
END IF
IF ( sw_diag(j_sw)%l_ls_cloud_weight_extinction ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, cloud_levels,          &
     sw_diag(j_sw)%ls_cloud_weight_extinction                   &
     )
END IF
IF ( sw_diag(j_sw)%l_cnv_cloud_extinction ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, cloud_levels,          &
     sw_diag(j_sw)%cnv_cloud_extinction                         &
     )
END IF
IF ( sw_diag(j_sw)%l_cnv_cloud_weight_extinction ) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, cloud_levels,          &
     sw_diag(j_sw)%cnv_cloud_weight_extinction                  &
     )
END IF

! Albedo diagnostics
IF ( sw_diag(j_sw)%l_direct_albedo) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp,es_space_interp,                         &
     sw_spectrum(j_sw)%basic%n_band,                            &
     sw_diag(j_sw)%direct_albedo                                &
     )
END IF
IF ( sw_diag(j_sw)%l_diffuse_albedo) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp,es_space_interp,                         &
     sw_spectrum(j_sw)%basic%n_band,                            &
     sw_diag(j_sw)%diffuse_albedo                               &
     )
END IF
IF ( sw_diag(j_sw)%l_vis_albedo_sc) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, ntiles,                &
     sw_diag(j_sw)%vis_albedo_sc                                &
     )
END IF
IF ( sw_diag(j_sw)%l_nir_albedo_sc) THEN
  CALL rad3d_inp(                                               &
     l_complete_north, l_complete_south, l_complete_deg,        &
     row_length, rows, off_x, off_y, first_row, last_row,       &
     first_data_interp, es_space_interp, ntiles,                &
     sw_diag(j_sw)%nir_albedo_sc                                &
     )
END IF

! Aerosol optical properties
IF ( sw_diag(j_sw)%l_aerosol_optical_depth ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels,         &
      sw_diag(j_sw)%aerosol_optical_depth(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_aerosol_scat_optical_depth ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels,         &
      sw_diag(j_sw)%aerosol_scat_optical_depth(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_aerosol_asymmetry_scat ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels,         &
      sw_diag(j_sw)%aerosol_asymmetry_scat(:,:,:,iwv) )
  END DO
END IF

! Band-by-band Fluxes
IF ( sw_diag(j_sw)%l_flux_up_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_up_band(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_flux_down_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_down_band(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_flux_up_clear_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_up_clear_band(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_flux_down_clear_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_down_clear_band(:,:,:,iwv) )
  END DO
END IF

! Clean-air Fluxes
IF ( sw_diag(j_sw)%l_flux_up_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_up_clean )
END IF
IF ( sw_diag(j_sw)%l_flux_down_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_down_clean )
END IF
IF ( sw_diag(j_sw)%l_flux_up_clear_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_up_clear_clean )
END IF
IF ( sw_diag(j_sw)%l_flux_down_clear_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, sw_diag(j_sw)%flux_down_clear_clean )
END IF

! GHG forcing Fluxes
IF ( sw_diag(j_sw)%l_flux_up_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     sw_diag(j_sw)%flux_up_forc )
END IF
IF ( sw_diag(j_sw)%l_flux_down_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     sw_diag(j_sw)%flux_down_forc )
END IF
IF ( sw_diag(j_sw)%l_flux_up_clear_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     sw_diag(j_sw)%flux_up_clear_forc )
END IF
IF ( sw_diag(j_sw)%l_flux_down_clear_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     sw_diag(j_sw)%flux_down_clear_forc )
END IF
IF ( sw_diag(j_sw)%l_flux_up_forc_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_up_forc_band(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_flux_down_forc_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_down_forc_band(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_flux_up_clear_forc_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_up_clear_forc_band(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_flux_down_clear_forc_band ) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      sw_diag(j_sw)%flux_down_clear_forc_band(:,:,:,iwv) )
  END DO
END IF
IF ( sw_diag(j_sw)%l_easyaerosol_extinction) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      sw_diag(j_sw)%easyaerosol_extinction(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF
IF ( sw_diag(j_sw)%l_easyaerosol_absorption) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      sw_diag(j_sw)%easyaerosol_absorption(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF
IF (sw_diag(j_sw)%l_easyaerosol_scattering) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      sw_diag(j_sw)%easyaerosol_scattering(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF
IF (sw_diag(j_sw)%l_easyaerosol_asytimscat) THEN
  DO iwv = 1, sw_spectrum(j_sw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      sw_diag(j_sw)%easyaerosol_asytimscat(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fill_missing_data_sw
END MODULE fill_missing_data_sw_mod
