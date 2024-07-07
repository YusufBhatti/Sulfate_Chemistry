! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!-----------------------------------------------------------------------
!
MODULE fill_missing_data_lw_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FILL_MISSING_DATA_LW_MOD'
CONTAINS

SUBROUTINE fill_missing_data_lw(                                  &
  off_x, off_y, row_length, rows, model_levels,                   &
  cloud_levels, n_aod_wavel,                                      &
  first_row,last_row,                                             &
  first_data_interp, ES_space_interp,                             &
  L_complete_North, L_complete_South, L_complete_deg,             &
  n_channel, j_lw, l_extra_top,                                   &
  LW_incs, olr, lw_down, LWsea, top_absorption )

!  Subroutine to fill in missing data
!
! Method:
!
!   Radiative fluxes may not have been calculated at all
!   points: we now fill in as required. This part of the
!   code was originally located in RAD_CTL2 (v6.1 and below)
!   but has been move into a subroutine in order to make
!   RAD_CTL2 more readable.


USE spec_sw_lw, ONLY: lw_spectrum
USE lw_diag_mod, ONLY: LW_diag
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE rad3d_inp_mod, ONLY: rad3d_inp
IMPLICIT NONE

!
! VARIABLES WITH INTENT IN
!
INTEGER ::                                                        &
  off_x                                                           &
, off_y                                                           &
, row_length                                                      &
, rows                                                            &
, model_levels                                                    &
, cloud_levels                                                    &
, n_aod_wavel

INTEGER ::                                                        &
  first_row                                                       &
, last_row                                                        &
, first_data_interp                                               &
, n_channel                                                       &
, j_lw

REAL ::                                                           &
   es_space_interp(4, row_length, rows)

LOGICAL ::                                                        &
  L_complete_North                                                &
, L_complete_South                                                &
, L_complete_deg                                                  &
, l_extra_top
!
! VARIABLES WITH INTENT IN/OUT
!
REAL ::                                                           &
  LW_incs(row_length, rows, 0:model_levels)                       &
, LWsea(row_length, rows)                                         &
, olr(row_length, rows)                                           &
, lw_down(row_length, rows)                                       &
, top_absorption(row_length, rows)

INTEGER :: iwv
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FILL_MISSING_DATA_LW'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!
! Primary Fields:
!
CALL rad3d_inp(                                                 &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp,                         &
    model_levels+1,                                             &
    LW_incs                                                     &
    )
CALL rad3d_inp(                                                 &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    olr                                                         &
    )
CALL rad3d_inp(                                                 &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    lw_down                                                     &
    )
CALL rad3d_inp(                                                 &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    LWsea                                                       &
    )
IF (l_extra_top) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    top_absorption                                              &
    )
END IF
!
! LW Diagnostics:
!
IF ( lw_diag(j_lw)%l_flux_up ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_up )
END IF
IF ( lw_diag(j_lw)%l_flux_down ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_down )
END IF
IF ( lw_diag(j_lw)%l_flux_up_clear ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_up_clear )
END IF
IF ( lw_diag(j_lw)%l_flux_down_clear ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_down_clear )
END IF
IF ( LW_diag(j_lw)%L_total_cloud_cover ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    LW_diag(j_lw)%total_cloud_cover                             &
    )
END IF
IF ( LW_diag(j_lw)%L_clear_olr ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    LW_diag(j_lw)%clear_olr                                     &
    )
END IF
IF ( LW_diag(j_lw)%L_surf_down_clr ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    LW_diag(j_lw)%surf_down_clr                                 &
    )
END IF
IF ( LW_diag(j_lw)%L_clear_hr ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp,                         &
    model_levels,                                               &
    LW_diag(j_lw)%clear_hr                                      &
    )
END IF
IF ( LW_diag(j_lw)%L_net_flux_trop ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    LW_diag(j_lw)%net_flux_trop                                 &
    )
END IF
IF ( LW_diag(j_lw)%L_down_flux_trop ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, 1,                      &
    LW_diag(j_lw)%down_flux_trop                                &
    )
END IF
IF ( LW_diag(j_lw)%L_total_cloud_on_levels ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, cloud_levels,           &
    LW_diag(j_lw)%total_cloud_on_levels                         &
    )
END IF
!
! Grid-box mean cloud diagnostics as seen by radiation:
!
IF ( LW_diag(j_lw)%L_ls_qcl_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ls_qcl_rad                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_ls_qcf_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ls_qcf_rad                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_cc_qcl_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%cc_qcl_rad                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_cc_qcf_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%cc_qcf_rad                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_ls_cl_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ls_cl_rad                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_ls_cf_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ls_cf_rad                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_cc_cl_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%cc_cl_rad                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_cc_cf_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%cc_cf_rad                                   &
    )
END IF
!
!   Cloud effective dimension diagnostics
!
IF ( LW_diag(j_lw)%L_ls_del_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ls_del_rad                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_ls_def_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ls_def_rad                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_cc_del_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%cc_del_rad                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_cc_def_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%cc_def_rad                                  &
    )
END IF
!
! Radiances
!
IF ( LW_diag(j_lw)%L_toa_radiance ) THEN
  CALL rad3d_inp(                                              &
     L_complete_North, L_complete_South, L_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, ES_space_interp, n_channel,            &
     LW_diag(j_lw)%toa_radiance                                &
    )
END IF
!
!   Absorptivity diagnostics:
!
IF ( LW_diag(j_lw)%L_cloud_absorptivity ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, cloud_levels,           &
    LW_diag(j_lw)%cloud_absorptivity                            &
    )
END IF
IF ( LW_diag(j_lw)%L_cloud_weight_absorptivity ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, cloud_levels,           &
    LW_diag(j_lw)%cloud_weight_absorptivity                     &
    )
END IF
IF ( LW_diag(j_lw)%L_ls_cloud_absorptivity ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, cloud_levels,           &
    LW_diag(j_lw)%ls_cloud_absorptivity                         &
    )
END IF
IF ( LW_diag(j_lw)%L_ls_cloud_weight_absorptivity ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, cloud_levels,           &
    LW_diag(j_lw)%ls_cloud_weight_absorptivity                  &
    )
END IF
IF ( LW_diag(j_lw)%L_cnv_cloud_absorptivity ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, cloud_levels,           &
    LW_diag(j_lw)%cnv_cloud_absorptivity                        &
    )
END IF
IF ( LW_diag(j_lw)%L_cnv_cloud_weight_absorptivity ) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, cloud_levels,           &
    LW_diag(j_lw)%cnv_cloud_weight_absorptivity                 &
    )
END IF

! CLASSIC aerosol optical depth diagnostics

IF ( LW_diag(j_lw)%L_aod_sulphate) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_sulphate                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_dust) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_dust                                      &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_seasalt) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_seasalt                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_soot) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_soot                                      &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_biomass) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_biomass                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_biogenic) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_biogenic                                  &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_ocff) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_ocff                                      &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_delta) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_delta                                     &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_nitrate) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_nitrate                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_prog_sulphate) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_prog_sulphate                             &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_prog_dust) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_prog_dust                                 &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_prog_seasalt) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_prog_seasalt                              &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_prog_soot) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_prog_soot                                 &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_prog_biomass) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_prog_biomass                              &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_prog_ocff) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_prog_ocff                                 &
    )
END IF
IF ( LW_diag(j_lw)%L_aod_prog_nitrate) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_prog_nitrate                              &
    )
END IF
!
! CLASSIC absorption aerosol optical depth diagnostics
!
IF ( LW_diag(j_lw)%L_aaod_sulphate) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_sulphate                                 &
    )
END IF
IF ( LW_diag(j_lw)%L_aaod_dust) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_dust                                     &
    )
END IF
IF ( LW_diag(j_lw)%L_aaod_seasalt) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_seasalt                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_aaod_soot) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_soot                                     &
    )
END IF
IF ( LW_diag(j_lw)%L_aaod_biomass) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_biomass                                   &
    )
END IF
IF ( LW_diag(j_lw)%L_aaod_biogenic) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_biogenic                                 &
    )
END IF
IF ( LW_diag(j_lw)%L_aaod_ocff) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_ocff                                      &
    )
END IF
IF ( LW_diag(j_lw)%L_aaod_nitrate) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_nitrate                                  &
    )
END IF
!
! UKCA aerosol optical depth diagnostics
!
IF (LW_diag(j_lw)%L_aod_ukca_ait_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_ukca_ait_sol                              &
    )
END IF
IF (LW_diag(j_lw)%L_aod_ukca_acc_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_ukca_acc_sol                              &
    )
END IF
IF (LW_diag(j_lw)%L_aod_ukca_cor_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_ukca_cor_sol                              &
    )
END IF
IF (LW_diag(j_lw)%L_aod_ukca_ait_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_ukca_ait_ins                              &
    )
END IF
IF (LW_diag(j_lw)%L_aod_ukca_acc_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_ukca_acc_ins                              &
    )
END IF
IF (LW_diag(j_lw)%L_aod_ukca_cor_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aod_ukca_cor_ins                              &
    )
END IF
!
! UKCA stratospheric aerosol optical depth diagnostics
!
IF (LW_diag(j_lw)%L_sod_ukca_ait_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%sod_ukca_ait_sol                              &
    )
END IF
IF (LW_diag(j_lw)%L_sod_ukca_acc_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%sod_ukca_acc_sol                              &
    )
END IF
IF (LW_diag(j_lw)%L_sod_ukca_cor_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%sod_ukca_cor_sol                              &
    )
END IF
IF (LW_diag(j_lw)%L_sod_ukca_ait_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%sod_ukca_ait_ins                              &
    )
END IF
IF (LW_diag(j_lw)%L_sod_ukca_acc_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%sod_ukca_acc_ins                              &
    )
END IF
IF (LW_diag(j_lw)%L_sod_ukca_cor_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%sod_ukca_cor_ins                              &
    )
END IF
!
! UKCA absorption aerosol optical depth diagnostics
!
IF (LW_diag(j_lw)%L_aaod_ukca_ait_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_ukca_ait_sol                             &
    )
END IF
IF (LW_diag(j_lw)%L_aaod_ukca_acc_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_ukca_acc_sol                             &
    )
END IF
IF (LW_diag(j_lw)%L_aaod_ukca_cor_sol) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_ukca_cor_sol                             &
    )
END IF
IF (LW_diag(j_lw)%L_aaod_ukca_ait_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_ukca_ait_ins                             &
    )
END IF
IF (LW_diag(j_lw)%L_aaod_ukca_acc_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_ukca_acc_ins                             &
    )
END IF
IF (LW_diag(j_lw)%L_aaod_ukca_cor_ins) THEN
  CALL rad3d_inp(                                               &
    L_complete_North, L_complete_South, L_complete_deg,         &
    row_length, rows, off_x, off_y, first_row, last_row,        &
    first_data_interp, ES_space_interp, n_aod_wavel,            &
    LW_diag(j_lw)%aaod_ukca_cor_ins                             &
    )
END IF
!
! Total aerosol optical properties
!
IF ( lw_diag(j_lw)%l_aerosol_optical_depth ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels,         &
      lw_diag(j_lw)%aerosol_optical_depth(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_aerosol_scat_optical_depth ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels,         &
      lw_diag(j_lw)%aerosol_scat_optical_depth(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_aerosol_asymmetry_scat ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels,         &
      lw_diag(j_lw)%aerosol_asymmetry_scat(:,:,:,iwv) )
  END DO
END IF

! Band-by-band Fluxes
IF ( lw_diag(j_lw)%l_flux_up_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_up_band(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_flux_down_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_down_band(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_flux_up_clear_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_up_clear_band(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_flux_down_clear_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_down_clear_band(:,:,:,iwv) )
  END DO
END IF

! Clean-air Fluxes
IF ( lw_diag(j_lw)%l_flux_up_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_up_clean )
END IF
IF ( lw_diag(j_lw)%l_flux_down_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_down_clean )
END IF
IF ( lw_diag(j_lw)%l_flux_up_clear_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_up_clear_clean )
END IF
IF ( lw_diag(j_lw)%l_flux_down_clear_clean ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp,                       &
     model_levels+1, lw_diag(j_lw)%flux_down_clear_clean )
END IF
!
! Convective core diagnostics
!
IF ( LW_diag(j_lw)%l_ccore_clt_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ccore_clt_rad                               &
    )
END IF
IF ( LW_diag(j_lw)%l_ccore_qcl_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ccore_qcl_rad                               &
    )
END IF
IF ( LW_diag(j_lw)%l_ccore_qcf_rad ) THEN
  CALL rad3d_inp(                                             &
    L_complete_North, L_complete_South, L_complete_deg,       &
    row_length, rows, off_x, off_y, first_row, last_row,      &
    first_data_interp, ES_space_interp, model_levels,         &
    LW_diag(j_lw)%ccore_qcf_rad                               &
    )
END IF

! GHG forcing Fluxes
IF ( lw_diag(j_lw)%l_flux_up_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     lw_diag(j_lw)%flux_up_forc )
END IF
IF ( lw_diag(j_lw)%l_flux_down_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     lw_diag(j_lw)%flux_down_forc )
END IF
IF ( lw_diag(j_lw)%l_flux_up_clear_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     lw_diag(j_lw)%flux_up_clear_forc )
END IF
IF ( lw_diag(j_lw)%l_flux_down_clear_forc ) THEN
  CALL rad3d_inp(                                              &
     l_complete_north, l_complete_south, l_complete_deg,       &
     row_length, rows, off_x, off_y, first_row, last_row,      &
     first_data_interp, es_space_interp, model_levels+1,       &
     lw_diag(j_lw)%flux_down_clear_forc )
END IF
IF ( lw_diag(j_lw)%l_flux_up_forc_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_up_forc_band(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_flux_down_forc_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_down_forc_band(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_flux_up_clear_forc_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_up_clear_forc_band(:,:,:,iwv) )
  END DO
END IF
IF ( lw_diag(j_lw)%l_flux_down_clear_forc_band ) THEN
  DO iwv = 1, lw_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                             &
      l_complete_north, l_complete_south, l_complete_deg,       &
      row_length, rows, off_x, off_y, first_row, last_row,      &
      first_data_interp, es_space_interp, model_levels+1,       &
      lw_diag(j_lw)%flux_down_clear_forc_band(:,:,:,iwv) )
  END DO
END IF
!
! CLASSIC 3D aerosol-radiation diagnostics
!
IF (LW_diag(j_lw)%n_clas_aerosol_ext > 0) THEN
  DO iwv = 1, LW_diag(j_lw)%n_clas_aerosol_ext
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%clas_aerosol_ext(:,:,:,iwv)                   &
      )
  END DO
END IF
IF (LW_diag(j_lw)%n_clas_aerosol_abs > 0) THEN
  DO iwv = 1, LW_diag(j_lw)%n_clas_aerosol_abs
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%clas_aerosol_abs(:,:,:,iwv)                   &
      )
  END DO
END IF
IF (LW_diag(j_lw)%n_clas_aerosol_sca > 0) THEN
  DO iwv = 1, LW_diag(j_lw)%n_clas_aerosol_sca
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%clas_aerosol_sca(:,:,:,iwv)                   &
      )
  END DO
END IF
!
! UKCA 3D aerosol-radiation diagnostics
!
IF (LW_diag(j_lw)%n_ukca_aerosol_ext > 0) THEN
  DO iwv = 1, LW_diag(j_lw)%n_ukca_aerosol_ext
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%ukca_aerosol_ext(:,:,:,iwv)                   &
      )
  END DO
END IF
IF (LW_diag(j_lw)%n_ukca_aerosol_abs > 0) THEN
  DO iwv = 1, LW_diag(j_lw)%n_ukca_aerosol_abs
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%ukca_aerosol_abs(:,:,:,iwv)                   &
      )
  END DO
END IF
IF (LW_diag(j_lw)%n_ukca_aerosol_sca > 0) THEN
  DO iwv = 1, LW_diag(j_lw)%n_ukca_aerosol_sca
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%ukca_aerosol_sca(:,:,:,iwv)                   &
      )
  END DO
END IF
IF (LW_diag(j_lw)%n_ukca_aerosol_gsca > 0) THEN
  DO iwv = 1, LW_diag(j_lw)%n_ukca_aerosol_gsca
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%ukca_aerosol_gsca(:,:,:,iwv)                  &
      )
  END DO
END IF
IF (LW_diag(j_lw)%l_easyaerosol_extinction) THEN
  DO iwv = 1, LW_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%easyaerosol_extinction(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF
IF (LW_diag(j_lw)%l_easyaerosol_absorption) THEN
  DO iwv = 1, LW_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%easyaerosol_absorption(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF
IF (LW_diag(j_lw)%l_easyaerosol_scattering) THEN
  DO iwv = 1, LW_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%easyaerosol_scattering(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF
IF (LW_diag(j_lw)%l_easyaerosol_asytimscat) THEN
  DO iwv = 1, LW_spectrum(j_lw)%basic%n_band
    CALL rad3d_inp(                                               &
      L_complete_North, L_complete_South, L_complete_deg,         &
      row_length, rows, off_x, off_y, first_row, last_row,        &
      first_data_interp, ES_space_interp, model_levels,           &
      LW_diag(j_lw)%easyaerosol_asytimscat(:,:,:,iwv)             &
      )
  END DO ! iwv
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE fill_missing_data_lw
END MODULE fill_missing_data_lw_mod
