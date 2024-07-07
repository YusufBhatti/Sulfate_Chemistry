! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! Update convection diagnostics
!

MODULE update_conv_diags_mod

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Description:
!   Updates convection diagnostics after call to glue each substep of
!   convection.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.2.
! ------------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='UPDATE_CONV_DIAGS_MOD'

CONTAINS

SUBROUTINE  update_conv_diags(rows, row_length,                              &
            n_conv_levels, call_number, n_conv_calls,                        &
          ntml,ntpar, freeze_lev, it_kterm_deep, it_kterm_shall, it_cg_term, &
          cumulus, l_shallow, l_congestus, l_congestus2,it_dp_cfl_limited,   &
          it_md_cfl_limited, it_mid_level,                                   &
          one_over_conv_calls, timestep_conv , wstar,                        &
          exner_theta_levels, z_rho,                                         &

         ! Input (values output by latest convection call)
         it_cape_out, it_conv_rain, it_conv_snow,                            &
         it_precip_dp, it_precip_sh, it_precip_md, it_precip_cg,             &
         ind_cape_reduced, cape_ts_used, it_ind_deep, it_ind_shall,          &
         it_wstar_up, it_mb1, it_mb2,                                        &
         dubydt_p, dvbydt_p,                                                 &
         it_up_flux, it_up_flux_half, it_dwn_flux, it_entrain_up,            &
         it_entrain_dwn, it_detrain_up, it_detrain_dwn,                      &
         it_conv_rain_3d, it_conv_snow_3d,                                   &
         it_wqt_flux,it_wthetal_flux,it_wthetav_flux,it_wql_flux,            &
         it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall, it_uw_mid, it_vw_mid, &
         it_mf_deep, it_mf_congest, it_mf_shall, it_mf_midlev,               &
         it_dt_deep, it_dt_congest, it_dt_shall, it_dt_midlev,               &
         it_dq_deep, it_dq_congest, it_dq_shall, it_dq_midlev,               &
         it_du_deep, it_du_congest, it_du_shall, it_du_midlev,               &
         it_dv_deep, it_dv_congest, it_dv_shall, it_dv_midlev,               &
         it_w2p,     it_dt_dd,      it_dq_dd,   it_du_dd, it_dv_dd,          &
         it_area_ud, it_area_dd,                                             &
         ! in/out diagnostics
         dubydt_pout, dvbydt_pout, conv_rain, conv_snow )


! Modules

USE atm_fields_bounds_mod, ONLY:                                            &
   tdims_s, pdims_s

USE cv_run_mod,  ONLY:                                                      &
   l_mom,                                                                   &
   i_convection_vn,                                              &
   i_convection_vn_5a,                                           &
   i_convection_vn_6a

USE cv_stash_flg_mod, ONLY:                                                 &
  flg_up_flx, flg_up_flx_half, flg_dwn_flx,                                 &
  flg_entr_up, flg_entr_dwn, flg_detr_up, flg_detr_dwn,                     &
  flg_conv_rain_3d, flg_conv_snow_3d,                                       &
  flg_uw_dp, flg_vw_dp, flg_uw_shall, flg_vw_shall, flg_uw_mid, flg_vw_mid, &
  flg_wqt_flux, flg_wql_flux, flg_wthetal_flux, flg_wthetav_flux,           &
  flg_mf_deep, flg_mf_congest, flg_mf_shall, flg_mf_midlev,                 &
  flg_dt_deep, flg_dt_congest, flg_dt_shall, flg_dt_midlev,                 &
  flg_dq_deep, flg_dq_congest, flg_dq_shall, flg_dq_midlev,                 &
  flg_du_deep, flg_du_congest, flg_du_shall, flg_du_midlev,                 &
  flg_dv_deep, flg_dv_congest, flg_dv_shall, flg_dv_midlev,                 &
  flg_deep_tops, flg_w_eqn,                                                 &
  flg_dt_dd, flg_dq_dd, flg_du_dd, flg_dv_dd, flg_area_ud, flg_area_dd

! Convective diagnostic output arrays
USE cv_diagnostic_array_mod, ONLY:                                      &
        precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out    &
       ,deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind   &
       ,ntml_diag ,ntpar_diag ,freeze_diag ,kterm_diag                  &
       ,wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb        &
       ,wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top       &
       ,cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag           &
       ,uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux       &
       ,dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn      &
       ,up_flux_half                                                    &
       ,mf_deep ,mf_congest ,mf_shall ,mf_midlev                        &
       ,dt_deep ,dt_congest ,dt_shall ,dt_midlev                        &
       ,dq_deep ,dq_congest ,dq_shall ,dq_midlev                        &
       ,du_deep ,du_congest ,du_shall ,du_midlev                        &
       ,dv_deep ,dv_congest ,dv_shall ,dv_midlev                        &
       ,wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh      &
       ,conv_rain_3d ,conv_snow_3d                                      &
       ,deep_cfl_limited, mid_cfl_limited, deep_tops, w2p               &
       ,dt_dd, dq_dd, du_dd, dv_dd, area_ud, area_dd

USE nlsizes_namelist_mod,   ONLY: model_levels

USE model_domain_mod, ONLY: model_type, mt_single_column

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! ------------------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN) ::  &
  row_length            & ! row length
 ,rows                  & ! rows
 ,n_conv_levels         & ! Number of convection levels
 ,call_number           & ! call number for convection sub-steps
 ,n_conv_calls            ! number of convection calls

INTEGER, INTENT(IN) ::           &
  ntml(row_length, rows)         & ! Top level of surface mixed layer
 ,ntpar(row_length, rows)        & ! Top level of initial parcel ascent
 ,freeze_lev(row_length, rows)   & ! freezing level
 ,it_kterm_deep(row_length,rows) & ! lev no for terminating deep convection
 ,it_kterm_shall(row_length,rows) & ! lev no for terminating shallow convection
 ,it_cg_term(row_length,rows)      ! lev no for terminating congestus

LOGICAL, INTENT(IN) ::            &
  cumulus(row_length, rows)       & ! Cumulus present if true (from BL)
 ,l_shallow(row_length, rows)     & ! Logical switch for shallow Cu
 ,l_congestus(row_length, rows)   & ! Logical switch for congestus
 ,l_congestus2(row_length, rows)  & ! congestus in descending air
 ,it_mid_level(row_length, rows)    ! Mid level present in column

REAL, INTENT(IN) ::      &
  one_over_conv_calls    & ! 1/n_conv_calls
 ,timestep_conv            ! convection tiemstep (s)

REAL, INTENT(IN) ::          &
  wstar(row_length, rows)    & ! Sub-cloud convective velocity scale (m/s)
 ,exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)  &
 ,z_rho(row_length, rows, model_levels)                ! rho levels (m)

REAL, INTENT(IN) ::                 &
  it_cape_out(row_length, rows)     &
 ,it_conv_rain(row_length, rows)    &
 ,it_conv_snow(row_length, rows)    &
 ,it_precip_dp(row_length, rows)    &  ! deep precip
 ,it_precip_sh(row_length, rows)    &  ! shallow precip
 ,it_precip_md(row_length, rows)    &  ! mid-level precip
 ,it_precip_cg(row_length, rows)    &  ! congestus precip
 ,it_wstar_up(row_length, rows)     &
 ,it_mb1(row_length, rows)          &
 ,it_mb2(row_length, rows)          &
 ,ind_cape_reduced(row_length,rows) &  ! indicator of reduced cape timescale
 ,cape_ts_used(row_length, rows)    &  ! Actual cape timescale used for deep(s)
 ,it_dp_cfl_limited(row_length, rows) &! indicator for CFL limited deep conv
 ,it_md_cfl_limited(row_length, rows) &! indicator for CFL limited mid conv
 ,it_ind_deep(row_length, rows)       &! indicator of real deep
 ,it_ind_shall(row_length, rows)       ! indicator of real shallow

REAL, INTENT(IN) ::                                &
  it_up_flux(row_length, rows, model_levels)       &
 ,it_up_flux_half(row_length, rows, model_levels)  &  !up flux on half levs.
 ,it_dwn_flux(row_length, rows, model_levels)      &
 ,it_entrain_up(row_length, rows, model_levels)    &
 ,it_detrain_up(row_length, rows, model_levels)    &
 ,it_entrain_dwn(row_length, rows, model_levels)   &
 ,it_detrain_dwn(row_length, rows, model_levels)   &
 ,it_conv_rain_3d (row_length, rows, model_levels) &
 ,it_conv_snow_3d (row_length, rows, model_levels) &
 ,it_wqt_flux(row_length, rows, model_levels)      &
 ,it_wql_flux(row_length, rows, model_levels)      &
 ,it_wthetal_flux(row_length, rows, model_levels)  &
 ,it_wthetav_flux(row_length, rows, model_levels)  &
 ,it_uw_dp(row_length, rows, model_levels)         &
 ,it_vw_dp(row_length, rows, model_levels)         &
 ,it_uw_shall(row_length, rows, model_levels)      &
 ,it_vw_shall(row_length, rows, model_levels)      &
 ,it_uw_mid(row_length, rows, model_levels)        &
 ,it_vw_mid(row_length, rows, model_levels)        &
 ,it_mf_deep(row_length, rows, model_levels)                &
 ,it_mf_congest (row_length, rows, model_levels)            &
 ,it_mf_shall(row_length, rows, model_levels)               &
 ,it_mf_midlev(row_length, rows, model_levels)              &
 ,it_dt_deep(row_length, rows, model_levels)                &
 ,it_dt_congest(row_length, rows, model_levels)             &
 ,it_dt_shall(row_length, rows, model_levels)               &
 ,it_dt_midlev(row_length, rows, model_levels)              &
 ,it_dq_deep(row_length, rows, model_levels)                &
 ,it_dq_congest(row_length, rows, model_levels)             &
 ,it_dq_shall(row_length, rows, model_levels)               &
 ,it_dq_midlev (row_length, rows, model_levels)             &
 ,it_du_deep (row_length, rows, model_levels)               &
 ,it_du_congest(row_length, rows, model_levels)             &
 ,it_du_shall(row_length, rows, model_levels)               &
 ,it_du_midlev (row_length, rows, model_levels)             &
 ,it_dv_deep (row_length, rows, model_levels)               &
 ,it_dv_congest(row_length, rows, model_levels)             &
 ,it_dv_shall (row_length, rows, model_levels)              &
 ,it_dv_midlev (row_length, rows, model_levels)             &
 ,it_w2p(row_length, rows, model_levels)                    &
 ,it_dt_dd (row_length, rows, model_levels)                 &
 ,it_dq_dd (row_length, rows, model_levels)                 &
 ,it_du_dd (row_length, rows, model_levels)                 &
 ,it_dv_dd (row_length, rows, model_levels)                 &
 ,it_area_ud (row_length, rows, model_levels)               &
 ,it_area_dd (row_length, rows, model_levels)

REAL, INTENT(IN) ::                         &
  dubydt_p(row_length, rows, model_levels)  &
 ,dvbydt_p(row_length, rows, model_levels)

! Total convective momentum transport tendencies on p-grid
REAL, INTENT(INOUT) :: dubydt_pout ( pdims_s%i_start:pdims_s%i_end,      &
                                     pdims_s%j_start:pdims_s%j_end,      &
                                     pdims_s%k_start:pdims_s%k_end )
REAL, INTENT(INOUT) :: dvbydt_pout ( pdims_s%i_start:pdims_s%i_end,      &
                                     pdims_s%j_start:pdims_s%j_end,      &
                                     pdims_s%k_start:pdims_s%k_end )


REAL,INTENT(INOUT) ::                     &
  conv_rain(row_length,rows)              & ! convective rainfall
, conv_snow(row_length,rows)                ! convective snowfall


! Local declarations:
INTEGER  ::         &
  i,j,k                ! loop counters

REAL     ::         &
  rsteps               ! 1./(Number of real shallow sub-steps)

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UPDATE_CONV_DIAGS'

!------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!------------------------------------------------------------------------------
! Diagnostics not on switches
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Single level fields
!------------------------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(NONE) PRIVATE(i, j, k, rsteps)                        &
!$OMP& SHARED(rows, row_length, n_conv_levels, call_number, n_conv_calls,     &
!$OMP&    ntml,ntpar, freeze_lev, it_kterm_deep, it_kterm_shall, it_cg_term,  &
!$OMP&    cumulus, l_shallow, l_congestus, l_congestus2,it_dp_cfl_limited,    &
!$OMP&    it_md_cfl_limited, it_mid_level,                                    &
!$OMP&    one_over_conv_calls, timestep_conv , wstar,                         &
!$OMP&    exner_theta_levels, z_rho,                                          &
!$OMP&    it_cape_out, it_conv_rain, it_conv_snow,                            &
!$OMP&    it_precip_dp, it_precip_sh, it_precip_md, it_precip_cg,             &
!$OMP&    ind_cape_reduced, cape_ts_used, it_ind_deep, it_ind_shall,          &
!$OMP&    it_wstar_up, it_mb1, it_mb2,                                        &
!$OMP&    dubydt_p, dvbydt_p,                                                 &
!$OMP&    it_up_flux, it_up_flux_half, it_dwn_flux, it_entrain_up,            &
!$OMP&    it_entrain_dwn, it_detrain_up, it_detrain_dwn,                      &
!$OMP&    it_conv_rain_3d, it_conv_snow_3d,                                   &
!$OMP&    it_wqt_flux,it_wthetal_flux,it_wthetav_flux,it_wql_flux,            &
!$OMP&    it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall, it_uw_mid, it_vw_mid, &
!$OMP&    it_mf_deep, it_mf_congest, it_mf_shall, it_mf_midlev,               &
!$OMP&    it_dt_deep, it_dt_congest, it_dt_shall, it_dt_midlev,               &
!$OMP&    it_dq_deep, it_dq_congest, it_dq_shall, it_dq_midlev,               &
!$OMP&    it_du_deep, it_du_congest, it_du_shall, it_du_midlev,               &
!$OMP&    it_dv_deep, it_dv_congest, it_dv_shall, it_dv_midlev,               &
!$OMP&    it_w2p,     it_dt_dd,      it_dq_dd,   it_du_dd, it_dv_dd,          &
!$OMP&    it_area_ud, it_area_dd,                                             &
!$OMP&    conv_rain, conv_snow,  tdims_s, l_mom, i_convection_vn,             &
!$OMP&    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                           &
!$OMP&    flg_entr_up, flg_entr_dwn, flg_detr_up, flg_detr_dwn,               &
!$OMP&    flg_conv_rain_3d, flg_conv_snow_3d,flg_vw_mid,                      &
!$OMP&    flg_uw_dp, flg_vw_dp, flg_uw_shall, flg_vw_shall, flg_uw_mid,       &
!$OMP&    flg_wqt_flux, flg_wql_flux, flg_wthetal_flux, flg_wthetav_flux,     &
!$OMP&    flg_mf_deep, flg_mf_congest, flg_mf_shall, flg_mf_midlev,           &
!$OMP&    flg_dt_deep, flg_dt_congest, flg_dt_shall, flg_dt_midlev,           &
!$OMP&    flg_dq_deep, flg_dq_congest, flg_dq_shall, flg_dq_midlev,           &
!$OMP&    flg_du_deep, flg_du_congest, flg_du_shall, flg_du_midlev,           &
!$OMP&    flg_dv_deep, flg_dv_congest, flg_dv_shall, flg_dv_midlev,           &
!$OMP&    flg_deep_tops, flg_w_eqn, flg_area_dd,                              &
!$OMP&    flg_dt_dd, flg_dq_dd, flg_du_dd, flg_dv_dd, flg_area_ud,            &
!$OMP&    precip_deep, precip_shall ,precip_mid ,precip_cong ,cape_out,       &
!$OMP&    deep_ind ,shallow_ind ,congestus_ind ,congestus_ind2 ,mid_ind,      &
!$OMP&    ntml_diag ,ntpar_diag ,freeze_diag ,kterm_diag,                     &
!$OMP&    wstar_up_diag ,wstar_dn_diag ,mb1_diag ,mb2_diag ,wqt_cb,           &
!$OMP&    wthetal_cb ,wqt_inv ,wthetal_inv ,sh_top ,sh_base ,cg_top,          &
!$OMP&    cg_base ,cg_term ,cape_ts_diag ,ind_cape_reduced_diag,              &
!$OMP&    uw_dp, vw_dp ,uw_shall ,vw_shall ,uw_mid ,vw_mid ,up_flux,          &
!$OMP&    dwn_flux ,entrain_up ,detrain_up ,entrain_dwn ,detrain_dwn,         &
!$OMP&    up_flux_half,                                                       &
!$OMP&    mf_deep ,mf_congest ,mf_shall ,mf_midlev,                           &
!$OMP&    dt_deep ,dt_congest ,dt_shall ,dt_midlev,                           &
!$OMP&    dq_deep ,dq_congest ,dq_shall ,dq_midlev,                           &
!$OMP&    du_deep ,du_congest ,du_shall ,du_midlev,                           &
!$OMP&    dv_deep ,dv_congest ,dv_shall ,dv_midlev,                           &
!$OMP&    wqt_flux_sh ,wql_flux_sh ,wthetal_flux_sh ,wthetav_flux_sh,         &
!$OMP&    dubydt_pout ,dvbydt_pout ,conv_rain_3d ,conv_snow_3d,               &
!$OMP&    deep_cfl_limited, mid_cfl_limited, deep_tops, w2p,                  &
!$OMP&    dt_dd, dq_dd, du_dd, dv_dd, area_ud, area_dd,model_levels,          &
!$OMP&    model_type )


!$OMP DO SCHEDULE(STATIC)
DO j = 1, rows
  DO i = 1, row_length
    ntml_diag(i,j)  = ntml_diag(i,j)                                         &
                         + REAL(ntml(i,j))*one_over_conv_calls
    ntpar_diag(i,j) = ntpar_diag(i,j)                                        &
                         + REAL(ntpar(i,j))*one_over_conv_calls
    freeze_diag(i,j) = freeze_diag(i,j)                                      &
                         + REAL(freeze_lev(i,j)) *one_over_conv_calls

    wstar_dn_diag(i,j) = wstar_dn_diag(i,j) + wstar(i,j)*one_over_conv_calls
  END DO
END DO
!$OMP END DO NOWAIT

! In the SCM diagnostic, value is not mean over the 2 substeps, just
! final substep value, in order to get same answer as before.
! So for now overwrite in SCM runs to preserve KGO.
IF ( model_type == mt_single_column .AND. call_number == n_conv_calls ) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = 1, rows
    DO i = 1, row_length
      freeze_diag(i,j) = REAL( freeze_lev(i,j) )
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

! dependent on type of convection
!$OMP DO SCHEDULE(STATIC)  
DO j = 1, rows
  DO i = 1, row_length

    IF (it_mid_level(i,j)) THEN
      mid_ind(i,j) = mid_ind(i,j)+one_over_conv_calls
    END IF

    ! deep_ind - now set dependent on whether deep convection really occurs
    ! in the deep convection scheme. Failed deep attempts will be set to zero
    ! and not counted.
    deep_ind(i,j) = deep_ind(i,j) + it_ind_deep(i,j)*one_over_conv_calls

    ! Dependent on whether shallow really occurs
    shallow_ind(i,j) = shallow_ind(i,j) + it_ind_shall(i,j)*one_over_conv_calls

    ! Only work out top and bottom for points where real shallow
    ! convection occurred
    IF (it_ind_shall(i,j) == 1.0) THEN
      ! Use actual termination level rather than ntpar (5A & 6A)
      k=it_kterm_shall(i,j)+1
      sh_top(i,j) =sh_top(i,j) + z_rho(i,j,k)

      k=ntml(i,j)+1
      sh_base(i,j)=sh_base(i,j)+ z_rho(i,j,k)

    END IF
  END DO
END DO
!$OMP END DO NOWAIT

! Final sub-step
IF (call_number == n_conv_calls) THEN
  ! correctly scale heights by number of sub-steps with real shallow

!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      IF (shallow_ind(i,j) > 0.0) THEN   ! there are shallow sub-steps
        rsteps=1.0/(shallow_ind(i,j)*n_conv_calls)
        sh_top(i,j) =sh_top(i,j)*rsteps
        sh_base(i,j)=sh_base(i,j)*rsteps
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF



! Congestus type only defined for 5A & 6A
IF (i_convection_vn == i_convection_vn_5a .OR.    &
    i_convection_vn == i_convection_vn_6a ) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      IF (l_congestus(i,j)) THEN
        congestus_ind(i,j) = congestus_ind(i,j) + one_over_conv_calls

        k=it_cg_term(i,j)+1
        cg_top(i,j) =cg_top(i,j) + z_rho(i,j,k) *one_over_conv_calls

        k=ntml(i,j)+1
        cg_base(i,j)=cg_base(i,j)+ z_rho(i,j,k) *one_over_conv_calls

      END IF
      IF (l_congestus2(i,j)) THEN
        congestus_ind2(i,j) = congestus_ind2(i,j)  +one_over_conv_calls
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF


!$OMP DO SCHEDULE(STATIC)  
DO j = 1, rows
  DO i = 1, row_length

    conv_rain(i,j) = conv_rain(i,j) + it_conv_rain(i,j) * one_over_conv_calls

    conv_snow(i,j) = conv_snow(i,j) + it_conv_snow(i,j) * one_over_conv_calls

    precip_deep(i,j) = precip_deep(i,j) + it_precip_dp(i,j)                  &
                                                  * one_over_conv_calls
    precip_shall(i,j)= precip_shall(i,j) + it_precip_sh(i,j)                 &
                                                  * one_over_conv_calls
    precip_mid(i,j)  = precip_mid(i,j) + it_precip_md(i,j)                    &
                                                  * one_over_conv_calls
    cape_out(i,j)   = cape_out(i,j) +it_cape_out(i,j)*one_over_conv_calls

    kterm_diag(i,j) = kterm_diag(i,j)+                        &
                        REAL(it_kterm_deep(i,j)) *one_over_conv_calls

    ind_cape_reduced_diag(i,j) = ind_cape_reduced_diag(i,j) + &
                         ind_cape_reduced(i,j) *one_over_conv_calls
    cape_ts_diag(i,j) = cape_ts_diag(i,j) +                   &
                         cape_ts_used(i,j) *one_over_conv_calls

    deep_cfl_limited(i,j) = deep_cfl_limited(i,j) +                         &
                                 it_dp_cfl_limited(i,j) *one_over_conv_calls
    mid_cfl_limited(i,j)  = mid_cfl_limited(i,j) +                          &
                                 it_md_cfl_limited(i,j) *one_over_conv_calls

  END DO
END DO
!$OMP END DO NOWAIT

! 5A & 6A only
IF (i_convection_vn == i_convection_vn_5a .OR.    &
    i_convection_vn == i_convection_vn_6a ) THEN

!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      wstar_up_diag(i,j) = wstar_up_diag(i,j)                                 &
                             +it_wstar_up(i,j)*one_over_conv_calls
      mb1_diag(i,j) = mb1_diag(i,j)                                           &
                         +it_mb1(i,j)*one_over_conv_calls
      mb2_diag(i,j) = mb2_diag(i,j)                                           &
                             +it_mb2(i,j)*one_over_conv_calls
      cg_term(i,j) = cg_term(i,j)+                                            &
                        REAL(it_cg_term(i,j)) *one_over_conv_calls

      precip_cong(i,j) = precip_cong(i,j) + it_precip_cg(i,j)                 &
                                                    * one_over_conv_calls
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

! ------------------------------------------------------------------------------
! Diagnostics on switches
! ------------------------------------------------------------------------------

IF (flg_w_eqn) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        w2p(i,j,k)  = w2p(i,j,k) + it_w2p(i,j,k) * one_over_conv_calls
      END DO      ! i (row_length)
    END DO      ! j (rows)
  END DO      ! k (model_levels)
!$OMP END DO NOWAIT
END IF

IF (flg_conv_rain_3d) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        conv_rain_3d(i,j,k) = conv_rain_3d(i,j,k)                             &
                               + it_conv_rain_3d(i,j,k) * one_over_conv_calls
      END DO      ! i (row_length)
    END DO      ! j (rows)
  END DO      ! k (model_levels)
!$OMP END DO NOWAIT
END IF

IF (flg_conv_snow_3d) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, model_levels
    DO j=1, rows
      DO i=1, row_length
        conv_snow_3d(i,j,k) = conv_snow_3d(i,j,k)                             &
                               + it_conv_snow_3d(i,j,k) * one_over_conv_calls
      END DO      ! i (row_length)
    END DO      ! j (rows)
  END DO      ! k (model_levels)
!$OMP END DO NOWAIT
END IF

IF (flg_deep_tops) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      ! Deep convection
      ! Altered condition as can get failed deep cases
      IF (it_ind_deep(i,j) == 1.0) THEN
        k=it_kterm_deep(i,j)
        IF (k > 0) THEN  ! in case still get a zero value
          deep_tops(i,j,k) = deep_tops(i,j,k) + one_over_conv_calls
        END IF
      END IF
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_mf_deep) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_deep(i,j,k) = mf_deep(i,j,k)+it_mf_deep(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_mf_congest) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_congest(i,j,k) = mf_congest(i,j,k)                                 &
                                    +it_mf_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_mf_shall) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_shall(i,j,k) = mf_shall(i,j,k)                                     &
                                      +it_mf_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_mf_midlev) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        mf_midlev(i,j,k) = mf_midlev(i,j,k)                                   &
                                      +it_mf_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dt_deep) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_deep(i,j,k) = dt_deep(i,j,k)                                       &
                       + ( it_dt_deep(i,j,k) * exner_theta_levels(i,j,k)      &
                       * one_over_conv_calls )
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dt_congest) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_congest(i,j,k) = dt_congest(i,j,k)                                 &
                          + (it_dt_congest(i,j,k) * exner_theta_levels(i,j,k) &
                          * one_over_conv_calls )
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dt_shall) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_shall(i,j,k) = dt_shall(i,j,k)                                     &
                        + (it_dt_shall(i,j,k) * exner_theta_levels(i,j,k)     &
                        * one_over_conv_calls )
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dt_midlev) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dt_midlev(i,j,k) = dt_midlev(i,j,k)                                   &
                         + (it_dt_midlev(i,j,k) * exner_theta_levels(i,j,k)   &
                         * one_over_conv_calls )
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dq_deep) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_deep(i,j,k) = dq_deep(i,j,k) +it_dq_deep(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dq_congest) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_congest(i,j,k) = dq_congest(i,j,k)                                 &
                                    +it_dq_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dq_shall) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_shall(i,j,k) = dq_shall(i,j,k)                                     &
                                      +it_dq_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dq_midlev) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dq_midlev(i,j,k) = dq_midlev(i,j,k)                                   &
                                     +it_dq_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_du_deep) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels+1
    DO j = 1, rows
      DO i = 1, row_length
        du_deep(i,j,k) =  du_deep(i,j,k) +it_du_deep(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_du_congest) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,n_conv_levels+1
    DO j = 1, rows
      DO i = 1, row_length
        du_congest(i,j,k) = du_congest(i,j,k)                                 &
                                    +it_du_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_du_shall) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        du_shall(i,j,k) = du_shall(i,j,k)                                     &
                                      +it_du_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_du_midlev) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        du_midlev(i,j,k) = du_midlev(i,j,k)                                   &
                                      +it_du_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dv_deep) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        dv_deep(i,j,k) = dv_deep(i,j,k) +it_dv_deep(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dv_congest) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        dv_congest(i,j,k) = dv_congest(i,j,k)                 &
                 +it_dv_congest(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dv_shall) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        dv_shall(i,j,k) = dv_shall(i,j,k)                     &
                  +it_dv_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_dv_midlev) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1,model_levels
    DO j = 1, rows
      DO i = 1, row_length
        dv_midlev(i,j,k) = dv_midlev(i,j,k)                   &
                +it_dv_midlev(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (l_mom) THEN
  ! Require fields whether diagnostics required or not now as
  ! dubydt_pout and dvbydt_pout hold total increments to U and V from
  ! convection and are used to update r_u and r_v in atmos_physics2
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dubydt_pout(i,j,k) =dubydt_pout(i,j,k)                &
                              +dubydt_p(i,j,k)*one_over_conv_calls
        dvbydt_pout(i,j,k) =dvbydt_pout(i,j,k)                &
                              +dvbydt_p(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_up_flx) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        up_flux(i,j,k)=up_flux(i,j,k)+it_up_flux(i,j,k)*      &
                                   one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_up_flx_half) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        up_flux_half(i,j,k)=up_flux_half(i,j,k)               &
            +it_up_flux_half(i,j,k)* one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_dwn_flx) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        dwn_flux(i,j,k)=dwn_flux(i,j,k)+it_dwn_flux(i,j,k)*   &
                                   one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_entr_up) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        entrain_up(i,j,k)=entrain_up(i,j,k)+                  &
               it_entrain_up(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_entr_dwn) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        entrain_dwn(i,j,k)=entrain_dwn(i,j,k)+                &
               it_entrain_dwn(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_detr_up) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        detrain_up(i,j,k)=detrain_up(i,j,k)+                  &
               it_detrain_up(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_detr_dwn) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        detrain_dwn(i,j,k)=detrain_dwn(i,j,k)+                &
               it_detrain_dwn(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_uw_dp) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        uw_dp(i,j,k)=uw_dp(i,j,k)+                            &
               it_uw_dp(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_vw_dp) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        vw_dp(i,j,k)=vw_dp(i,j,k)+                            &
               it_vw_dp(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_uw_shall) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        uw_shall(i,j,k)=uw_shall(i,j,k)+                      &
               it_uw_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_vw_shall) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        vw_shall(i,j,k)=vw_shall(i,j,k)+                      &
               it_vw_shall(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_uw_mid) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        uw_mid(i,j,k)=uw_mid(i,j,k)+                          &
                           it_uw_mid(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_vw_mid) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        vw_mid(i,j,k)=vw_mid(i,j,k)+                      &
                           it_vw_mid(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_wqt_flux) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        wqt_flux_sh(i,j,k)=wqt_flux_sh(i,j,k)+                &
               it_wqt_flux(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      k=ntml(i,j)+1
      wqt_cb(i,j)=wqt_cb(i,j)+                              &
             it_wqt_flux(i,j,k)*one_over_conv_calls
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      k=ntpar(i,j)+1
      wqt_inv(i,j)=wqt_inv(i,j)+                            &
             it_wqt_flux(i,j,k)*one_over_conv_calls
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_wql_flux) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        wql_flux_sh(i,j,k)=wql_flux_sh(i,j,k)+                &
               it_wql_flux(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_wthetal_flux) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        wthetal_flux_sh(i,j,k)=wthetal_flux_sh(i,j,k)+        &
               it_wthetal_flux(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      k=ntml(i,j)+1
      wthetal_cb(i,j)=wthetal_cb(i,j)+                      &
             it_wthetal_flux(i,j,k)*one_over_conv_calls
    END DO
  END DO
!$OMP END DO NOWAIT
!$OMP DO SCHEDULE(STATIC)  
  DO j = 1, rows
    DO i = 1, row_length
      k=ntpar(i,j)+1
      wthetal_inv(i,j)=wthetal_inv(i,j)+                    &
             it_wthetal_flux(i,j,k)*one_over_conv_calls
    END DO
  END DO
!$OMP END DO NOWAIT
END IF
IF (flg_wthetav_flux) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k = 1, n_conv_levels
    DO j = 1, rows
      DO i = 1, row_length
        wthetav_flux_sh(i,j,k)=wthetav_flux_sh(i,j,k)+         &
               it_wthetav_flux(i,j,k)*one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

! Downdraught diagnostics

IF (flg_dt_dd) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, n_conv_levels
    DO j=1, rows
      DO i=1, row_length
        dt_dd(i,j,k) = dt_dd(i,j,k)  + it_dt_dd(i,j,k) * one_over_conv_calls
      END DO      ! i (row_length)
    END DO      ! j (rows)
  END DO      ! k
!$OMP END DO NOWAIT
END IF

IF (flg_dq_dd) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, n_conv_levels
    DO j=1, rows
      DO i=1, row_length
        dq_dd(i,j,k) = dq_dd(i,j,k)  + it_dq_dd(i,j,k) * one_over_conv_calls
      END DO      ! i (row_length)
    END DO      ! j (rows)
  END DO      ! k
!$OMP END DO NOWAIT
END IF

IF (flg_du_dd) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, n_conv_levels
    DO j=1, rows
      DO i=1, row_length
        du_dd(i,j,k) = du_dd(i,j,k)  + it_du_dd(i,j,k) * one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_dv_dd) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, n_conv_levels
    DO j=1, rows
      DO i=1, row_length
        dv_dd(i,j,k) = dv_dd(i,j,k)  + it_dv_dd(i,j,k) * one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_area_ud) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, n_conv_levels
    DO j=1, rows
      DO i=1, row_length
        area_ud(i,j,k) = area_ud(i,j,k)                                      &
                                  + it_area_ud(i,j,k) * one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

IF (flg_area_dd) THEN
!$OMP DO SCHEDULE(STATIC)  
  DO k=1, n_conv_levels
    DO j=1, rows
      DO i=1, row_length
        area_dd(i,j,k) = area_dd(i,j,k)                                      &  
                                  + it_area_dd(i,j,k) * one_over_conv_calls
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!$OMP END PARALLEL

!-------------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
END SUBROUTINE update_conv_diags

END MODULE update_conv_diags_mod
