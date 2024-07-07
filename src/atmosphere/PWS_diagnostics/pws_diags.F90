! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE pws_diags_driver_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_DIAGS_DRIVER_MOD'

CONTAINS
!
! Purpose: Driver routine for pws diagnostics (migrated from FieldCalc utility)
!
! Programming standard: Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics
!
SUBROUTINE pws_diags_driver(stash_len,icode,cmessage)

USE atm_step_local, ONLY: stashwork20
USE um_stashcode_mod, ONLY: stashcode_pws_sec, stashcode_pws_windspeed10m, &
    stashcode_bl_sec, stashcode_bl_windspeed_10mb,                         &
    stashcode_pws_windspeedplev,                                           &
    stashcode_conv_sec, stashcode_conv_icao_base, stashcode_conv_icao_top, &
    stashcode_pws_divergence, stashcode_pws_rel_vorticity,                 &
    stashcode_pws_snow_prob, stashcode_pws_conv_cld_dep,                   &
    stashcode_pws_conv_icao_base, stashcode_pws_conv_icao_top,             &
    stashcode_pws_precip_sym, stashcode_bl_1p5m_temp,                      &
    stashcode_pws_max_wind_ub, stashcode_pws_max_wind_vb,                  &
    stashcode_pws_max_wind_pb, stashcode_pws_max_wind_base,                &
    stashcode_pws_max_wind_top,stashcode_pws_max_wind_icao,                &
    stashcode_pws_thickness500, stashcode_pws_thickness850,                &
    stashcode_proc_phys_sec, stashcode_phy_geopht,                         &
    stashcode_proc_dyn_sec, stashcode_dyn_wind_ub, stashcode_dyn_wind_vb,  &
    stashcode_pws_tropopause_ht, stashcode_pws_tropopause_temp,            &
    stashcode_pws_tropopause_press, stashcode_pws_tropopause_icao,         &
    stashcode_pws_freezing_ht, stashcode_pws_freezing_press,               &
    stashcode_pws_freezing_icao, stashcode_pws_isotherm_ms20_ht,           &
    stashcode_pws_isotherm_ms20_press, stashcode_pws_isotherm_ms20_icao,   &
    stashcode_pws_isotherm_ms70_ht, stashcode_pws_isotherm_ms70_press,     &
    stashcode_pws_isotherm_ms70_icao, stashcode_theta, stashcode_exner,    &
    stashcode_aero_total_dust, stashcode_pws_dustconc_surf,                &
    stashcode_pws_dustconc_5000, stashcode_aerosol_sec,                    &
    stashcode_pws_1p5m_vis_tot, stashcode_pws_1p5m_vis_dust,               &
    stashcode_bl_1p5m_temp, stashcode_bl_1p5m_vis_tot,                     &
    stashcode_pws_thermal_advec, stashcode_pws_zenithdelay,                &
    stashcode_bl_1p5m_vis_tot,                                             &
    stashcode_pws_cat_turb, stashcode_pws_max_cat,                         &
    stashcode_pws_max_cat_press, stashcode_pws_mtn_wave_turb,              &
    stashcode_pws_contrail_bot, stashcode_pws_contrail_top,                &
    stashcode_pws_cloudturb_pot_diag, stashcode_pws_icing_pot_diag,        &
    stashcode_pws_wafc_caturb,                                             &
    stashcode_pws_panofsky_turb, stashcode_pws_inv_richardson,             &
    stashcode_pws_ellrodt1_turb, stashcode_pws_upd_helicity_5k

USE atm_fields_bounds_mod
USE atm_fields_mod
USE missing_data_mod,      ONLY: imdi, rmdi
USE level_heights_mod,     ONLY: r_theta_levels
USE planet_constants_mod,  ONLY: planet_radius

USE jules_snow_mod, ONLY: nsmax
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE free_tracers_inputs_mod, ONLY: a_tracer_last, a_tracer_first
USE nlsizes_namelist_mod, ONLY:                                          &
    land_field, model_levels, n_cca_lev, n_rows, ntiles,                 &
    river_row_length, river_rows, row_length, rows, sm_levels,           &
    theta_off_size, tr_levels, tr_ukca, tr_vars
USE um_parvars
USE submodel_mod, ONLY: atmos_sm, atmos_im
USE stash_array_mod, ONLY: si, sf, len_stlist, stindex, stlist,           &
                           num_stash_levels, stash_levels
USE horiz_grid_mod,   ONLY: cartesian_grid
USE dump_headers_mod, ONLY: rh_z_top_theta, a_realhd, rh_rotlat, rh_rotlong
USE lbc_mod, ONLY: lenrima
USE nlsizes_namelist_mod, ONLY: len_tot, n_obj_d1_max, len_fixhd,         & 
                                len_dumphist, len1_lookup, mpp_len1_lookup

USE pws_diags_mod, ONLY: pws_wind_speed_10mb, flag_windspeed_10m,   &
                         pws_diags_dealloc,                         &
                         pws_precip_sym, flag_precip_sym,           &
                         pws_snow_prob, flag_snow_prob,             &
                         pws_conv_cld_dep,pws_conv_icao_base,       &
                         pws_conv_icao_top,flag_conv_cld_dep,       &
                         flag_divergence, flag_rel_vorticity,       &
                         flag_conv_cld_base, flag_conv_cld_top,     &
                         pws_wind_ub_200, pws_wind_vb_200,          &
                         pws_wind_ub_250, pws_wind_vb_250,          &
                         pws_wind_ub_300, pws_wind_vb_300,          &
                         pws_max_wind_ub, pws_max_wind_vb,          &
                         pws_max_wind_pb, pws_max_wind_base,        &
                         pws_max_wind_top,pws_max_wind_icao,        &
                         flag_max_wind_ub, flag_max_wind_vb,        &
                         flag_max_wind_pb, flag_max_wind_base,      &
                         flag_max_wind_top, flag_max_wind_icao,     &
                         flag_windspeed_plev, flag_thermal_advec,   &
                         pws_thickness,                             &
                         flag_thickness_500, flag_thickness_850,    &
                         pws_tropopause_ht, flag_tropopause_ht,     &
                         pws_tropopause_temp, flag_tropopause_temp, &
                         pws_tropopause_press,flag_tropopause_press,&
                         pws_tropopause_icao,flag_tropopause_icao,  &
                         pws_freezing_ht, pws_freezing_press,       &
                         flag_freezing_ht, flag_freezing_press,     &
                         flag_freezing_icao, pws_freezing_icao,     &
                         flag_isotherm_ms20_ht,                     &
                         flag_isotherm_ms20_press,                  &
                         flag_isotherm_ms20_icao,                   &
                         flag_isotherm_ms70_ht,                     &
                         flag_isotherm_ms70_press,                  &
                         flag_isotherm_ms70_icao,                   &
                         pws_dustconc_surf, pws_zen_tot_delay,      &
                         pws_dustconc_5000, pws_dustconc_tot,       &
                         flag_dustconc_surf, flag_dustconc_5000,    &
                         pws_bl_1p5m_temp, pws_bl_1p5m_vis_tot,     &
                         flag_1p5m_vis_tot, flag_1p5m_vis_dust,     &
                         pws_cat_turb, flag_cat_turb,               &
                         cat_press, cat_press_levs, uwind, vwind,   &
                         pws_max_cat, flag_max_cat,                 &
                         pws_max_cat_press, flag_max_cat_press,     &
                         pws_1p5m_vis_dust, pws_1p5m_vis_tot,       &
                         pws_contrail_bot, flag_contrail_bot,       &
                         pws_contrail_top, flag_contrail_top,       &
                         pws_cloudturb_pot_diag, flag_cloudturb_pot,&
                         flag_zenithdelay,                          &
                         pot_press, pot_press_levs,                 &
                         tropo_model_top, heightcut_top,            &
                         pws_icing_pot_diag, flag_icing_pot,        &
                         icing_pot_press_levs,                      &
                         pws_wafc_cat_diag, flag_wafc_cat,          &
                         wafc_cat_press_levs,                       &
                         pws_panofsky_turb, flag_panofsky_turb,     &
                         panofsky_turb_press_levels,                &
                         pws_inv_richardson, flag_inv_richardson,   &
                         inv_richardson_press_levels,               &
                         pws_mtn_wave_turb, flag_mtn_wave_turb,     &
                         pws_ellrodt1_turb, flag_ellrodt1_turb,     &
                         ellrodt1_turb_press_levels,                &
                         pws_upd_helicity_5k, flag_upd_helicity_5k

USE pws_max_winds_topbase_mod, ONLY: pws_max_winds_topbase
USE pws_tropoht_mod, ONLY: pws_tropoht
USE pws_cloudturb_pot_mod, ONLY: pws_cloudturb_pot
USE pws_contrail_mod, ONLY: pws_contrail
USE pws_colson_panofsky_turb_mod, ONLY: pws_colson_panofsky_turb
USE pws_icing_pot_mod,  ONLY: pws_icing_pot
USE pws_wafc_cat_mod, ONLY: pws_wafc_cat
USE pws_inv_richardson_number_mod, ONLY: pws_inv_richardson_number
USE pws_mtn_turb_mod, ONLY: pws_mtn_turb
USE pws_ellrod1_mod, ONLY: pws_ellrod1
USE pws_vertlev_choices_mod, ONLY: MX_NumLevs_Gl, pws_vertlev_choices_init

USE pws_isotherm_mod, ONLY: pws_isotherm
USE water_constants_mod, ONLY: tm  ! freezing/melting temp of water.

USE pws_cat_mod, ONLY: pws_cat
USE pws_dust_diag_mod, ONLY: pws_dust_diag
USE pws_vis_diag_mod, ONLY: pws_vis_diag, pws_vis2_diag
USE dust_parameters_mod, ONLY: i_dust, i_dust_off, l_twobin_dust
USE pws_updh_diag_mod, ONLY: pws_updh_diag

USE t_vert_interp_to_p_mod, ONLY: t_vert_interp_to_p
USE model_domain_mod, ONLY: model_type, mt_global, l_regular

! lookup_addresses in control/dump_io
! blev - zsea, bhlev - c, lblev - level number
USE lookup_addresses, ONLY: item_code, blev, bhlev, lblev 

USE umPrintMgr ! ummessage

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

USE zenith_delay_mod, ONLY: zenith_delay
USE uv_plevsB_mod, ONLY: uv_plevsB
USE therm_advec_mod, ONLY: therm_advec
USE temp_plevs_mod, ONLY: temperature_plevs

USE pws_thickness_diag_mod, ONLY: pws_thickness_diag
USE pws_snow_prob_diag_mod, ONLY: pws_snow_prob_diag
USE pws_rel_vortic_diag_mod, ONLY: pws_rel_vortic_diag
USE pws_precip_sym_diag_mod, ONLY: pws_precip_sym_diag
USE pws_divergence_diag_mod, ONLY: pws_divergence_diag

IMPLICIT NONE


INTEGER :: stash_len
INTEGER :: istind     ! stash index id
INTEGER :: nlevlist   ! stash levels list number
INTEGER :: nplevs_pws ! number of p levels for diags on p levels
INTEGER, PARAMETER :: ereport_pws_stash = -666
! icode for ereport that is reserved for when a stash item is
! missing a dependency.
LOGICAL :: missing_stash_dependency

! Output p levels for diags on p levels
REAL, ALLOCATABLE :: plevs_pws(:) 

! Wind speed on p levels B grid
REAL, ALLOCATABLE :: pws_windspeed_plev(:,:,:)
! u,v  on p levels B grid
REAL, ALLOCATABLE :: pws_uwind_plev    (:,:,:)
REAL, ALLOCATABLE :: pws_vwind_plev    (:,:,:)
! divergence, rel vorticity on P levels B-grid
REAL, ALLOCATABLE :: pws_divergence    (:,:,:)
REAL, ALLOCATABLE :: pws_rel_vorticity (:,:,:)

! Thermal advection on p levels B-grid
REAL, ALLOCATABLE :: pws_thermal_advec(:,:,:)
! Temperature advection on p levels B-grid
REAL, ALLOCATABLE :: pws_temp_plev    (:,:,:)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: lev_num
INTEGER :: i, j, k
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName = 'PWS_DIAGS_DRIVER'

LOGICAL :: rotated

INTEGER :: im_index

REAL :: pressure_pa


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! set up level dept variables for use by this routine.
CALL pws_vertlev_choices_init(a_realhd(rh_z_top_theta))

ALLOCATE (STASHwork20(stash_len))

STASHwork20(:) = 0.0

icode = 0 
missing_stash_dependency = .FALSE.
im_index = 1

! Check for rotated grids
! Non-rotated global model has pole at (90,0)
! Non-rotated LAM has pole at either (90,0) or (90,180)
rotated = (a_realhd(rh_rotlat) /= 90.0 .OR.                                   &
            (.NOT. (a_realhd(rh_rotlong) == 0.0                               &
               .OR. a_realhd(rh_rotlong) == 180.0))                           &
          .AND. .NOT. cartesian_grid)

IF (flag_windspeed_10m) THEN

  IF (.NOT. sf(stashcode_bl_windspeed_10mb,stashcode_bl_sec) ) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20003 PWS windspeed at 10m not available unless '         &
         //'stash 3227 windspeed at 10m B grid activated'
    CALL ereport(RoutineName,icode,cmessage)
  ELSE

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_windspeed10m,                  &
                            stashcode_pws_sec,im_index)), pws_wind_speed_10mb,&
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity,  &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_windspeed10m, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, 10m wind speed) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

END IF


IF (flag_windspeed_plev) THEN

  ! position in stash list array
  istind=stindex(1, stashcode_pws_windspeedplev, stashcode_pws_sec, im_index)

  NPLevs_pws = 0

  IF (istind >  0) THEN
    ! stash levels list number
    nlevlist = -stlist(10,istind)

    ! Number of output levels
    NPLevs_pws = stash_levels(1,nlevlist)

    IF (NPLevs_pws > 0) THEN

      ALLOCATE (PLevs_pws(NPLevs_pws))

      ! Allocate B grid arrays
      ALLOCATE(pws_windspeed_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
      pws_windspeed_plev(:,:,:) = 0.0
      ALLOCATE(pws_uwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
      pws_uwind_plev(:,:,:) = 0.0
      ALLOCATE(pws_vwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
      pws_vwind_plev(:,:,:) = 0.0    

      ! Output p levels
      DO k = 1, NPLevs_pws
        ! Need to divide by 1000 to convert to hPa
        PLevs_pws(k) = stash_levels(k+1,nlevlist)/1000.0
      END DO
    ELSE
      icode = -1  ! Warning
      cmessage=": PWS wind speed on p levels called but no levels list "
      CALL ereport(routinename,icode,cmessage)
    END IF
  END IF

  IF (icode==0) THEN

    ! Interp u, v on model levels C grid to p levels B grid
    CALL uv_plevsB(NPLevs_pws, PLevs_pws, u, v, pws_uwind_plev, pws_vwind_plev)

    ! Compute wind speeds
    pws_windspeed_plev(:,:,:) = pws_uwind_plev(:,:,:) * pws_uwind_plev(:,:,:)+ &
                                pws_vwind_plev(:,:,:) * pws_vwind_plev(:,:,:)

    pws_windspeed_plev(:,:,:) = SQRT( pws_windspeed_plev(:,:,:) )

    ! Copy wind speeds into stashwork array
    DO  k = 1, NPLevs_pws                    !loop over output p levels
      DO j = vdims%j_start, vdims%j_end
        DO i = udims%i_start, udims%i_end
          stashwork20(                                                       & 
              si(stashcode_pws_windspeedplev, stashcode_pws_sec, im_index) + & 
                        (i-udims%i_start) + (j-vdims%j_start)*row_length   + &
                        (k-1)*n_rows*row_length)                             &
                                                  = pws_windspeed_plev(i,j,k)
        END DO
      END DO
    END DO

  END IF

  IF (ALLOCATED(PLevs_pws))          DEALLOCATE(PLevs_pws)
  IF (ALLOCATED(pws_uwind_plev))     DEALLOCATE(pws_uwind_plev)
  IF (ALLOCATED(pws_vwind_plev))     DEALLOCATE(pws_vwind_plev)
  IF (ALLOCATED(pws_windspeed_plev)) DEALLOCATE(pws_windspeed_plev)  

END IF


IF (flag_thermal_advec) THEN

  IF (rotated .OR. .NOT. l_regular) THEN

    icode = -1        !Warning
    cmessage='Stash 20031 thermal advection pressure levels unavailable ' &
    //'on rotated grids or variable res grids'

    CALL ereport(RoutineName,icode,cmessage)
 
  ELSE

    ! position in stash list array
    istind=stindex(1, stashcode_pws_thermal_advec, stashcode_pws_sec, im_index)

    NPLevs_pws = 0

    IF (istind >  0) THEN
      ! stash levels list number
      nlevlist = -stlist(10,istind)

      ! Number of output levels
      NPLevs_pws = stash_levels(1,nlevlist)

      IF (NPLevs_pws > 0) THEN

        ALLOCATE (PLevs_pws(NPLevs_pws))
 
        ! Allocate B grid arrays
        ALLOCATE(pws_thermal_advec &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, NPLevs_pws))
        pws_thermal_advec(:,:,:) = 0.0
        ALLOCATE(pws_uwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
        pws_uwind_plev(:,:,:) = 0.0
        ALLOCATE(pws_vwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
        pws_vwind_plev(:,:,:) = 0.0 
        ALLOCATE(pws_temp_plev &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, NPLevs_pws))
        pws_temp_plev(:,:,:) = 0.0

        ! Output p levels
        DO k = 1, NPLevs_pws
          ! Need to divide by 1000 to convert to hPa
          PLevs_pws(k) = stash_levels(k+1,nlevlist)/1000.0
        END DO
      ELSE
        icode = -1 ! Warning
        cmessage=":PWS thermal advection on p levels called but no levels list "
        CALL ereport(routinename,icode,cmessage)
      END IF  
    
      IF (icode==0) THEN

        ! Interp u, v on model levels C grid to p levels B grid
        CALL uv_plevsB(NPLevs_pws,PLevs_pws,u,v,pws_uwind_plev,pws_vwind_plev)

        ! Compute temperature on p levels
        CALL temperature_plevs(NPlevs_pws, plevs_pws, pws_temp_plev)

        ! Compute thermal advection on p levels
        CALL therm_advec (NPlevs_pws, pws_uwind_plev, pws_vwind_plev, &
                                    pws_temp_plev, pws_thermal_advec)

        ! Copy thermal advection into stashwork array
        DO  k = 1, NPLevs_pws                    !loop over output p levels
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              stashwork20(                                                    & 
              si(stashcode_pws_thermal_advec, stashcode_pws_sec, im_index) + & 
                        (i-tdims%i_start) + (j-tdims%j_start)*row_length   + &
                        (k-1)*rows*row_length)                               &
                                                  = pws_thermal_advec(i,j,k)
            END DO
          END DO
        END DO

      END IF

    END IF

    IF (ALLOCATED(PLevs_pws))          DEALLOCATE(PLevs_pws)
    IF (ALLOCATED(pws_thermal_advec))  DEALLOCATE(pws_thermal_advec)
    IF (ALLOCATED(pws_uwind_plev))     DEALLOCATE(pws_uwind_plev)
    IF (ALLOCATED(pws_vwind_plev))     DEALLOCATE(pws_vwind_plev)
    IF (ALLOCATED(pws_temp_plev ))     DEALLOCATE(pws_temp_plev )

  END IF

END IF


IF (flag_precip_sym) THEN

  IF (.NOT. sf(stashcode_bl_1p5m_temp,stashcode_bl_sec) ) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20014 PWS precip symbol not available unless '         &
         //'stash 3236 temp at 1p5m activated'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

    CALL pws_precip_sym_diag(STASHwork20,                                   &
                         icode,cmessage)

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_precip_sym,                  &
                            stashcode_pws_sec,im_index)), pws_precip_sym,   &
                            row_length, rows, 0,0,0,0, at_extremity,  &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_precip_sym, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, precip symbol) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

END IF



IF (flag_max_wind_ub .OR. flag_max_wind_vb .OR. flag_max_wind_pb            &
    .OR. flag_max_wind_base .OR. flag_max_wind_top                          &
    .OR. flag_max_wind_icao) THEN

  CALL pws_max_winds_topbase(u,v,p,  MX_NumLevs_Gl,                         &
                             icode,cmessage)


  IF (flag_max_wind_ub) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_max_wind_ub,                 &
                            stashcode_pws_sec,im_index)), pws_max_wind_ub,  &
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity,&
                            atmos_im, stashcode_pws_sec,                    &
                            stashcode_pws_max_wind_ub, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, max_wind_ub) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

  IF (flag_max_wind_vb) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_max_wind_vb,                 &
                            stashcode_pws_sec,im_index)), pws_max_wind_vb,  &
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity,&
                            atmos_im, stashcode_pws_sec,                    &
                            stashcode_pws_max_wind_vb, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, max_wind_vb) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

  IF (flag_max_wind_pb) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_max_wind_pb,                 &
                            stashcode_pws_sec,im_index)), pws_max_wind_pb,  &
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity,&
                            atmos_im, stashcode_pws_sec,                    &
                            stashcode_pws_max_wind_pb, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, max_wind_pb) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

  IF (flag_max_wind_base) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_max_wind_base,               &
                            stashcode_pws_sec,im_index)), pws_max_wind_base,&
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity,&
                            atmos_im, stashcode_pws_sec,                    &
                            stashcode_pws_max_wind_base, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, max_wind_base) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

  IF (flag_max_wind_top) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_max_wind_top,                &
                            stashcode_pws_sec,im_index)), pws_max_wind_top, &
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity,&
                            atmos_im, stashcode_pws_sec,                    &
                            stashcode_pws_max_wind_top, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, max_wind_top) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

  IF (flag_max_wind_icao) THEN
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_max_wind_icao,               &
                            stashcode_pws_sec,im_index)), pws_max_wind_icao,&
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity,&
                            atmos_im, stashcode_pws_sec,                    &
                            stashcode_pws_max_wind_icao, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, max_wind_icao) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF


END IF


IF (flag_thickness_500 .OR. flag_thickness_850 .OR. flag_snow_prob) THEN

  IF (.NOT. sf(stashcode_phy_geopht,stashcode_proc_phys_sec) ) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20001/2/28 PWS thickness, snow prob not available ' &
    //'unless stash 16202 geopotential heights on pressure levels activated'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

    IF (flag_thickness_500) THEN
      CALL pws_thickness_diag(20001,                                   &
                         icode,cmessage)

! DEPENDS ON: copydiag
      CALL copydiag(stashwork20(si(stashcode_pws_thickness500,             &
                            stashcode_pws_sec,im_index)), pws_thickness,   &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_thickness500, icode, cmessage)

      IF (icode  >   0) THEN
        cmessage=": error in copydiag(PWS diags, 1000-500mb thickness) "
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF

    IF (flag_thickness_850) THEN
      CALL pws_thickness_diag(20002,                                   &
                         icode,cmessage)

! DEPENDS ON: copydiag
      CALL copydiag(stashwork20(si(stashcode_pws_thickness850,             &
                            stashcode_pws_sec,im_index)), pws_thickness,   &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_thickness850, icode, cmessage)

      IF (icode  >   0) THEN
        cmessage=": error in copydiag(PWS diags, 1000-850mb thickness) "
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF

    IF (flag_snow_prob) THEN
      CALL pws_snow_prob_diag(20028,                                   &
                         icode,cmessage)

! DEPENDS ON: copydiag
      CALL copydiag(stashwork20(si(stashcode_pws_snow_prob,                &
                            stashcode_pws_sec,im_index)), pws_snow_prob,   &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_snow_prob, icode, cmessage)

      IF (icode  >   0) THEN
        cmessage=": error in copydiag(PWS diags, snow probability) "
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF

  END IF

END IF


IF (flag_contrail_bot .OR. flag_contrail_top) THEN

  CALL pws_contrail(theta, p_theta_levels,exner_theta_levels,pstar)

END IF

IF (flag_contrail_bot) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_contrail_bot,              &
                            stashcode_pws_sec,im_index)), pws_contrail_bot,   &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_contrail_bot, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, contrail lowest level) "
      CALL ereport(routinename,icode,cmessage)
    END IF

END IF

IF (flag_contrail_top) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_contrail_top,              &
                            stashcode_pws_sec,im_index)), pws_contrail_top,   &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_contrail_top, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, contrail lowest level) "
      CALL ereport(routinename,icode,cmessage)
    END IF

END IF

IF (flag_icing_pot)   THEN

  CALL pws_icing_pot(p_theta_levels,theta,exner_theta_levels,  &
                           cf_bulk, q)

  DO k = 1, icing_pot_press_levs

    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        stashwork20(si(stashcode_pws_icing_pot_diag,stashcode_pws_sec,  &
                 im_index)+ (i-tdims%i_start)+(j-tdims%j_start)*row_length+ &
                  (k-1)*rows*row_length) = pws_icing_pot_diag(i,j,k)
      END DO
    END DO

  END DO

END IF



IF (flag_tropopause_ht.OR.flag_tropopause_temp.OR.flag_tropopause_press.OR. &
        flag_tropopause_icao ) THEN
 
  tropo_model_top = imdi

  ! Initialise tropo_model_top to a model level close to 22km
  DO k=1, model_levels
    IF ( (r_theta_levels(tdims%i_start,tdims%j_start,k+1)-planet_radius    &
          > heightcut_top) ) THEN
      tropo_model_top = k
      EXIT
    END IF
  END DO

  ! if model lid is too low we do not bother to hunt for tropopause.
  IF (tropo_model_top /= imdi) THEN   
    CALL pws_tropoht(theta, p_theta_levels, exner_theta_levels)
  END IF

END IF

IF (flag_tropopause_ht) THEN

  IF (tropo_model_top /= imdi) THEN 
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_tropopause_ht,              &
                            stashcode_pws_sec,im_index)), pws_tropopause_ht,   &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_tropopause_ht, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, tropopause height) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  ELSE
    icode = -1        !Warning
    cmessage='Stash 20064 tropopause height is unavailable ' &
    //'as model lid may be too low for tropopause algorithm'

    CALL ereport(RoutineName,icode,cmessage)
  END IF

END IF

IF (flag_tropopause_temp) THEN

  IF (tropo_model_top /= imdi) THEN 

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_tropopause_temp,              &
                            stashcode_pws_sec,im_index)), pws_tropopause_temp,&
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_tropopause_temp, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, tropopause temperature) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  ELSE
    icode = -1        !Warning
    cmessage='Stash 20065 tropopause temp is unavailable ' &
    //'as model lid may be too low for tropopause algorithm'

    CALL ereport(RoutineName,icode,cmessage)
  END IF
END IF

IF (flag_tropopause_press) THEN

  IF (tropo_model_top /= imdi) THEN 
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_tropopause_press,              &
                            stashcode_pws_sec,im_index)), pws_tropopause_press,&
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_tropopause_press, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, tropopause pressure) "
      CALL ereport(routinename,icode,cmessage)
    END IF
  ELSE
    icode = -1        !Warning
    cmessage='Stash 20066 tropopause press is unavailable ' &
    //'as model lid may be too low for tropopause algorithm'

    CALL ereport(RoutineName,icode,cmessage)
  END IF

END IF

IF (flag_tropopause_icao) THEN

  IF (tropo_model_top /= imdi) THEN 
! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_tropopause_icao,              &
                            stashcode_pws_sec,im_index)), pws_tropopause_icao,&
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_tropopause_icao, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, tropopause pressure) "
      CALL ereport(routinename,icode,cmessage)
    END IF
  ELSE
    icode = -1        !Warning
    cmessage='Stash 20067 tropopause ICAO height is unavailable ' &
    //'as model lid may be too low for tropopause algorithm'

    CALL ereport(RoutineName,icode,cmessage)
  END IF

END IF

IF (flag_wafc_cat)   THEN

  IF (rotated .OR. .NOT. l_regular) THEN

    icode = -1        !Warning
    cmessage='Stash 20047 WAFC Cat potential is unavailable ' &
        //'on rotated grids or variable res grids'
    
    CALL ereport(RoutineName,icode,cmessage)
    
  ELSE

    CALL pws_wafc_cat()

    DO k = 1, wafc_cat_press_levs
      DO j = vdims%j_start, vdims%j_end
        DO i = udims%i_start, udims%i_end
          stashwork20(si(stashcode_pws_wafc_caturb,stashcode_pws_sec,       &
                 im_index)+ (i-udims%i_start)+(j-vdims%j_start)*row_length+ &
                 (k-1)*n_rows*row_length) = pws_wafc_cat_diag(i,j,k)
        END DO
      END DO
    END DO

  END IF ! rotated and l_regular
END IF ! flag_wafc_cat

IF (flag_freezing_ht .OR. flag_freezing_press .OR.                         &
    flag_freezing_icao ) THEN

  CALL pws_isotherm(tm, orography,pstar,theta,p_theta_levels,              &
                    exner_theta_levels,flag_freezing_icao)

  IF (flag_freezing_ht) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_freezing_ht,                &
                            stashcode_pws_sec,im_index)), pws_freezing_ht, &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_freezing_ht, icode, cmessage)

  END IF

  IF (flag_freezing_press) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_freezing_press,                &
                            stashcode_pws_sec,im_index)), pws_freezing_press, &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_freezing_press, icode, cmessage)

  END IF

  IF (flag_freezing_icao) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_freezing_icao,                 &
                            stashcode_pws_sec,im_index)), pws_freezing_icao,  &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_freezing_icao, icode, cmessage)

  END IF


END IF


IF (flag_isotherm_ms20_ht .OR. flag_isotherm_ms20_press .OR.                  &
    flag_isotherm_ms20_icao) THEN

  CALL pws_isotherm(tm-20.0, orography,pstar,theta,p_theta_levels,            &
                    exner_theta_levels,flag_isotherm_ms20_icao)

  IF (flag_isotherm_ms20_ht) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_isotherm_ms20_ht,              &
                        stashcode_pws_sec,im_index)), pws_freezing_ht,        &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_isotherm_ms20_ht, icode, cmessage)

  END IF

  IF (flag_isotherm_ms20_press) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_isotherm_ms20_press,           &
                      stashcode_pws_sec,im_index)), pws_freezing_press,       &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_isotherm_ms20_press, icode, cmessage)

  END IF

  IF (flag_isotherm_ms20_icao) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_isotherm_ms20_icao,            &
                       stashcode_pws_sec,im_index)), pws_freezing_icao,       &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_isotherm_ms20_icao, icode, cmessage)

  END IF


END IF


IF (flag_isotherm_ms70_ht .OR. flag_isotherm_ms70_press .OR.                  &
    flag_isotherm_ms70_icao) THEN

  CALL pws_isotherm(tm-70.0, orography,pstar,theta,p_theta_levels,            &
                    exner_theta_levels,flag_isotherm_ms70_icao)

  IF (flag_isotherm_ms70_ht) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_isotherm_ms70_ht,              &
                        stashcode_pws_sec,im_index)), pws_freezing_ht,        &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_isotherm_ms70_ht, icode, cmessage)

  END IF

  IF (flag_isotherm_ms70_press) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_isotherm_ms70_press,           &
                      stashcode_pws_sec,im_index)), pws_freezing_press,       &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_isotherm_ms70_press, icode, cmessage)

  END IF

  IF (flag_isotherm_ms70_icao) THEN

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_isotherm_ms70_icao,            &
                       stashcode_pws_sec,im_index)), pws_freezing_icao,       &
                            row_length, rows, 0,0,0,0, at_extremity,          &
                            atmos_im, stashcode_pws_sec,                      &
                            stashcode_pws_isotherm_ms70_icao, icode, cmessage)

  END IF


END IF

IF (flag_cloudturb_pot)   THEN

  DO k = 1, pot_press_levs

    pressure_pa = pot_press(k)*100.0   ! convert to Pascals

    CALL pws_cloudturb_pot(p_theta_levels,theta,exner_theta_levels,  &
                           cf_bulk, pressure_pa,k)


    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        stashwork20(si(stashcode_pws_cloudturb_pot_diag,stashcode_pws_sec,  &
                 im_index)+ (i-tdims%i_start)+(j-tdims%j_start)*row_length+ &
                  (k-1)*rows*row_length) = pws_cloudturb_pot_diag(i,j,k)
      END DO
    END DO

  END DO

END IF


! ----------------------------------------------------------------------
! STASH item 16 : cat on pressure levels
! ----------------------------------------------------------------------

IF (flag_cat_turb .OR. flag_max_cat .OR. flag_max_cat_press) THEN

  IF (rotated .OR. .NOT. l_regular) THEN

     icode = -1        !Warning
    cmessage='Stash 20016-18 cat turb diags are unavailable '      &
    //'as they are invalid on rotated grids or variable res grids.'

    CALL ereport(RoutineName,icode,cmessage)

  ELSE IF (.NOT. sf(stashcode_dyn_wind_ub,stashcode_proc_dyn_sec) .OR. &
      .NOT. sf(stashcode_dyn_wind_vb,stashcode_proc_dyn_sec) ) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20016-18 cat turb are unavailable ' &
    //'unless stash 15201/2 u,v-winds on pressure levels also requested'

    CALL ereport(RoutineName,icode,cmessage)

  ELSE

    DO k = 1, cat_press_levs
    
      IF (NINT(cat_press(k)) == 300) THEN
        uwind(:,:)=pws_wind_ub_300(:,:)
        vwind(:,:)=pws_wind_vb_300(:,:)
      ELSE IF (NINT(cat_press(k)) == 250) THEN
        uwind(:,:)=pws_wind_ub_250(:,:)
        vwind(:,:)=pws_wind_vb_250(:,:)
      ELSE IF (NINT(cat_press(k)) == 200) THEN
        uwind(:,:)=pws_wind_ub_200(:,:)
        vwind(:,:)=pws_wind_vb_200(:,:)
      ELSE
        icode = -1        ! Warning
        cmessage='Stash 20016 cat turb not available ' &
        //'on requested model level; only on 300,250 and 200hPa'

        CALL ereport(RoutineName,icode,cmessage)
      END IF

      pressure_pa = cat_press(k)*100.0   ! convert to Pascals

      CALL pws_cat(u,v,p, uwind,vwind,pressure_pa,k)

      DO j = vdims%j_start, vdims%j_end
        DO i = udims%i_start, udims%i_end
          stashwork20(si(stashcode_pws_cat_turb,stashcode_pws_sec,            &
                   im_index)+ (i-udims%i_start)+(j-vdims%j_start)*row_length+ &
                   (k-1)*n_rows*row_length) = pws_cat_turb(i,j,k)
        END DO
      END DO

    END DO

    IF (flag_max_cat .OR. flag_max_cat_press) THEN

      pws_max_cat(:,:) = rmdi
      pws_max_cat_press (:,:) = rmdi

      DO k = 1, cat_press_levs
        WHERE ( pws_cat_turb(:,:,k) > pws_max_cat(:,:) )
          pws_max_cat(:,:) = pws_cat_turb(:,:,k)
          pws_max_cat_press (:,:) = cat_press(k)*100.0 ! Pa (SI unit)
        END WHERE
      END DO
    
      IF (flag_max_cat) THEN

        CALL copydiag(stashwork20(si(stashcode_pws_max_cat,                  &
                            stashcode_pws_sec,im_index)), pws_max_cat,       &
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity, &
                            atmos_im, stashcode_pws_sec,                     &
                            stashcode_pws_max_cat, icode, cmessage)

      END IF

      IF (flag_max_cat_press) THEN

        CALL copydiag(stashwork20(si(stashcode_pws_max_cat_press,            &
                            stashcode_pws_sec,im_index)), pws_max_cat_press, &
                            udims%i_len, vdims%j_len, 0,0,0,0, at_extremity, &
                            atmos_im, stashcode_pws_sec,                     &
                            stashcode_pws_max_cat_press, icode, cmessage)

      END IF

    END IF

  END IF

END IF


IF (flag_dustconc_surf) THEN

  IF (.NOT. sf(stashcode_aero_total_dust,stashcode_aerosol_sec) ) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20058 PWS surface dust conc not available unless '         &
         //'stash 17257 total dust conc activated'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

    ! Surface dust conc : convert from microg/m3 to g/m3
    pws_dustconc_surf(:,:) = pws_dustconc_surf(:,:) * 1.0e-6

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_dustconc_surf,                  &
                            stashcode_pws_sec,im_index)), pws_dustconc_surf,   &
                            row_length, rows, 0,0,0,0, at_extremity,           &
                            atmos_im, stashcode_pws_sec,                       &
                            stashcode_pws_dustconc_surf, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, surface dust conc) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

END IF

IF (flag_dustconc_5000) THEN

  IF (.NOT. sf(stashcode_aero_total_dust,stashcode_aerosol_sec) ) THEN
    icode = ereport_pws_stash        ! Warning
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20059 PWS dust conc 2000-5000ft not available unless '     &
         //'stash 17257 total dust conc activated'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

    CALL pws_dust_diag(pws_dustconc_tot, pws_dustconc_5000)

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_dustconc_5000,                  &
                            stashcode_pws_sec,im_index)), pws_dustconc_5000,   &
                            row_length, rows, 0,0,0,0, at_extremity,           &
                            atmos_im, stashcode_pws_sec,                       &
                            stashcode_pws_dustconc_5000, icode, cmessage)

    IF (icode  >   0) THEN
      cmessage=": error in copydiag(PWS diags, 2000-5000ft dust conc) "
      CALL ereport(routinename,icode,cmessage)
    END IF

  END IF

END IF

IF (flag_1p5m_vis_tot .OR. flag_1p5m_vis_dust) THEN

  IF (.NOT. (sf(stashcode_bl_1p5m_temp,stashcode_bl_sec) .AND. &
             sf(stashcode_bl_1p5m_vis_tot,stashcode_bl_sec) .AND. &
                                        .NOT.(i_dust==i_dust_off) )) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20070/20071 PWS 1.5m visibilities not available unless '   &
         //'aerosol dust option active and '                                   &
         //'stash 3236 1.5m temp and 3281 1.5m total vis active'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

    IF (l_twobin_dust) THEN
      CALL pws_vis2_diag(pws_bl_1p5m_vis_tot, pws_bl_1p5m_temp, &
                         pws_1p5m_vis_dust, pws_1p5m_vis_tot,   &
                         flag_1p5m_vis_tot, flag_1p5m_vis_dust)
    ELSE 
      ! 6 bin dust
      CALL pws_vis_diag(pws_bl_1p5m_vis_tot, pws_bl_1p5m_temp, &
                        pws_1p5m_vis_dust, pws_1p5m_vis_tot,   &
                        flag_1p5m_vis_tot, flag_1p5m_vis_dust)
    END IF
    
    IF (flag_1p5m_vis_tot) THEN
! DEPENDS ON: copydiag
      CALL copydiag(stashwork20(si(stashcode_pws_1p5m_vis_tot,                 &
                            stashcode_pws_sec,im_index)), pws_1p5m_vis_tot, &
                            row_length, rows, 0,0,0,0, at_extremity,           &
                            atmos_im, stashcode_pws_sec,                       &
                            stashcode_pws_1p5m_vis_tot, icode, cmessage)

      IF (icode  >   0) THEN
        cmessage=": error in copydiag(PWS diags, 1.5m total visibility) "
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF

    IF (flag_1p5m_vis_dust) THEN
! DEPENDS ON: copydiag
      CALL copydiag(stashwork20(si(stashcode_pws_1p5m_vis_dust,                &
                            stashcode_pws_sec,im_index)), pws_1p5m_vis_dust,   &
                            row_length, rows, 0,0,0,0, at_extremity,           &
                            atmos_im, stashcode_pws_sec,                       &
                            stashcode_pws_1p5m_vis_dust, icode, cmessage)

      IF (icode  >   0) THEN
        cmessage=": error in copydiag(PWS diags, 1.5m dust visibility) "
        CALL ereport(routinename,icode,cmessage)
      END IF
    END IF

  END IF

END IF

IF (flag_conv_cld_dep) THEN

  IF (.NOT. sf(stashcode_conv_icao_base,stashcode_conv_sec) .OR. &
      .NOT. sf(stashcode_conv_icao_top,stashcode_conv_sec) ) THEN
    icode = ereport_pws_stash        ! Warning
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20012 convective cloud depth not available ' &
    //'unless stash 5210/211 conv icao base/top also chosen'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

    pws_conv_cld_dep(:,:) = pws_conv_icao_top(:,:)-pws_conv_icao_base(:,:)

    WHERE (pws_conv_icao_top(:,:) == rmdi .OR. pws_conv_icao_base(:,:) == rmdi)
      pws_conv_cld_dep(:,:) = rmdi
    END WHERE

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_conv_cld_dep,               &
                     stashcode_pws_sec,im_index)), pws_conv_cld_dep(:,:),  &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_conv_cld_dep, icode, cmessage)

  END IF

END IF

IF (flag_conv_cld_base) THEN

  IF (.NOT. sf(stashcode_conv_icao_base,stashcode_conv_sec)) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20039 convective cloud base not available ' &
    //'unless stash 5210 conv icao base also chosen'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_conv_icao_base,             &
                     stashcode_pws_sec,im_index)), pws_conv_icao_base(:,:),&
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_conv_icao_base, icode, cmessage)

  END IF

END IF

IF (flag_conv_cld_top) THEN

  IF (.NOT. sf(stashcode_conv_icao_top,stashcode_conv_sec)) THEN
    icode = ereport_pws_stash
    missing_stash_dependency = .TRUE.
    cmessage='Stash 20040 convective cloud top not available ' &
    //'unless stash 5211 conv icao top also chosen'

    CALL ereport(RoutineName,icode,cmessage)
  ELSE

! DEPENDS ON: copydiag
    CALL copydiag(stashwork20(si(stashcode_pws_conv_icao_top,              &
                     stashcode_pws_sec,im_index)), pws_conv_icao_top(:,:), &
                            row_length, rows, 0,0,0,0, at_extremity,       &
                            atmos_im, stashcode_pws_sec,                   &
                            stashcode_pws_conv_icao_top, icode, cmessage)

  END IF

END IF

IF (flag_divergence .OR. flag_rel_vorticity) THEN

  IF (rotated .OR. .NOT. l_regular) THEN

    icode = -1        ! Warning
    cmessage='Stash 20005/6 divergence, rel vorticity on pressure levels ' &
    //'unavailable on rotated grids or variable res grids'

    CALL ereport(RoutineName,icode,cmessage)

  ELSE

    IF (flag_divergence) THEN

      ! position in stash list array
      istind=stindex(1, stashcode_pws_divergence, stashcode_pws_sec, im_index)

      NPLevs_pws = 0

      IF (istind >  0) THEN
        ! stash levels list number
        nlevlist = -stlist(10,istind)

        ! Number of output levels
        NPLevs_pws = stash_levels(1,nlevlist)

        IF (NPLevs_pws > 0) THEN

          ALLOCATE (PLevs_pws(NPLevs_pws))

          ! Allocate B grid arrays
          ALLOCATE(pws_divergence &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
          pws_divergence(:,:,:) = 0.0
          ALLOCATE(pws_uwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
          pws_uwind_plev(:,:,:) = 0.0
          ALLOCATE(pws_vwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
          pws_vwind_plev(:,:,:) = 0.0

          ! Output p levels
          DO k = 1, NPLevs_pws
            ! Need to divide by 1000 to convert to hPa
            PLevs_pws(k) = stash_levels(k+1,nlevlist)/1000.0
          END DO
        ELSE
          icode = -1 ! Warning
          cmessage=":PWS Divergence on p levels called but no levels list"
          CALL ereport(routinename,icode,cmessage)
        END IF 

        IF (icode==0) THEN

          ! Interp u, v on model levels C grid to p levels B grid
          CALL uv_plevsB(NPLevs_pws,PLevs_pws,u,v,pws_uwind_plev,pws_vwind_plev)


          CALL pws_divergence_diag(pws_uwind_plev, pws_vwind_plev,        &
              NPlevs_pws, pws_divergence(:,:,:), icode, cmessage)

          DO  k=1, NPLevs_pws                    !loop over output p levels
            DO j = vdims%j_start, vdims%j_end
              DO i = udims%i_start, udims%i_end
                stashwork20(si(stashcode_pws_divergence,stashcode_pws_sec, &
                im_index)+ (i-udims%i_start)+(j-vdims%j_start)*row_length+ &
                (k-1)*n_rows*row_length) = pws_divergence(i,j,k)
              END DO
            END DO
          END DO
        END IF
      END IF
      IF (ALLOCATED(PLevs_pws))          DEALLOCATE(PLevs_pws)
      IF (ALLOCATED(pws_divergence))     DEALLOCATE(pws_divergence)
      IF (ALLOCATED(pws_uwind_plev))     DEALLOCATE(pws_uwind_plev)
      IF (ALLOCATED(pws_vwind_plev))     DEALLOCATE(pws_vwind_plev)
    END IF

    IF (flag_rel_vorticity) THEN

      ! position in stash list array
      istind=stindex(1,stashcode_pws_rel_vorticity,stashcode_pws_sec,im_index)

      NPLevs_pws = 0

      IF (istind >  0) THEN
        ! stash levels list number
        nlevlist = -stlist(10,istind)

        ! Number of output levels
        NPLevs_pws = stash_levels(1,nlevlist)

        IF (NPLevs_pws > 0) THEN

          ALLOCATE (PLevs_pws(NPLevs_pws))

          ! Allocate B grid arrays
          ALLOCATE(pws_rel_vorticity &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
          pws_rel_vorticity(:,:,:) = 0.0
          ALLOCATE(pws_uwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
          pws_uwind_plev(:,:,:) = 0.0
          ALLOCATE(pws_vwind_plev &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end, NPLevs_pws))
          pws_vwind_plev(:,:,:) = 0.0

          ! Output p levels
          DO k = 1, NPLevs_pws
            ! Need to divide by 1000 to convert to hPa
            PLevs_pws(k) = stash_levels(k+1,nlevlist)/1000.0
          END DO
        ELSE
          icode = -1 ! Warning
          cmessage=":PWS Rel Vorticity on p levels called but no levels list"
          CALL ereport(routinename,icode,cmessage)
        END IF 

        IF (icode==0) THEN

          ! Interp u, v on model levels C grid to p levels B grid
          CALL uv_plevsB(NPLevs_pws,PLevs_pws,u,v,pws_uwind_plev,pws_vwind_plev)

          CALL pws_rel_vortic_diag(pws_uwind_plev, pws_vwind_plev,        &
              NPLevs_pws, pws_rel_vorticity(:,:,:), icode, cmessage)

          DO  k=1, NPLevs_pws                    !loop over output p levels
            DO j = vdims%j_start, vdims%j_end
              DO i = udims%i_start, udims%i_end
                stashwork20(si(stashcode_pws_rel_vorticity,stashcode_pws_sec, &
                im_index)+ (i-udims%i_start)+(j-vdims%j_start)*row_length+ &
                (k-1)*n_rows*row_length) = pws_rel_vorticity(i,j,k)
              END DO
            END DO
          END DO

        END IF

      END IF

      IF (ALLOCATED(PLevs_pws))          DEALLOCATE(PLevs_pws)
      IF (ALLOCATED(pws_rel_vorticity))  DEALLOCATE(pws_rel_vorticity)
      IF (ALLOCATED(pws_uwind_plev))     DEALLOCATE(pws_uwind_plev)
      IF (ALLOCATED(pws_vwind_plev))     DEALLOCATE(pws_vwind_plev)

    END IF

  END IF

END IF

IF (flag_zenithdelay) THEN
 
  CALL zenith_delay (pws_zen_tot_delay)

! DEPENDS ON: copydiag
  CALL copydiag(stashwork20(si(stashcode_pws_zenithdelay,                    &
                            stashcode_pws_sec,im_index)), pws_zen_tot_delay, &
                            tdims%i_len, tdims%j_len, 0,0,0,0, at_extremity, &
                            atmos_im, stashcode_pws_sec,                     &
                            stashcode_pws_zenithdelay, icode, cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(PWS diags, zenith total delay) "
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF


IF (flag_panofsky_turb) THEN
  CALL pws_colson_panofsky_turb(u, v, theta, p_theta_levels)

  DO k = 1, panofsky_turb_press_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        stashwork20(si(stashcode_pws_panofsky_turb, stashcode_pws_sec,    &
                   im_index) + (i-tdims%i_start) + &
                   (j-tdims%j_start)*row_length+ &
                   (k-1)*rows*row_length) = pws_panofsky_turb(i,j,k)
      END DO
    END DO
  END DO

END IF

IF (flag_ellrodt1_turb) THEN
  CALL pws_ellrod1(u, v, p_theta_levels, p)

  DO k = 1, ellrodt1_turb_press_levels
    DO j = vdims%j_start, vdims%j_end
      DO i = udims%i_start, udims%i_end
        stashwork20(si(stashcode_pws_ellrodt1_turb, stashcode_pws_sec,    &
                   im_index) + (i-udims%i_start) + &
                   (j-vdims%j_start)*row_length+ &
                   (k-1)*n_rows*row_length) = pws_ellrodt1_turb(i,j,k)
      END DO
    END DO
  END DO

END IF

IF (flag_inv_richardson) THEN
  CALL pws_inv_richardson_number(u, v, theta, p_theta_levels)

  DO k = 1, inv_richardson_press_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        stashwork20(si(stashcode_pws_inv_richardson, stashcode_pws_sec,    &
                   im_index) + (i-tdims%i_start) + &
                   (j-tdims%j_start)*row_length+ &
                   (k-1)*rows*row_length) = pws_inv_richardson(i,j,k)
      END DO
    END DO
  END DO

END IF

IF (flag_mtn_wave_turb) THEN

  CALL pws_mtn_turb( icode, cmessage)

! DEPENDS ON: copydiag
  CALL copydiag(stashwork20(si(stashcode_pws_mtn_wave_turb,           &
                stashcode_pws_sec,im_index)), pws_mtn_wave_turb(:,:), &
                row_length, rows, 0,0,0,0, at_extremity,              &
                atmos_im, stashcode_pws_sec,                          &
                stashcode_pws_mtn_wave_turb, icode, cmessage)

END IF

IF (flag_upd_helicity_5k) THEN

! Stash 20080 PWS upd helicity 2-5km has gathered stash 30001,2,3 wind 
! components automatically by forcing routine st_diag3 to be called.

  CALL pws_updh_diag(pws_upd_helicity_5k, icode, cmessage)

! DEPENDS ON: copydiag
  CALL copydiag(stashwork20(si(stashcode_pws_upd_helicity_5k,                &
                          stashcode_pws_sec,im_index)), pws_upd_helicity_5k, &
                          row_length, rows, 0,0,0,0, at_extremity,           &
                          atmos_im, stashcode_pws_sec,                       &
                          stashcode_pws_upd_helicity_5k, icode, cmessage)

  IF (icode  >   0) THEN
    cmessage=": error in copydiag(PWS diags, 2-5km upd helicity) "
    CALL ereport(routinename,icode,cmessage)
  END IF

END IF


! Raise an error if missing stash dependencies
IF (missing_stash_dependency) THEN
  icode = -ereport_pws_stash
  WRITE(cmessage,'(A,I0)') 'Missing dependencies in the PWS diagnostics. '// &
                           'See messages with warning number ',              &
                           ereport_pws_stash 
  CALL ereport(RoutineName,icode,cmessage)
END IF

! DEPENDS ON: stash
CALL stash(atmos_sm, atmos_im, stashcode_pws_sec, stashwork20, &
             icode,cmessage)

CALL pws_diags_dealloc()

DEALLOCATE (STASHwork20)


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN

END SUBROUTINE pws_diags_driver

END MODULE pws_diags_driver_mod
