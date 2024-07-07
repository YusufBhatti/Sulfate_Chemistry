! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Purpose: Interface for pws diagnostics (migrated from FieldCalc utility)
!
! Programming standard: Unified Model Documentation Paper No. 3
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: PWS_diagnostics

MODULE pws_diags_mod

IMPLICIT NONE

! Arrays containing PWS diagnostics and their respective flags are
! declared in this module.
!
! Diag arrays should follow the naming convention
! TYPE, ALLOCATABLE :: pws_[name](:,:)
!
! Diag flags should follow the naming convention
! LOGICAL :: flag_[name] = .FALSE.
!
! Wind speed at 10m B-grid
REAL, ALLOCATABLE :: pws_wind_speed_10mb(:,:)
LOGICAL           :: flag_windspeed_10m = .FALSE.
! Precipitation symbol
REAL, ALLOCATABLE :: pws_precip_sym(:,:)
REAL, ALLOCATABLE :: pws_precip_sym_ls_rain(:,:)
REAL, ALLOCATABLE :: pws_precip_sym_ls_snow(:,:)
REAL, ALLOCATABLE :: pws_precip_sym_conv_rain(:,:)
REAL, ALLOCATABLE :: pws_precip_sym_conv_snow(:,:)
REAL, ALLOCATABLE :: pws_precip_sym_t1p5m(:,:)
LOGICAL           :: flag_precip_sym = .FALSE.
! Snow probability
REAL, ALLOCATABLE :: pws_snow_prob(:,:)
LOGICAL           :: flag_snow_prob = .FALSE.
! Convective cloud depth, base and top (ICAO ht)
REAL, ALLOCATABLE :: pws_conv_cld_dep(:,:)
REAL, ALLOCATABLE :: pws_conv_icao_base(:,:)
REAL, ALLOCATABLE :: pws_conv_icao_top(:,:)
LOGICAL           :: flag_conv_cld_dep = .FALSE.
LOGICAL           :: flag_conv_cld_base = .FALSE.
LOGICAL           :: flag_conv_cld_top = .FALSE.

! Flags for Divergence and Relative Vorticity on various pressure levels
LOGICAL           :: flag_divergence = .FALSE.
LOGICAL           :: flag_rel_vorticity = .FALSE.
! Maxwinds
REAL, ALLOCATABLE :: pws_max_wind_ub(:,:)
REAL, ALLOCATABLE :: pws_max_wind_vb(:,:)
REAL, ALLOCATABLE :: pws_max_wind_pb(:,:)
REAL, ALLOCATABLE :: pws_max_wind_base(:,:)
REAL, ALLOCATABLE :: pws_max_wind_top(:,:)
REAL, ALLOCATABLE :: pws_max_wind_icao(:,:)
LOGICAL           :: flag_max_wind_ub = .FALSE.
LOGICAL           :: flag_max_wind_vb = .FALSE.
LOGICAL           :: flag_max_wind_pb = .FALSE.
LOGICAL           :: flag_max_wind_base = .FALSE.
LOGICAL           :: flag_max_wind_top = .FALSE.
LOGICAL           :: flag_max_wind_icao = .FALSE.
! Output p levels for diags on p levels 
LOGICAL           :: flag_windspeed_plev = .FALSE.
! Thermal advection on p levels B-grid
LOGICAL           :: flag_thermal_advec = .FALSE.
! Thickness 1000-500mb and/or 1000-850mb
REAL, ALLOCATABLE :: pws_thickness(:,:)
REAL, ALLOCATABLE :: pws_geopht_1000(:,:)
REAL, ALLOCATABLE :: pws_geopht_850(:,:)
REAL, ALLOCATABLE :: pws_geopht_500(:,:)
LOGICAL           :: flag_thickness_500 = .FALSE.
LOGICAL           :: flag_thickness_850 = .FALSE.
! Tropopause diagnostics
REAL, ALLOCATABLE :: pws_tropopause_ht(:,:)
REAL, ALLOCATABLE :: pws_tropopause_press(:,:)
REAL, ALLOCATABLE :: pws_tropopause_temp(:,:)
REAL, ALLOCATABLE :: pws_tropopause_icao(:,:)

LOGICAL           ::  flag_tropopause_ht = .FALSE.
LOGICAL           ::  flag_tropopause_temp = .FALSE.
LOGICAL           ::  flag_tropopause_press = .FALSE.
LOGICAL           ::  flag_tropopause_icao = .FALSE.

! Icing potential
REAL, ALLOCATABLE :: pws_icing_pot_diag(:,:,:)
REAL, ALLOCATABLE :: icing_pot_press(:)
LOGICAL           :: flag_icing_pot = .FALSE.
INTEGER           :: icing_pot_press_levs

! Freezing Levels
REAL, ALLOCATABLE :: pws_freezing_ht(:,:)
REAL, ALLOCATABLE :: pws_freezing_press(:,:)
REAL, ALLOCATABLE :: pws_freezing_icao(:,:)
! zenith total delay
REAL, ALLOCATABLE :: pws_zen_tot_delay(:,:)
REAL, ALLOCATABLE :: h_rho_levels(:,:,:)
REAL, ALLOCATABLE :: h_theta_levels(:,:,:)
LOGICAL           :: flag_zenithdelay = .FALSE.

LOGICAL           :: flag_freezing_ht = .FALSE.
LOGICAL           :: flag_freezing_press = .FALSE.
LOGICAL           :: flag_freezing_icao = .FALSE.
LOGICAL           :: flag_isotherm_ms20_ht = .FALSE.
LOGICAL           :: flag_isotherm_ms20_press = .FALSE.
LOGICAL           :: flag_isotherm_ms20_icao = .FALSE.
LOGICAL           :: flag_isotherm_ms70_ht = .FALSE.
LOGICAL           :: flag_isotherm_ms70_press = .FALSE.
LOGICAL           :: flag_isotherm_ms70_icao = .FALSE.

! In cloud turb potential.
REAL, ALLOCATABLE :: pws_cloudturb_pot_diag(:,:,:)
LOGICAL           :: flag_cloudturb_pot = .FALSE.
INTEGER           :: pot_press_levs
REAL, ALLOCATABLE :: pot_press(:)

! WAFC caturb potential
REAL, ALLOCATABLE :: pws_wafc_cat_diag(:,:,:)
REAL, ALLOCATABLE :: wafc_cat_press(:)
LOGICAL           :: flag_wafc_cat = .FALSE.
INTEGER           :: wafc_cat_press_levs

REAL, ALLOCATABLE :: pws_gwd_stress_lev_u(:,:,:)
REAL, ALLOCATABLE :: pws_gwd_stress_lev_v(:,:,:)

! Cat turb
REAL, ALLOCATABLE :: pws_wind_ub_200(:,:)
REAL, ALLOCATABLE :: pws_wind_vb_200(:,:)
REAL, ALLOCATABLE :: pws_wind_ub_250(:,:)
REAL, ALLOCATABLE :: pws_wind_vb_250(:,:)
REAL, ALLOCATABLE :: pws_wind_ub_300(:,:)
REAL, ALLOCATABLE :: pws_wind_vb_300(:,:)
REAL, ALLOCATABLE :: uwind(:,:)
REAL, ALLOCATABLE :: vwind(:,:)
REAL, ALLOCATABLE :: pws_cat_turb(:,:,:)
REAL, ALLOCATABLE :: pws_max_cat(:,:)
REAL, ALLOCATABLE :: pws_max_cat_press(:,:)
LOGICAL           :: flag_cat_turb = .FALSE.
LOGICAL           :: flag_max_cat = .FALSE.
LOGICAL           :: flag_max_cat_press = .FALSE.
INTEGER           :: cat_press_levs
REAL, ALLOCATABLE :: cat_press(:)

! Mountain wave turbulence
REAL, ALLOCATABLE :: pws_mtn_wave_turb(:,:)
LOGICAL           :: flag_mtn_wave_turb = .FALSE.

! Dust concentrations
REAL, ALLOCATABLE :: pws_dustconc_tot (:,:,:)
REAL, ALLOCATABLE :: pws_dustconc_surf(:,:)
REAL, ALLOCATABLE :: pws_dustconc_5000(:,:)
LOGICAL           :: flag_dustconc_surf = .FALSE.
LOGICAL           :: flag_dustconc_5000 = .FALSE.
! Model level 1 dust mass-mixing-ratios (in bins), within timestep
! After emission:
REAL, ALLOCATABLE :: pws_dustmmr1_em (:,:)
REAL, ALLOCATABLE :: pws_dustmmr2_em (:,:)
REAL, ALLOCATABLE :: pws_dustmmr3_em (:,:)
REAL, ALLOCATABLE :: pws_dustmmr4_em (:,:)
REAL, ALLOCATABLE :: pws_dustmmr5_em (:,:)
REAL, ALLOCATABLE :: pws_dustmmr6_em (:,:)
LOGICAL           :: flag_dustmmr_em = .FALSE.

! Contrail forecasts
REAL, ALLOCATABLE :: pws_contrail_bot(:,:)
REAL, ALLOCATABLE :: pws_contrail_top(:,:)
LOGICAL           :: flag_contrail_bot = .FALSE.
LOGICAL           :: flag_contrail_top = .FALSE.

! Visibility-related arrays
REAL, ALLOCATABLE :: pws_bl_1p5m_vis_tot(:,:)
REAL, ALLOCATABLE :: pws_bl_1p5m_temp(:,:)
REAL, ALLOCATABLE :: pws_1p5m_vis_tot(:,:)
REAL, ALLOCATABLE :: pws_1p5m_vis_dust(:,:)
LOGICAL           :: flag_1p5m_vis_tot = .FALSE.
LOGICAL           :: flag_1p5m_vis_dust = .FALSE.

! Colson-Panofsky Turbulence-related arrays
REAL, ALLOCATABLE :: pws_panofsky_turb(:,:,:)
REAL, ALLOCATABLE :: panofsky_turb_press(:)
LOGICAL           :: flag_panofsky_turb = .FALSE.
INTEGER           :: panofsky_turb_press_levels

! Ellrodt1 Turbulence-related arrays
REAL, ALLOCATABLE :: pws_ellrodt1_turb(:,:,:)
REAL, ALLOCATABLE :: ellrodt1_turb_press(:)
LOGICAL           :: flag_ellrodt1_turb = .FALSE.
INTEGER           :: ellrodt1_turb_press_levels

! Inverse Richardson Number related arrays
REAL, ALLOCATABLE :: pws_inv_richardson(:,:,:)
REAL, ALLOCATABLE :: inv_richardson_press(:)
LOGICAL           :: flag_inv_richardson = .FALSE.
INTEGER           :: inv_richardson_press_levels

! Updraught Helicity and related velocity arrays
REAL, ALLOCATABLE :: pws_upd_helicity_5k(:,:)
REAL, ALLOCATABLE :: pws_pcd_mlev_u(:,:,:)
REAL, ALLOCATABLE :: pws_pcd_mlev_v(:,:,:)
REAL, ALLOCATABLE :: pws_pcd_mlev_w(:,:,:)
LOGICAL           :: flag_upd_helicity_5k = .FALSE.

! cut off limits to be used in tropopause calculations.
REAL, PARAMETER :: heightcut_top = 22000.0 !  arbritary limits for high
REAL, PARAMETER :: heightcut_bot = 4500.0  !  and low trop levels for search
REAL, PARAMETER :: tempcut= 243.0          !  max temp allowed for tropopause
INTEGER         :: tropo_model_top         !  model level that is close to 
                                           !  heightcut_top. It is imdi when
                                           !  model level not identified.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='PWS_DIAGS_MOD'

CONTAINS

SUBROUTINE pws_diags_alloc( &
                      icode)

USE atm_fields_bounds_mod, ONLY : tdims, udims, vdims
USE stash_array_mod, ONLY: sf, stlist, stindex, stash_levels
USE um_stashcode_mod, ONLY:                                            &
     stashcode_pws_sec, stashcode_pws_windspeed10m,                    &
     stashcode_pws_windspeedplev, stashcode_pws_thermal_advec,         &
     stashcode_pws_divergence,stashcode_pws_rel_vorticity,             &
     stashcode_pws_snow_prob, stashcode_pws_conv_cld_dep,              &
     stashcode_pws_conv_icao_base, stashcode_pws_conv_icao_top,        &
     stashcode_pws_precip_sym,stashcode_pws_max_wind_ub,               &
     stashcode_pws_max_wind_vb, stashcode_pws_max_wind_pb,             &
     stashcode_pws_max_wind_base, stashcode_pws_max_wind_top,          &
     stashcode_pws_max_wind_icao, stashcode_pws_thickness500,          &
     stashcode_pws_thickness850, stashcode_pws_tropopause_ht,          &
     stashcode_pws_tropopause_temp, stashcode_pws_tropopause_press,    &
     stashcode_pws_tropopause_icao, stashcode_pws_freezing_ht,         &
     stashcode_pws_freezing_press, stashcode_pws_freezing_icao,        &
     stashcode_pws_isotherm_ms20_ht, stashcode_pws_isotherm_ms20_press,&
     stashcode_pws_isotherm_ms20_icao, stashcode_pws_isotherm_ms70_ht, &
     stashcode_pws_isotherm_ms70_press,stashcode_pws_isotherm_ms70_icao,&
     stashcode_pws_dustconc_surf, stashcode_pws_zenithdelay,           &
     stashcode_pws_dustconc_5000, stashcode_pws_1p5m_vis_tot,          &
     stashcode_pws_1p5m_vis_dust, stashcode_pws_thickness850,          &
     stashcode_pws_cat_turb, stashcode_pws_max_cat,                    &
     stashcode_pws_max_cat_press,                                      &
     stashcode_pws_contrail_bot, stashcode_pws_contrail_top,           &
     stashcode_pws_cloudturb_pot_diag, stashcode_pws_icing_pot_diag,   &
     stashcode_pws_wafc_caturb, stashcode_pws_mtn_wave_turb,           &
     stashcode_gwd_stress_lev_u, stashcode_gwd_stress_lev_v,           &
     stashcode_pws_cloudturb_pot_diag, stashcode_pws_panofsky_turb,    &
     stashcode_pws_inv_richardson, stashcode_pws_ellrodt1_turb,        &
     stashcode_pws_upd_helicity_5k

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
USE missing_data_mod, ONLY: rmdi
USE dust_parameters_mod, ONLY: l_twobin_dust, pwsdiag_sfc_em

IMPLICIT NONE

INTEGER :: icode

INTEGER :: im_index,isl,ni,k

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PWS_DIAGS_ALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

im_index = 1

flag_windspeed_10m = sf(stashcode_pws_windspeed10m, stashcode_pws_sec)

IF (flag_windspeed_10m) THEN

  ALLOCATE(pws_wind_speed_10mb &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  pws_wind_speed_10mb(:,:) = rmdi

END IF


flag_precip_sym    = sf(stashcode_pws_precip_sym, stashcode_pws_sec)

IF (flag_precip_sym) THEN

  ALLOCATE(pws_precip_sym &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_precip_sym_ls_rain &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_precip_sym_ls_snow &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_precip_sym_conv_rain &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_precip_sym_conv_snow &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_precip_sym_t1p5m &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))

  pws_precip_sym(:,:) = 0.0 !  set to no rain and then alter if precip present.

END IF

flag_max_wind_ub   = sf(stashcode_pws_max_wind_ub,  stashcode_pws_sec)
flag_max_wind_vb   = sf(stashcode_pws_max_wind_vb,  stashcode_pws_sec)
flag_max_wind_pb   = sf(stashcode_pws_max_wind_pb,  stashcode_pws_sec)
flag_max_wind_base = sf(stashcode_pws_max_wind_base,stashcode_pws_sec)
flag_max_wind_top  = sf(stashcode_pws_max_wind_top, stashcode_pws_sec)
flag_max_wind_icao = sf(stashcode_pws_max_wind_icao,stashcode_pws_sec)

IF (flag_max_wind_ub .OR. flag_max_wind_vb .OR. flag_max_wind_pb       &
     .OR. flag_max_wind_base .OR. flag_max_wind_top                    &
     .OR. flag_max_wind_icao ) THEN
  ALLOCATE(pws_max_wind_ub                                              &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
  ALLOCATE(pws_max_wind_vb                                              &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
  ALLOCATE(pws_max_wind_pb                                              &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
  pws_max_wind_ub(:,:) = rmdi
  pws_max_wind_vb(:,:) = rmdi
  pws_max_wind_pb(:,:) = rmdi

  IF (flag_max_wind_base .OR. flag_max_wind_top) THEN
    ALLOCATE(pws_max_wind_base                                          &
            (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
    ALLOCATE(pws_max_wind_top                                           &
            (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
    pws_max_wind_base(:,:) = rmdi
    pws_max_wind_top(:,:)  = rmdi
  END IF

  IF (flag_max_wind_icao) THEN
    ALLOCATE(pws_max_wind_icao                                          &
            (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
    pws_max_wind_icao(:,:) = rmdi
  END IF
END IF

flag_windspeed_plev= sf(stashcode_pws_windspeedplev,stashcode_pws_sec)

flag_thermal_advec = sf(stashcode_pws_thermal_advec,stashcode_pws_sec)

flag_thickness_500 = sf(stashcode_pws_thickness500, stashcode_pws_sec)

IF (flag_thickness_500) THEN
  ALLOCATE(pws_geopht_500 &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_geopht_500(:,:) = rmdi
END IF

flag_thickness_850 = sf(stashcode_pws_thickness850, stashcode_pws_sec)
flag_snow_prob    = sf(stashcode_pws_snow_prob, stashcode_pws_sec)
! Snow probability calculation (Boyden Method) also needs 
! 850mb & 1000mb geopotential heights

IF (flag_thickness_850 .OR. flag_snow_prob) THEN
  ALLOCATE(pws_geopht_850 &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_geopht_850(:,:) = rmdi
END IF

IF (flag_thickness_850 .OR. flag_thickness_500 .OR. flag_snow_prob) THEN
  ALLOCATE(pws_thickness &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_geopht_1000 &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_thickness(:,:) = rmdi
  pws_geopht_1000(:,:) = rmdi
END IF

IF (flag_snow_prob) THEN
  ALLOCATE(pws_snow_prob &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_snow_prob(:,:) = rmdi
END IF

flag_conv_cld_dep = sf(stashcode_pws_conv_cld_dep,stashcode_pws_sec)
flag_conv_cld_base= sf(stashcode_pws_conv_icao_base,stashcode_pws_sec)
flag_conv_cld_top = sf(stashcode_pws_conv_icao_top,stashcode_pws_sec)
IF (flag_conv_cld_dep .OR. flag_conv_cld_base .OR. &
                           flag_conv_cld_top) THEN
  ALLOCATE(pws_conv_cld_dep &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_conv_icao_base &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_conv_icao_top &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_conv_cld_dep(:,:) = rmdi
  pws_conv_icao_base(:,:) = rmdi
  pws_conv_icao_top(:,:) = rmdi
END IF


flag_cat_turb      = sf(stashcode_pws_cat_turb, stashcode_pws_sec)
flag_max_cat       = sf(stashcode_pws_max_cat, stashcode_pws_sec)
flag_max_cat_press = sf(stashcode_pws_max_cat_press, stashcode_pws_sec)

IF (flag_cat_turb .OR. flag_max_cat .OR. flag_max_cat_press) THEN

  isl=stindex(1,16,20,im_index)

  IF (isl >  0) THEN
    ni = -stlist(10,isl)
    cat_press_levs = stash_levels(1,ni)

    ALLOCATE (cat_press(cat_press_levs))

    
    DO k = 1,cat_press_levs
      cat_press(k) = stash_levels(k+1,ni)/1000.0
      ! ***** levels are stored as integers so divide by a thousand **
    END DO
  ELSE
    cat_press_levs = 1
  END IF

! ----------Extract required pressures for Potn_vort on press ----


  ALLOCATE(pws_wind_ub_250 &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
  ALLOCATE(pws_wind_vb_250 &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
  ALLOCATE(pws_wind_ub_300 &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
  ALLOCATE(pws_wind_vb_300 &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  ALLOCATE(pws_cat_turb  &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end,cat_press_levs))

  ALLOCATE(pws_max_cat  &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  ALLOCATE(pws_max_cat_press  &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  ALLOCATE(uwind &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  ALLOCATE(vwind &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  pws_cat_turb(:,:,:)=0.0
  pws_max_cat(:,:) = 0.0
  pws_max_cat_press(:,:) = 0.0
  pws_wind_ub_250(:,:) = 0.0
  pws_wind_vb_250(:,:) = 0.0
  pws_wind_ub_300(:,:) = 0.0
  pws_wind_vb_300(:,:) = 0.0
  uwind(:,:) = 0.0
  vwind(:,:) = 0.0


END IF

flag_divergence    = sf(stashcode_pws_divergence, stashcode_pws_sec)
flag_rel_vorticity = sf(stashcode_pws_rel_vorticity,stashcode_pws_sec)
flag_mtn_wave_turb = sf(stashcode_pws_mtn_wave_turb,stashcode_pws_sec)


IF (flag_divergence .OR. flag_rel_vorticity .OR. flag_cat_turb .OR.  &
    flag_mtn_wave_turb) THEN

  ALLOCATE(pws_wind_ub_200 &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))
  ALLOCATE(pws_wind_vb_200 &
          (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  pws_wind_ub_200(:,:) = 0.0
  pws_wind_vb_200(:,:) = 0.0

END IF


flag_tropopause_ht =sf(stashcode_pws_tropopause_ht, stashcode_pws_sec)
flag_tropopause_temp =sf(stashcode_pws_tropopause_temp, stashcode_pws_sec)
flag_tropopause_press =sf(stashcode_pws_tropopause_press, stashcode_pws_sec)
flag_tropopause_icao =sf(stashcode_pws_tropopause_icao, stashcode_pws_sec)

IF (flag_tropopause_ht.OR.flag_tropopause_temp.OR.flag_tropopause_press.OR. &
                               flag_tropopause_icao   ) THEN
  ALLOCATE(pws_tropopause_ht &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_tropopause_temp &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_tropopause_press &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_tropopause_icao &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_tropopause_ht(:,:) = rmdi
  pws_tropopause_press(:,:) = rmdi
  pws_tropopause_temp(:,:) = rmdi
  pws_tropopause_icao(:,:) = rmdi
END IF

flag_icing_pot = sf(stashcode_pws_icing_pot_diag,stashcode_pws_sec)

IF (flag_icing_pot) THEN

  isl=stindex(1,stashcode_pws_icing_pot_diag,20,im_index)

  IF (isl >  0) THEN
    ni = -stlist(10,isl)
    icing_pot_press_levs = stash_levels(1,ni)
    ALLOCATE(icing_pot_press(icing_pot_press_levs))

    DO k = 1,icing_pot_press_levs
      icing_pot_press(k) = stash_levels(k+1,ni)/1000.0
      ! ***** levels are stored as integers so divide by a thousand **
      ! next multiply by 100 to convert to pascals as stash uses hPa.
      icing_pot_press(k) = icing_pot_press(k)*100.0
    END DO
  ELSE
    ALLOCATE(icing_pot_press(icing_pot_press_levs))
    icing_pot_press_levs = 1
  END IF

  ALLOCATE( pws_icing_pot_diag                                         &
          ( tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end,         &
                                           icing_pot_press_levs))
  pws_icing_pot_diag(:,:,:) = rmdi

END IF


flag_cloudturb_pot = sf(stashcode_pws_cloudturb_pot_diag,stashcode_pws_sec) 

isl=stindex(1,stashcode_pws_cloudturb_pot_diag,20,im_index)

IF (flag_cloudturb_pot) THEN

  IF (isl >  0) THEN
    ni = -stlist(10,isl)
    pot_press_levs = stash_levels(1,ni)
    ALLOCATE(pot_press(pot_press_levs))

    DO k = 1,pot_press_levs
      pot_press(k) = stash_levels(k+1,ni)/1000.0
      ! ***** levels are stored as integers so divide by a thousand **
    END DO
  ELSE
    pot_press_levs = 1
  END IF

  ALLOCATE( pws_cloudturb_pot_diag                                         &
        (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end, pot_press_levs))
  pws_cloudturb_pot_diag(:,:,:) = 0.0

END IF


flag_mtn_wave_turb = sf(stashcode_pws_mtn_wave_turb,stashcode_pws_sec)

IF (flag_mtn_wave_turb) THEN

  ALLOCATE( pws_mtn_wave_turb  &
        (udims%i_start:udims%i_end, vdims%j_start:vdims%j_end))

  pws_mtn_wave_turb(:,:) = 0.0 

END IF


flag_contrail_bot     = sf(stashcode_pws_contrail_bot,  stashcode_pws_sec)
flag_contrail_top     = sf(stashcode_pws_contrail_top,  stashcode_pws_sec)

IF (flag_contrail_bot .OR. flag_contrail_top) THEN

  ALLOCATE(pws_contrail_bot &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  ALLOCATE(pws_contrail_top &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_contrail_bot(:,:) = 0.0
  pws_contrail_top(:,:) = 0.0

END IF

flag_wafc_cat = sf(stashcode_pws_wafc_caturb,stashcode_pws_sec)

IF (flag_wafc_cat) THEN

  isl=stindex(1,stashcode_pws_wafc_caturb,stashcode_pws_sec,im_index)

  IF (isl >  0) THEN
    ni = -stlist(10,isl)
    wafc_cat_press_levs = stash_levels(1,ni)
    ALLOCATE(wafc_cat_press(wafc_cat_press_levs))

    DO k = 1,wafc_cat_press_levs
      wafc_cat_press(k) = stash_levels(k+1,ni)/1000.0
      ! ***** levels are stored as integers so divide by a thousand **
    END DO
  ELSE
    wafc_cat_press_levs = 1
    ALLOCATE(wafc_cat_press(wafc_cat_press_levs))
  END IF

  ALLOCATE( pws_wafc_cat_diag                                         &
          ( udims%i_start:udims%i_end, vdims%j_start:vdims%j_end,         &
                                           wafc_cat_press_levs))
  pws_wafc_cat_diag(:,:,:) = rmdi
END IF

IF (flag_wafc_cat .OR. flag_mtn_wave_turb) THEN
  ALLOCATE( pws_gwd_stress_lev_u  &
        (udims%i_start:udims%i_end, udims%j_start:udims%j_end, 0:udims%k_end))
  ALLOCATE( pws_gwd_stress_lev_v  &
        (vdims%i_start:vdims%i_end, vdims%j_start:vdims%j_end, 0:vdims%k_end))
 
  pws_gwd_stress_lev_u(:,:,:) = rmdi
  pws_gwd_stress_lev_v(:,:,:) = rmdi

END IF



flag_freezing_ht      = sf(stashcode_pws_freezing_ht,  stashcode_pws_sec)
flag_freezing_press   = sf(stashcode_pws_freezing_press,  stashcode_pws_sec)
flag_freezing_icao    = sf(stashcode_pws_freezing_icao,  stashcode_pws_sec)
flag_isotherm_ms20_ht = sf(stashcode_pws_isotherm_ms20_ht,                   &
                                               stashcode_pws_sec)
flag_isotherm_ms20_press  = sf(stashcode_pws_isotherm_ms20_press,            &
                                               stashcode_pws_sec)
flag_isotherm_ms20_icao   = sf(stashcode_pws_isotherm_ms20_icao,             &
                                               stashcode_pws_sec)
flag_isotherm_ms70_ht     = sf(stashcode_pws_isotherm_ms70_ht,               &
                                               stashcode_pws_sec)
flag_isotherm_ms70_press  = sf(stashcode_pws_isotherm_ms70_press,            &
                                               stashcode_pws_sec)
flag_isotherm_ms70_icao   = sf(stashcode_pws_isotherm_ms70_icao,             &
                                               stashcode_pws_sec)

IF (flag_freezing_ht         .OR. flag_freezing_press      .OR.             &
    flag_freezing_icao       .OR. flag_isotherm_ms20_ht    .OR.             &
    flag_isotherm_ms20_press .OR. flag_isotherm_ms20_icao  .OR.             &
    flag_isotherm_ms70_ht    .OR. flag_isotherm_ms70_press .OR.             &
    flag_isotherm_ms70_icao ) THEN

  ALLOCATE(pws_freezing_ht &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_freezing_ht(:,:) = rmdi
  ALLOCATE(pws_freezing_press &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_freezing_press(:,:) = rmdi
  ALLOCATE(pws_freezing_icao &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_freezing_icao(:,:) = rmdi

END IF


flag_dustconc_surf = sf(stashcode_pws_dustconc_surf, stashcode_pws_sec)
flag_dustconc_5000 = sf(stashcode_pws_dustconc_5000, stashcode_pws_sec)

IF (flag_dustconc_surf .OR. flag_dustconc_5000) THEN
  ALLOCATE(pws_dustconc_tot &
          (tdims%i_start:tdims%i_end, &
           tdims%j_start:tdims%j_end, &
                       1:tdims%k_end))      ! 1:model_levels
END IF

IF (flag_dustconc_surf) THEN
  ALLOCATE(pws_dustconc_surf &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_dustconc_surf(:,:) = rmdi
END IF

IF (flag_dustconc_5000) THEN
  ALLOCATE(pws_dustconc_5000 &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_dustconc_5000(:,:) = rmdi
END IF

flag_1p5m_vis_tot = sf(stashcode_pws_1p5m_vis_tot, stashcode_pws_sec)

IF (flag_1p5m_vis_tot) THEN
  ALLOCATE(pws_1p5m_vis_tot &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_1p5m_vis_tot(:,:) = rmdi
END IF

flag_1p5m_vis_dust = sf(stashcode_pws_1p5m_vis_dust, stashcode_pws_sec)

IF (flag_1p5m_vis_dust) THEN
  ALLOCATE(pws_1p5m_vis_dust &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_1p5m_vis_dust(:,:) = rmdi
END IF

IF (flag_1p5m_vis_tot .OR. flag_1p5m_vis_dust) THEN
  ALLOCATE(pws_bl_1p5m_vis_tot &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_bl_1p5m_vis_tot(:,:) = rmdi
  ALLOCATE(pws_bl_1p5m_temp &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_bl_1p5m_temp(:,:) = rmdi
END IF 

! set up the store for the dust MMRs, if requried:
IF (flag_1p5m_vis_tot .OR. flag_1p5m_vis_dust .OR. flag_dustconc_surf) THEN
  IF (pwsdiag_sfc_em > 0.0) THEN
    flag_dustmmr_em = .TRUE.
    IF (l_twobin_dust) THEN
       ALLOCATE(pws_dustmmr1_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
       ALLOCATE(pws_dustmmr2_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))       
    ELSE
       ALLOCATE(pws_dustmmr1_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
       ALLOCATE(pws_dustmmr2_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))       
       ALLOCATE(pws_dustmmr3_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
       ALLOCATE(pws_dustmmr4_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))       
       ALLOCATE(pws_dustmmr5_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
       ALLOCATE(pws_dustmmr6_em &
            (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))       
    END IF
  END IF
END IF


flag_zenithdelay = sf(stashcode_pws_zenithdelay, stashcode_pws_sec)

IF (flag_zenithdelay) THEN

  ALLOCATE(pws_zen_tot_delay(tdims%i_start:tdims%i_end, &
                             tdims%j_start:tdims%j_end))
  pws_zen_tot_delay(:,:) = 0.0

END IF


flag_panofsky_turb = sf(stashcode_pws_panofsky_turb, stashcode_pws_sec)
IF (flag_panofsky_turb) THEN

  isl = stindex(1, stashcode_pws_panofsky_turb, stashcode_pws_sec, im_index)
  IF (isl > 0) THEN
    ni = - stlist(10, isl)
    panofsky_turb_press_levels = stash_levels(1, ni)
    ALLOCATE(panofsky_turb_press(panofsky_turb_press_levels))
    DO k = 1, panofsky_turb_press_levels
      panofsky_turb_press(k) = stash_levels(k+1,ni)/1000.0
      ! Levels are integers so divide by one thousand. Then multiply by 100
      ! to convert hPa to Pa.
      panofsky_turb_press(k) = panofsky_turb_press(k)  * 100.0
    END DO
  ELSE
    panofsky_turb_press_levels = 1
    ALLOCATE(panofsky_turb_press(panofsky_turb_press_levels))
  END IF
  
  ALLOCATE(pws_panofsky_turb(tdims%i_start:tdims%i_end,         &
                             tdims%j_start:tdims%j_end,         &
                             panofsky_turb_press_levels))
  pws_panofsky_turb(:,:,:) = rmdi

END IF

flag_ellrodt1_turb = sf(stashcode_pws_ellrodt1_turb, stashcode_pws_sec)
IF (flag_ellrodt1_turb) THEN

  isl = stindex(1, stashcode_pws_ellrodt1_turb, stashcode_pws_sec, im_index)
  IF (isl > 0) THEN
    ni = - stlist(10, isl)
    ellrodt1_turb_press_levels = stash_levels(1, ni)
    ALLOCATE(ellrodt1_turb_press(ellrodt1_turb_press_levels))
    DO k = 1, ellrodt1_turb_press_levels
      ellrodt1_turb_press(k) = stash_levels(k+1,ni)/1000.0
      ! Levels are integers so divide by one thousand. Then multiply by 100
      ! to convert hPa to Pa.
      ellrodt1_turb_press(k) = ellrodt1_turb_press(k)  * 100.0
    END DO
  ELSE
    ellrodt1_turb_press_levels = 1
    ALLOCATE(ellrodt1_turb_press(ellrodt1_turb_press_levels))
  END IF
  
  ALLOCATE(pws_ellrodt1_turb(udims%i_start:udims%i_end,         &
                             vdims%j_start:vdims%j_end,         &
                             ellrodt1_turb_press_levels))
  pws_ellrodt1_turb(:,:,:) = rmdi

END IF

flag_inv_richardson = sf(stashcode_pws_inv_richardson, stashcode_pws_sec)
IF (flag_inv_richardson) THEN

  isl = stindex(1, stashcode_pws_inv_richardson, stashcode_pws_sec, im_index)
  IF (isl > 0) THEN
    ni = - stlist(10, isl)
    inv_richardson_press_levels = stash_levels(1, ni)
    ALLOCATE(inv_richardson_press(inv_richardson_press_levels))
    DO k = 1, inv_richardson_press_levels
      inv_richardson_press(k) = stash_levels(k+1,ni)/1000.0
      ! Levels are integers so divide by one thousand. Then multiply by 100
      ! to convert hPa to Pa.
      inv_richardson_press(k) = inv_richardson_press(k)  * 100.0
    END DO
  ELSE
    inv_richardson_press_levels = 1
    ALLOCATE(inv_richardson_press(inv_richardson_press_levels))
  END IF
  
  ALLOCATE(pws_inv_richardson(tdims%i_start:tdims%i_end,         &
                              tdims%j_start:tdims%j_end,         &
                              inv_richardson_press_levels))
  pws_inv_richardson(:,:,:) = rmdi

END IF

flag_upd_helicity_5k = sf(stashcode_pws_upd_helicity_5k, stashcode_pws_sec)

IF (flag_upd_helicity_5k) THEN
  ALLOCATE(pws_pcd_mlev_u &
          (tdims%i_start:tdims%i_end, &
           tdims%j_start:tdims%j_end, &
                       1:tdims%k_end))      ! 1:model_levels
  ALLOCATE(pws_pcd_mlev_v &
          (tdims%i_start:tdims%i_end, &
           tdims%j_start:tdims%j_end, &
                       1:tdims%k_end))      ! 1:model_levels
  ALLOCATE(pws_pcd_mlev_w &
          (tdims%i_start:tdims%i_end, &
           tdims%j_start:tdims%j_end, &
                       1:tdims%k_end))      ! 1:model_levels

  ALLOCATE(pws_upd_helicity_5k &
          (tdims%i_start:tdims%i_end, tdims%j_start:tdims%j_end))
  pws_upd_helicity_5k(:,:) = rmdi
END IF


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE pws_diags_alloc

!-------------------------------------------------------------

SUBROUTINE pws_diags_dealloc()

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PWS_DIAGS_DEALLOC'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (ALLOCATED(pws_wind_speed_10mb))    DEALLOCATE(pws_wind_speed_10mb)
IF (ALLOCATED(pws_precip_sym))         DEALLOCATE(pws_precip_sym)
IF (ALLOCATED(pws_precip_sym_ls_rain)) DEALLOCATE(pws_precip_sym_ls_rain)
IF (ALLOCATED(pws_precip_sym_ls_snow)) DEALLOCATE(pws_precip_sym_ls_snow)
IF (ALLOCATED(pws_precip_sym_conv_rain)) DEALLOCATE(pws_precip_sym_conv_rain)
IF (ALLOCATED(pws_precip_sym_conv_snow)) DEALLOCATE(pws_precip_sym_conv_snow)
IF (ALLOCATED(pws_precip_sym_t1p5m))   DEALLOCATE(pws_precip_sym_t1p5m)
IF (ALLOCATED(pws_max_wind_ub))        DEALLOCATE(pws_max_wind_ub)
IF (ALLOCATED(pws_max_wind_vb))        DEALLOCATE(pws_max_wind_vb)
IF (ALLOCATED(pws_max_wind_pb))        DEALLOCATE(pws_max_wind_pb)
IF (ALLOCATED(pws_max_wind_base))      DEALLOCATE(pws_max_wind_base)
IF (ALLOCATED(pws_max_wind_top))       DEALLOCATE(pws_max_wind_top)
IF (ALLOCATED(pws_max_wind_icao))      DEALLOCATE(pws_max_wind_icao)
IF (ALLOCATED(pws_thickness))          DEALLOCATE(pws_thickness)
IF (ALLOCATED(pws_geopht_1000))        DEALLOCATE(pws_geopht_1000)
IF (ALLOCATED(pws_geopht_500))         DEALLOCATE(pws_geopht_500)
IF (ALLOCATED(pws_geopht_850))         DEALLOCATE(pws_geopht_850)
IF (ALLOCATED(pws_snow_prob))          DEALLOCATE(pws_snow_prob)
IF (ALLOCATED(pws_conv_cld_dep))       DEALLOCATE(pws_conv_cld_dep)
IF (ALLOCATED(pws_conv_icao_base))     DEALLOCATE(pws_conv_icao_base)
IF (ALLOCATED(pws_conv_icao_top))      DEALLOCATE(pws_conv_icao_top)
IF (ALLOCATED(uwind))                  DEALLOCATE(uwind)
IF (ALLOCATED(vwind))                  DEALLOCATE(vwind)
IF (ALLOCATED(pws_wind_ub_300))        DEALLOCATE(pws_wind_ub_300)
IF (ALLOCATED(pws_wind_vb_300))        DEALLOCATE(pws_wind_vb_300)
IF (ALLOCATED(pws_wind_ub_250))        DEALLOCATE(pws_wind_ub_250)
IF (ALLOCATED(pws_wind_vb_250))        DEALLOCATE(pws_wind_vb_250)
IF (ALLOCATED(pws_wind_ub_200))        DEALLOCATE(pws_wind_ub_200)
IF (ALLOCATED(pws_wind_vb_200))        DEALLOCATE(pws_wind_vb_200)
IF (ALLOCATED(pws_max_cat_press))      DEALLOCATE(pws_max_cat_press)
IF (ALLOCATED(pws_max_cat))            DEALLOCATE(pws_max_cat)
IF (ALLOCATED(cat_press))              DEALLOCATE(cat_press)
IF (ALLOCATED(pws_cat_turb))           DEALLOCATE(pws_cat_turb)
IF (ALLOCATED(pws_contrail_top))       DEALLOCATE(pws_contrail_top)
IF (ALLOCATED(pws_contrail_bot))       DEALLOCATE(pws_contrail_bot)
IF (ALLOCATED(pot_press))              DEALLOCATE(pot_press)
IF (ALLOCATED(pws_cloudturb_pot_diag)) DEALLOCATE(pws_cloudturb_pot_diag)
IF (ALLOCATED(wafc_cat_press))         DEALLOCATE(wafc_cat_press)
IF (ALLOCATED(pws_wafc_cat_diag))      DEALLOCATE(pws_wafc_cat_diag)
IF (ALLOCATED(pws_gwd_stress_lev_u))   DEALLOCATE(pws_gwd_stress_lev_u)
IF (ALLOCATED(pws_gwd_stress_lev_v))   DEALLOCATE(pws_gwd_stress_lev_v)
IF (ALLOCATED(pws_mtn_wave_turb))      DEALLOCATE(pws_mtn_wave_turb)
IF (ALLOCATED(pws_tropopause_icao))    DEALLOCATE(pws_tropopause_icao)
IF (ALLOCATED(pws_tropopause_press))   DEALLOCATE(pws_tropopause_press)
IF (ALLOCATED(pws_tropopause_temp))    DEALLOCATE(pws_tropopause_temp)
IF (ALLOCATED(pws_tropopause_ht))      DEALLOCATE(pws_tropopause_ht)
IF (ALLOCATED(pws_freezing_ht))        DEALLOCATE(pws_freezing_ht)
IF (ALLOCATED(pws_freezing_press))     DEALLOCATE(pws_freezing_press)
IF (ALLOCATED(pws_freezing_icao))      DEALLOCATE(pws_freezing_icao)
IF (ALLOCATED(icing_pot_press))        DEALLOCATE(icing_pot_press)
IF (ALLOCATED(pws_icing_pot_diag))     DEALLOCATE(pws_icing_pot_diag)
IF (ALLOCATED(pws_dustconc_surf))      DEALLOCATE(pws_dustconc_surf)
IF (ALLOCATED(pws_dustconc_5000))      DEALLOCATE(pws_dustconc_5000)
IF (ALLOCATED(pws_dustconc_tot))       DEALLOCATE(pws_dustconc_tot)
IF (ALLOCATED(pws_1p5m_vis_tot))       DEALLOCATE(pws_1p5m_vis_tot)
IF (ALLOCATED(pws_1p5m_vis_dust))      DEALLOCATE(pws_1p5m_vis_dust)
IF (ALLOCATED(pws_bl_1p5m_vis_tot))    DEALLOCATE(pws_bl_1p5m_vis_tot)
IF (ALLOCATED(pws_dustmmr6_em))        DEALLOCATE(pws_dustmmr6_em)
IF (ALLOCATED(pws_dustmmr5_em))        DEALLOCATE(pws_dustmmr5_em)
IF (ALLOCATED(pws_dustmmr4_em))        DEALLOCATE(pws_dustmmr4_em)
IF (ALLOCATED(pws_dustmmr3_em))        DEALLOCATE(pws_dustmmr3_em)
IF (ALLOCATED(pws_dustmmr2_em))        DEALLOCATE(pws_dustmmr2_em)
IF (ALLOCATED(pws_dustmmr1_em))        DEALLOCATE(pws_dustmmr1_em)
IF (ALLOCATED(pws_bl_1p5m_temp))       DEALLOCATE(pws_bl_1p5m_temp)
IF (ALLOCATED(pws_zen_tot_delay))      DEALLOCATE(pws_zen_tot_delay)
IF (ALLOCATED(pws_panofsky_turb))      DEALLOCATE(pws_panofsky_turb)
IF (ALLOCATED(panofsky_turb_press))    DEALLOCATE(panofsky_turb_press)
IF (ALLOCATED(pws_ellrodt1_turb))      DEALLOCATE(pws_ellrodt1_turb)
IF (ALLOCATED(ellrodt1_turb_press))    DEALLOCATE(ellrodt1_turb_press)
IF (ALLOCATED(pws_inv_richardson))     DEALLOCATE(pws_inv_richardson)
IF (ALLOCATED(inv_richardson_press))   DEALLOCATE(inv_richardson_press)
IF (ALLOCATED(pws_upd_helicity_5k))    DEALLOCATE(pws_upd_helicity_5k)
IF (ALLOCATED(pws_pcd_mlev_w))         DEALLOCATE(pws_pcd_mlev_w)
IF (ALLOCATED(pws_pcd_mlev_v))         DEALLOCATE(pws_pcd_mlev_v)
IF (ALLOCATED(pws_pcd_mlev_u))         DEALLOCATE(pws_pcd_mlev_u)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE pws_diags_dealloc


END MODULE pws_diags_mod
