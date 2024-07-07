
! *********************************COPYRIGHT*********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt which you should
! have received as part of this distribution.
! *********************************COPYRIGHT*********************************
!
!  Routine to run the SCM.
!
! Subroutine Interface:
SUBROUTINE scm_main                                                           &
  ( vert_lev, nfor, l_ts_log, ntrop, sec_day, land_points, nsprog             &
  , co2_dim_len, co2_dim_row, cloud_levels                                    &
  , tr_levels, tr_vars, tr_ukca, st_levels, sm_levels, bl_levels              &
  , ozone_levels, ntiles, nice, nice_use, l_netcdf_obs )

  ! Physics modules
  !---------------------------------------------------------------------------
USE jules_snow_mod,    ONLY: nsmax
USE ancil_info,        ONLY: l_soil_point, l_lice_point
USE rad_input_mod,     ONLY: l_radiation,                                   &
                             l_use_dust,                                    &
                             l_use_soot_direct, l_use_soot_indirect,        &
                             l_use_bmass_direct, l_use_bmass_indirect,      &
                             l_use_ocff_direct, l_use_ocff_indirect,        &
                             l_use_nitrate_direct, l_use_nitrate_indirect,  &
                             l_rad_deg,                                     &
                             a_sw_radstep_prog, a_sw_radstep_diag,          &
                             lrad_ccrad,                                    &
                             co2_mmr,                                       &
                             it_rad1
USE set_rad_steps_mod, ONLY: set_a_radstep,                                 &
                             l_rad_step_prog, l_rad_step_diag
USE tuning_segments_mod, ONLY: a_lw_segments, a_lw_seg_size, a_sw_seg_size

USE sw_control_struct
USE lw_control_struct
USE ukca_radaer_struct_mod, ONLY: ukca_radaer_struct
USE mphys_inputs_mod, ONLY: l_psd, l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,   &
                            l_use_ocff_autoconv,  l_use_soot_autoconv,      &
                            l_use_bmass_autoconv, l_use_nitrate_autoconv
USE stochastic_physics_run_mod, ONLY: rhcrit_max, rhcrit_min                &
  , m_ci, m_ci_max, m_ci_min, g0_rp, par_mezcla, lambda_min_rp              &
  , ricrit_rp, a_ent_1_rp, g1_rp

USE cloud_inputs_mod, ONLY: rhcrit, i_cld_area, i_cld_vn
USE pc2_constants_mod, ONLY: acf_off, acf_brooks, i_cld_pc2
USE river_inputs_mod, ONLY: river_vel, river_mcoef, l_inland, river_step,   &
                            i_river_vn
USE mphys_bypass_mod
USE g_wave_input_mod, ONLY: Gsharp, fbcd,l_smooth,l_nonhydro,l_dynbeta      &
  , l_gw_heating, l_use_ussp
USE turb_diff_mod, ONLY: qlimit, l_qpos, l_leonard_term
USE turb_diff_ctl_mod, ONLY:                                                &
    visc_m, visc_h, rneutml_sq, max_diff, delta_smag, shear
USE arcl_mod,  ONLY: npd_arcl_compnts, npd_arcl_species, ip_arcl_sulp,      &
                     ip_arcl_dust, ip_arcl_sulp, ip_arcl_dust, ip_arcl_sslt,&
                     ip_arcl_blck, ip_arcl_biom, ip_arcl_ocff, ip_arcl_dlta

USE ukca_option_mod,  ONLY: &
    l_ukca_chem,            &
    l_ukca_set_trace_gases, &
    l_ukca_strat,           &
    l_ukca_strattrop,       &
    l_ukca_prescribech4,    &
    l_ukca_aie1,            &
    l_ukca_aie2,            &
    l_ukca_radaer,          &
    l_ukca_radaer_sustrat
                           
USE glomap_clim_option_mod, ONLY: &
    l_glomap_clim_aie1,           &
    l_glomap_clim_aie2,           &
    l_glomap_clim_radaer,         &
    l_glomap_clim_radaer_sustrat

USE jules_surface_types_mod, ONLY: npft
USE nlsizes_namelist_mod, ONLY: model_levels
USE def_easyaerosol, ONLY: t_easyaerosol_rad, t_easyaerosol_cdnc

! SCM modules
!---------------------------------------------------------------------------
USE netcdf
USE netcdf_obs, ONLY: get_netcdf_obs
USE mcc_data, ONLY: get_mcc, alloc_mcc, dealloc_mcc, mcc_o3, mcc_th,        &
  mcc_th_t, mcc_q, mcc_rh_p
USE global_scmop
USE scm_utils, ONLY: scm_message, alloc_common_scm, dealloc_common_scm,     &
  time_info, scm_timestep_count, daycount, scm_timestep, nlnd, z_th, z_rh,  &
  l_emcorr_opt, nprimvars, time_string, resdump
USE restart_dump_mod, ONLY: restart_dump
USE dumpinit_mod, ONLY: dumpinit
USE s_main_force, timestep_nml=>timestep
USE scm_cntl_mod
USE s_scmop_mod, ONLY:                                                      &
    t_inst, t_avg, t_max, t_min, t_acc, t_div, t_mult, t_acc_div            &
  , t_acc_mult, t_const, only_radsteps                                      &
  , d_sl, d_soilt, d_bl, d_wet, d_all, d_soilm, d_tile, d_vis, d_point      &
  , d_allxtra, d_land, d_cloud, d_mlsnow, default_streams                   &
  , SCMDiag_gen, SCMDiag_rad, SCMDiag_bl, SCMDiag_surf, SCMDiag_land        &
  , SCMDiag_sea, SCMDiag_lsp, SCMDiag_conv, SCMDiag_lscld, SCMDiag_pc2      &
  , SCMDiag_forc, SCMDiag_incs, SCMDiag_gwd, SCMDiag_MLSnow
USE scmoutput_mod, ONLY: scmoutput

! Constants modules
!---------------------------------------------------------------------------
USE visbty_constants_mod, ONLY: n_vis_thresh, vis_thresh
USE planet_constants_mod, ONLY: cp, omega, two_omega
USE water_constants_mod,  ONLY: lc
USE conversions_mod,      ONLY: pi_over_180, pi, isec_per_day, rsec_per_day
USE errormessagelength_mod, ONLY: errormessagelength

! UM Modules
!---------------------------------------------------------------------------
USE vertnamelist_mod, ONLY:                                                 &
    first_constant_r_rho_level, z_top_of_model, eta_theta, eta_rho          &
  , vertlevs

USE Control_Max_Sizes, ONLY:  max_121_rows                                  &
  , max_updiff_levels, max_bl_levels, max_req_thpv_levs, max_sponge_width   &
  , max_look

! Trig arrays (requires dynamics A12_2A code)
USE trignometric_mod, ONLY:                                                 &
    cos_theta_latitude, sec_theta_latitude, FV_cos_theta_latitude           &
  , sin_theta_latitude, cos_theta_longitude, sin_theta_longitude            &
  , true_longitude, true_latitude
USE cderived_mod, ONLY: delta_lambda, delta_phi
! In the UM, Coriolis array is allocated in set_trigs.f90, which isn't
! called in the SCM.  So need to allocate it in scm_main, so use it:
USE dyn_coriolis_mod, ONLY:                                                 &
    f3_at_u

! Model level heights from centre of Earth  (requires dynamics A12_2A code)
USE level_heights_mod, ONLY:                                                &
    eta_theta_levels, eta_rho_levels, r_theta_levels, r_rho_levels,         &
    z_top_theta, r_layer_centres, r_layer_boundaries, d_layer

USE timestep_mod, ONLY: timestep_number, timestep, recip_timestep

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


USE ereport_mod, ONLY: ereport
USE UM_ParVars, ONLY: change_decomposition, at_extremity,                   &
  halo_i, halo_j, offx, offy
USE decomp_db, ONLY: decompose
USE decomp_params, ONLY: decomp_smexe
USE trophgt1_mod, ONLY: z_min_trop, z_max_trop

USE dynamics_input_mod, ONLY:                                &
    NumCycles, l_lbc_old
USE lam_config_inputs_mod, ONLY: n_rims_to_do
USE dynamics_testing_mod, ONLY:                              &
    L_dry
USE dust_parameters_mod, ONLY: l_dust, l_dust_diag
USE run_aerosol_mod, ONLY:  &
     l_soot, l_biomass, l_ocff, l_nitrate, so2_high_level
USE eng_corr_inputs_mod, ONLY: lflux_reset
USE murk_inputs_mod, ONLY: l_murk, l_murk_bdry, l_murk_rad
USE cv_run_mod, ONLY: l_ccrad, l_3d_cca, l_mr_conv,                    &
    l_conv_prog_group_1, l_conv_prog_group_2, l_conv_prog_group_3,     &
    l_conv_prog_precip
USE cosp_input_mod, ONLY: l_cosp
USE jules_radiation_mod, ONLY: l_snow_albedo
USE jules_vegetation_mod, ONLY: l_triffid, l_phenol
USE nlstcall_mod, ONLY: lcal360

USE atm_fields_bounds_mod, ONLY: pdims, tdims

USE umPrintMgr, ONLY: umPrint, umMessage, newline

USE model_domain_mod, ONLY: l_regular, model_type, mt_single_column

USE ls_acf_brooks_mod, ONLY: ls_acf_brooks
USE pc2_homog_plus_turb_mod, ONLY: pc2_homog_plus_turb
USE pc2_initiation_ctl_mod, ONLY: pc2_initiation_ctl

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_wat_new     => qsat_wat,                       &
                    l_new_qsat_cntl !Currently defaults to FALSE

USE gen_phys_inputs_mod, ONLY: l_mr_physics

IMPLICIT NONE

!
! Description:
!   Subroutine that runs the single column model.  This replaces the previous
!   version of scm_main that was the top level of the model.
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
! Code description:
!   FORTRAN 90 This code is written to UM programming standards version 8.2

! Arguments with Intent In. ie: Input variables.

INTEGER, INTENT(IN) :: nfor        ! Number terms for observational forcing
INTEGER, INTENT(IN) :: ntrop       ! Max number of levels in the troposphere
INTEGER, INTENT(IN) :: co2_dim_len
INTEGER, INTENT(IN) :: co2_dim_row
INTEGER, INTENT(IN) :: cloud_levels
INTEGER, INTENT(IN) :: tr_levels
INTEGER, INTENT(IN) :: tr_vars
INTEGER, INTENT(IN) :: tr_ukca
INTEGER, INTENT(IN) :: st_levels
INTEGER, INTENT(IN) :: sm_levels
INTEGER, INTENT(IN) :: bl_levels
INTEGER, INTENT(IN) :: ozone_levels
INTEGER, INTENT(IN) :: ntiles
INTEGER, INTENT(IN) :: nice
INTEGER, INTENT(IN) :: nice_use
INTEGER, INTENT(IN) :: sec_day
INTEGER, INTENT(IN) :: land_points ! No of land points - can be 0
INTEGER, INTENT(IN) :: nsprog      ! No. of single level prognostics


CHARACTER(LEN=200), INTENT(IN) :: vert_lev
CHARACTER(LEN=200)             :: sw_spec_file
CHARACTER(LEN=200)             :: lw_spec_file
CHARACTER(LEN=200)             :: FILE

LOGICAL,          INTENT(IN) :: &
  l_ts_log                      &! Output timestep information
, l_netcdf_obs                   ! Using netcdf driver file for forcing

!=====================================================================



! Local variables

INTEGER  :: istatus

INTEGER, PARAMETER :: height_gen_method = 2

INTEGER  :: icode
INTEGER  :: n_cca_lev
INTEGER  :: errorstatus
REAL     :: a2out(row_length,rows,model_levels)

! Arrays for PC2 initiation control arg list.
REAL ::                                     &
  t_work   (row_length, rows, model_levels) &
, q_work   (row_length, rows, model_levels) &
, qcl_work (row_length, rows, model_levels) &
, qcf_work (row_length, rows, model_levels) &
, cf_work  (row_length, rows, model_levels) &
, cfl_work (row_length, rows, model_levels) &
, cff_work (row_length, rows, model_levels)


!-----------------------------------------------------------------------------
! Dummy/Fixed variables
!-----------------------------------------------------------------------------
! These either have not impact on the SCM, or must be a set value
! for use with the SCM.

INTEGER :: CycleNo = 1  ! Define CycleNo default

INTEGER, PARAMETER ::  &
  n_proc   = 1         &! Total number of processors
, n_procx  = 1         &! Number of processors in longitude
, n_procy  = 1          ! Number of processors in latitude

!------------------------------------------------------------


REAL, PARAMETER ::         &
  lat_rot_NP        = 0.0  &! Diagnostic variable
, long_rot_NP       = 0.0   ! Diagnostic variable

REAL, PARAMETER ::         &! STASH space for
  stashwork1        = 0.0  &! Section 1  (sw)
, stashwork2        = 0.0  &! Section 2  (lw)
, stashwork3        = 0.0  &! Section 3  (b layer)
, stashwork4        = 0.0  &! Section 4  (ls precipitation)
, stashwork5        = 0.0  &! Section 5  (convection)
, stashwork6        = 0.0  &! Section 6  (gravity wave drag)
, stashwork8        = 0.0  &! Section 8  (hydrology)
, stashwork9        = 0.0  &! Section 9  (LS Cloud)
, stashwork14       = 0.0  &! Section 14 (energy correction)
, stashwork19       = 0.0  &! Section 19 (Veg)
, stashwork21       = 0.0  &! Section 21 (electric)
, stashwork26       = 0.0   ! Section 26 (River Routing)


! Duplication of variables needed to pass UKCA gases to radiation scheme.
! This option will not work in the SCM.
INTEGER, PARAMETER         :: ngrgas     = 8
INTEGER                    :: grgas_addr(ngrgas) = -1


INTEGER ::            &
  land_pts_trif  = 1  &! To dimension land fields
, npft_trif      = 1  &! To dimension pft  fields
, co2_dim_lev    = 1   ! 3-D CO2 field passed to NI_rad_ctl


! Declare pointers for dummy arrays
REAL, POINTER ::           &
  lambda_p                 &! VarRes hor. co-ordinate information
, lambda_u                 &! VarRes hor. co-ordinate information
, dlambda_p                &! VarRes hor. co-ordinate information
, wt_lambda_p              &! VarRes hor. co-ordinate information
, wt_lambda_u               ! VarRes hor. co-ordinate information


INTEGER, POINTER ::        &
  g_p_field                &! Size of global atmos field
, g_r_field                &! Size of global river field
, a_steps_since_riv        &! Number of physics timsteps since
                            ! last call to river routing
, river_row_length         &! Local  river row length
, river_rows               &! Local  river rows
, global_river_row_length  &! Global river row length
, global_river_rows         ! Global river rows


REAL, POINTER :: &
  inlandout_atm   (:)             &! Inland basin flow (kg/m2/s)
, tot_surf_runoff (:)             &! Accumulated runoff over river
, tot_sub_runoff  (:)             &! routing timestep (Kg/m2/s)
, xpa  (:)                        &! Atmosphere TP long coordinates
, xua  (:)                        &! Atmosphere U  long coordinates
, xva  (:)                        &! Atmosphere V  long coordinates
, ypa  (:)                        &! Atmosphere TP lat  coordinates
, yua  (:)                        &! Atmosphere U  lat  coordinates
, yva  (:)                         ! Atmosphere V  lat  coordinates


REAL, POINTER :: &
  acc_lake_evap(:,:) &! Accumulated lake evap over river routing timestep
, r_area(:,:)        &! Accumulated areas file
, slope(:,:)         &
, flowobs1(:,:)      &! Initialisation for flows
, r_inext(:,:)       &! x-coordinate of downstream grid point
, r_jnext(:,:)       &! y-coordinate of downstream grid point
, r_land(:,:)        &! Land/river/sea
, substore (:,:)     &! Routing sub_surface store (mm)
, surfstore(:,:)     &! Routing surface store (mm)
, flowin (:,:)       &! Rurface lateral inflow (mm)
, bflowin(:,:)       &! Rub-surface lateral inflow (mm)
, trivdir(:,:)       &! River direction file
, trivseq(:,:)       &! River sequence file
, twatstor (:,:)     &! Water storage
, soil_clay(:,:)     &! Fields for mineral dust
, soil_silt(:,:)     &! source flux calculations
, soil_sand(:,:)     &
, dust_mrel1(:,:)    &
, dust_mrel2(:,:)    &
, dust_mrel3(:,:)    &
, dust_mrel4(:,:)    &
, dust_mrel5(:,:)    &
, dust_mrel6(:,:)    &
, phi_p (:,:)        &! VarRes horizontal co-ordinate information
, phi_v (:,:)        &! VarRes horizontal co-ordinate information
, dphi_p (:,:)       &! VarRes horizontal co-ordinate information
, wt_phi_p (:,:)     &! VarRes horizontal co-ordinate information
, wt_phi_v (:,:)      ! VarRes horizontal co-ordinate information


REAL, POINTER :: &
  dust_div1(:,:,:)     &! Tracer variables
, dust_div2(:,:,:)     &
, dust_div3(:,:,:)     &
, dust_div4(:,:,:)     &
, dust_div5(:,:,:)     &
, dust_div6(:,:,:)

! UKCA_RADAER structure, for consistent call to Atmos_Physics1
TYPE (ukca_radaer_struct) :: ukca_radaer

! EasyAerosol structures, for consistent call to Atmos_Physics1
TYPE (t_easyaerosol_rad)  :: easyaerosol_sw
TYPE (t_easyaerosol_rad)  :: easyaerosol_lw
TYPE (t_easyaerosol_cdnc) :: easyaerosol_cdnc

! Declare dummy target arrays

REAL, TARGET ::                          &
  rdum0                                  &
, rdum1(row_length)                      &
, rdum2(row_length+1)                    &
, rdum3(row_length, rows)                &
, rdum4(row_length, rows, tdims%k_start:tdims%k_end)  &
, rdum5(row_length, rows, bl_levels)     &
, rdum6(row_length, rows, bl_levels-1)   &
, rdum7(row_length, rows, nice)          &
, dum_land(land_points)

INTEGER, TARGET :: idum0

REAL ::                                  &
  frac_control    (land_points,ntype)    &
, es_space_interp (4,row_length, rows)   &
, o3_trop_level   (row_length,rows)      &
, o3_trop_height  (row_length,rows)      &
, t_trop_level    (row_length,rows)      &
, t_trop_height   (row_length,rows)


! CDNC from UKCA-MODE, for consistent call to Atmos_Physics1
! This is currently a dummy array but could be used
! in future if UKCA is to be used with the SCM

! Declare dummy array and dimensions

INTEGER, PARAMETER :: cdnc_dim1 = 1
INTEGER, PARAMETER :: cdnc_dim2 = 1
INTEGER, PARAMETER :: cdnc_dim3 = 1

REAL :: cdnc_ukca_dummy( cdnc_dim1, cdnc_dim2, cdnc_dim3 )

!-----------------------------------------------------------------------------
! &INOBSFOR information
!-----------------------------------------------------------------------------
! We need to keep the variables that are read in from &INOBSFOR and reset
! them to their correct values at the start of each timestep so use dummy
! equivalent variables.

REAL ::                                    &
  tls     (row_length,rows,model_levels)   &
, uls     (row_length,rows,model_levels)   &
, vls     (row_length,rows,model_levels)   &
, wls     (row_length,rows,0:model_levels) &
, qls     (row_length,rows,model_levels)   &
, qcl_inc (row_length,rows,model_levels)   &
, qcf_inc (row_length,rows,model_levels)



!-----------------------------------------------------------------------------
! Time information
!-----------------------------------------------------------------------------
INTEGER ::     &
  ihour        &! Current hour
, imin         &! Current min
, isec = 0      ! Current sec, not implemented

!-----------------------------------------------------------------------------
! Initial day of the year and initial time (in seconds) in that day
!-----------------------------------------------------------------------------
INTEGER ::      &
  dayno_init    &! Initial day
, time_initi     ! Initial time in seconds (integer)

REAL ::         &
  time_init      ! Initial time in seconds (real)


!-----------------------------------------------------------------------------
! Primary Model Variables plus P_star (UMDP No1)
!-----------------------------------------------------------------------------

INTEGER ::                &
  iccb(row_length, rows)  &! Convective cloud base
, icct(row_length, rows)   ! Convective cloud top

REAL ::                                 &
  canopy_gb(land_points)                &! Canopy water content(kg/m2)
, cca(row_length,rows,model_levels)     &! Convective cloud amount
, cca_dp(row_length,rows,model_levels)  &! Deep Conv cloud amount
, cca_md(row_length,rows,model_levels)  &! Mid-level Conv cloud amount
, cca_sh(row_length,rows,model_levels)  &! Shallow Conv cloud amount
, q(row_length,rows,tdims%k_start:tdims%k_end)       & 
                                                ! Specific humidity (kg/kg)
, qcf(row_length,rows,tdims%k_start:tdims%k_end)     &
                                                ! Cloud ice content (kg/kg)
, qcl(row_length,rows,tdims%k_start:tdims%k_end)     &
                                                ! Cloud water content(kg/kg)
, qcf2(row_length,rows,tdims%k_start:tdims%k_end)    &
                                                ! 2nd ice water content (kg/kg)
, qrain(row_length,rows,tdims%k_start:tdims%k_end)   &
                                                ! Rain water content (kg/kg)
, qgraup(row_length,rows,tdims%k_start:tdims%k_end)  &
                                                ! Graupel water content(kg/kg)
, mix_v(row_length,rows,tdims%k_start:tdims%k_end)   &! Vapour mixing ratio
, mix_cf(row_length,rows,tdims%k_start:tdims%k_end)  &! Cloud ice mixing ratio
, mix_cl(row_length,rows,tdims%k_start:tdims%k_end)  &! Cloud water mixing ratio
, mix_cf2(row_length,rows,tdims%k_start:tdims%k_end) &! 2nd ice mixing ratio
, mix_rain(row_length,rows,tdims%k_start:tdims%k_end)&! Rain mixing ratio
, mix_graup(row_length,rows,tdims%k_start:tdims%k_end) &! Graupel mixing ratio
, bl_w_var(row_length, rows, model_levels) ! Turbulence for mixed-phase

REAL ::                                            &
  exner_rho_levels(row_length,rows,model_levels+1) &
, exner_theta_levels(row_length,rows,tdims%k_start:tdims%k_end) &
, exner_prime(row_length,rows,model_levels)         ! Increment to exner

REAL ::                                         &
  dtheta_dr_term(row_length,rows,tdims%k_start:tdims%k_end) &
, rho(row_length,rows,model_levels)             &
, rho_n(row_length,rows,model_levels)           &
, p(row_length,rows,model_levels+1)             &! Pressure on rho levels
, p_theta_levels(row_length,rows,tdims%k_start:tdims%k_end)  &
                                                 ! Pressure on theta levels
, rp(row_length,rows,model_levels+1)            &! 1/p on rho levels
, rp_theta(row_length,rows,model_levels)         ! 1/p on theta levels

REAL ::                          &
  p_star(row_length,rows)        &! Surface pressure
, smc(land_points)               &! Soil moisture content (Kg/m^2)
, smcl(land_points,sm_levels)    &! Soil moisture content in layers (Kg/m^2)
, sthf(land_points,sm_levels)    &! Frozen soil moisture content of each
                                  ! layer as a fraction of saturation.
, sthu(land_points,sm_levels)     ! Unfrozen soil moisture content of each
                                  ! layer as a fraction of saturation.

!
!---------------------------------------------------------------------
!     Water
!---------------------------------------------------------------------
REAL ::                                 &
  snodep(row_length,rows)               &! Snow depth (Kg/m^2)
, t(row_length,rows,model_levels)       &! Temperature (K)  
, t_deep_soil(land_points,st_levels)    &! Deep soil temperatures (K)
                                         ! top level not included, =surface
, tsi(row_length,rows)                  &! Temperature of sea-ice
, tstar(row_length,rows)                &! Surface temperature (K)
, u(row_length,rows,model_levels)       &! Zonal wind (m/s)
, v(row_length,rows,model_levels)       &! Meridional wind (m/s)
, w(row_length,rows,0:model_levels)     &! Vertical velocity (m/s)
, z0msea(row_length,rows)               &! Sea surface roughness length
, chloro_sea(row_length,rows)           &! Sea nr surface chlorophyll
, w_adv(row_length,rows,0:model_levels)  ! Advective w component of wind

REAL :: &
  ozone_tracer(row_length,rows,tdims%k_start:tdims%k_end)  ! Ozone tracer

REAL ::                         &
  deep_flag(row_length,rows)    &! Indicator of how long since last
                                 ! deep convective event.
, past_precip(row_length,rows)  &! Decayed memory of past total
                                 ! convective precipitation
, past_conv_ht(row_length,rows)  ! Memory of past height of convection

! Convection prognostics
REAL :: conv_prog_1(row_length, rows, tdims%k_start:tdims%k_end)
REAL :: conv_prog_2(row_length, rows, tdims%k_start:tdims%k_end)
REAL :: conv_prog_3(row_length, rows, tdims%k_start:tdims%k_end)
REAL :: conv_prog_precip(row_length, rows,tdims%k_start:tdims%k_end)

! convective cold pool (CCP) prognostics
REAL :: ux_ccp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: uy_ccp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: um_ccp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: g_ccp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: h_ccp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: riso_ccp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
REAL :: rdir_ccp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL :: totalppn(row_length, rows)
REAL ::                                            &
  theta_star     (row_length, rows, tdims%k_start:tdims%k_end)  &
, qcl_star       (row_length, rows, tdims%k_start:tdims%k_end)  &
, qcf_star       (row_length, rows, tdims%k_start:tdims%k_end)  &
, qcf2_star      (row_length, rows, tdims%k_start:tdims%k_end)  &
, qrain_star     (row_length, rows, tdims%k_start:tdims%k_end)  &
, qgraup_star    (row_length, rows, tdims%k_start:tdims%k_end)  &
, mix_v_star     (row_length, rows, model_levels)  &
, mix_cl_star    (row_length, rows, model_levels)  &
, mix_cf_star    (row_length, rows, model_levels)  &
, mix_cf2_star   (row_length, rows, model_levels)  &
, mix_rain_star  (row_length, rows, model_levels)  &
, mix_graup_star (row_length, rows, model_levels)  &
, cf_star        (row_length, rows, tdims%k_start:tdims%k_end)  &
, cfl_star       (row_length, rows, tdims%k_start:tdims%k_end)  &
, cff_star       (row_length, rows, tdims%k_start:tdims%k_end)  &
, uinc_geo       (row_length, rows, model_levels)  &
, vinc_geo       (row_length, rows, model_levels)  &
, theta_eg       (row_length, rows, tdims%k_start:tdims%k_end)

!-----------------------------------------------------------------------------
! CCRad Prognostics
!-----------------------------------------------------------------------------
INTEGER :: &
  lcbase (row_length, rows) ! Model level of lowest convective cloud base
                            ! in profile(theta level) passed to radiation

REAL :: &
  ccw (row_length, rows, model_levels)
                            ! Convective Cloud Water profile passed to
                            ! radiation scheme (theta levels). (kg/kg)
                            ! NOTE: May be different to ccw produced by
                            !       Convection Scheme it Radiative
                            !       convective cloud decay is used.


!-----------------------------------------------------------------------------
! PC2 cloud and condensation scheme
!-----------------------------------------------------------------------------
REAL ::                                        &
  q_forcing   (row_length, rows, model_levels) &
, qcl_forcing (row_length, rows, model_levels) &
, t_forcing   (row_length, rows, model_levels) &
, p_forcing   (row_length, rows, model_levels) &
, q_earliest  (row_length, rows, model_levels) &
, qcl_earliest(row_length, rows, model_levels) &
, t_earliest  (row_length, rows, model_levels)


!-----------------------------------------------------------------------------
! Water
!-----------------------------------------------------------------------------
REAL ::                     &
  ls_rain(row_length, rows) &! Large scale rainfall rate (Kg/m^2)
, ls_snow(row_length, rows) &! Large scale snowfall rate (Kg/m^2/s)
, ls_graup(row_length, rows) ! Large scale graupel fall rate (Kg/m^2/s)

REAL ::                                                               &
  ls_rainfrac(land_points)
                         ! Rain fraction array on land points to be
                         ! passed to atmos_physics2

!---------------------------------------------------------------------
! JULES2 snow scheme prognostics
!---------------------------------------------------------------------
REAL ::                             &
  snowdepth(    land_points,ntiles) &! Snow depth on ground on tiles (m)
, rho_snow_grnd(land_points,ntiles) &! Snowpack bulk density (kg/m3)
! , nsnow(        land_points,ntiles) &! No. of snow layers on ground on tiles
                                       ! NOTE: this is converted to an integer

  , ds(     land_points,ntiles,nsmax) &! Snow layer thickness (m)
  , sice(   land_points,ntiles,nsmax) &! Snow layer ice  mass on tiles (Kg/m2)
  , sliq(   land_points,ntiles,nsmax)  ! Snow layer liq. mass on tiles (Kg/m2)

REAL ::                                &
  rho_snow  (land_points,ntiles,nsmax) &! Snow layer densities (kg/m3)
, tsnowlayer(land_points,ntiles,nsmax) &! Snow layer temperature (K)
, rgrainl   (land_points,ntiles,nsmax)  ! Snow layer grain size on tiles
                                        ! (microns)

!-----------------------------------------------------------------------------
! Large scale statistical forcing
!-----------------------------------------------------------------------------
! Random generator variables

INTEGER ::            &
  iseed                ! Seed for random number generator

INTEGER ::            &
  dayno_wint           ! Day number relative to winter solstice

REAL ::                                &
  ad(row_length, rows,model_levels-1)  &! Term a of equation 2.22
                                        !  for dewpoint depression
, alfada(row_length,rows)              &! Amplitude and mean of seasonal
, alfadb(row_length,rows)              &!  variation of tuning factor
, at(row_length,rows,model_levels-1)   &! Term a of equation 2.22
                                        !  for dewpoint depression
, atime, btime                         &! Constants for calculating annual
                                        ! cycle
, avn(row_length,rows,model_levels-1)  &! Term a of equation 2.22
, aw(row_length,rows,ntrop-1)          &!  for horiz. and vert. vel.
, cdbar(row_length,rows,model_levels)  &! Mean and SD of random variable
, cdsd(row_length,rows,model_levels)   &!  for dew pt. depression
, ctbar(row_length,rows,model_levels)  &! Mean and SD of random variable
, ctsd(row_length,rows,model_levels)   &!  for temp.
, cvnbar(row_length,rows,model_levels) &! Mean and SD of random variable
, cvnsd(row_length,rows,model_levels)  &!  for velocity VN
, cwbar(row_length,rows,ntrop)         &! Mean and SD of random variable
, cwsd(row_length,rows,ntrop)          &!  for vertical velocity
, dbar(row_length,rows,model_levels)   &! Mean dewpoint depression at
                                        !  daycount days from winter
                                        !  solstice (K)
, dbara(row_length,rows,model_levels)  &! Amplitude and mean of seasonal
, dbarb(row_length,rows,model_levels)  &!  variation of mean dew pt.
                                        ! depression (K)
, ddash(row_length,rows,model_levels)  &! Dew pt depression correction (K)
, deltan(row_length,rows)              &! Radius of area (m)
, dgrada(row_length,rows,model_levels) &! Amplitude and mean of seasonal
, dgradb(row_length,rows,model_levels) &!  variation of dew pt. depression
                                        !  gradient (K/km)
, dsd(row_length,rows,model_levels)    &! SD dewpoint depression at daycount
                                        !  days from winter solstice (K)
, pa(row_length,rows,model_levels+1)   &! Amplitude and mean of seasonal
, pb(row_length,rows,model_levels+1)   &!  variation of pressure (Pa)
, px(row_length,rows,ntrop)            &! Reciprocal log functions for
, py(row_length,rows,ntrop-1)          &!  calc. of vert. advection
                                        !  used in eqns 2.12 and 2.13
, qr(row_length,rows,model_levels,2)   &! Randomly sampled specific
                                        !  humidity (kg/kg)
, tbar(row_length,rows,model_levels)   &! Mean temperature at daycount days
                                        !  from winter solstice (K)
, tbara(row_length,rows,model_levels)  &! Amplitude and mean of seasonal
, tbarb(row_length,rows,model_levels)  &!  variation of temp. (K)
, tdash(row_length,rows,model_levels)  &! Temp correction (K)
, tgrada(row_length,rows,model_levels) &! Amplitude and mean of seasonal
, tgradb(row_length,rows,model_levels) &!  variation of temp. gradient
                                        !  (k/Km)
, tr(row_length,rows,model_levels,2)   &! InOut Randomly sampled temp. (K)
, tsd(row_length,rows,model_levels)    &! SD temp. at daycount days
                                        ! from winter solstice (K)
, tsda(row_length,rows,model_levels)   &! Amplitude and mean of seasonal
, tsdb(row_length,rows,model_levels)   &!  variation of SD of temp. (K)
, vnbar(row_length,rows,model_levels)  &! Mean velocity VN at daycount days
                                        !  from winter solstice
, vnbara(row_length,rows,model_levels) &! Amplitude and mean of seasonal
, vnbarb(row_length,rows,model_levels) &!  variation of velocity VN (m/s)
, vnr(row_length,rows,model_levels,2)  &! InOut Randomly sampled horizontal
                                        !       velocity (m/s)
, vnsd(row_length,rows,model_levels)   &! SD velocity VN at daycount days
                                        !  from winter solstice (m/s)
, vnsda(row_length,rows,model_levels)  &! Amplitude and mean of seasonal
, vnsdb(row_length,rows,model_levels)  &!  variation of SD of velocity VN
                                        !  (m/s)
, vpbar(row_length,rows,model_levels)  &! Mean velocity VP at daycount days
                                        !  from winter solstice
, vpbara(row_length,rows,model_levels) &! Amplitude and mean of seasonal
, vpbarb(row_length,rows,model_levels) &!  variation of velocity VP (m/s)
, vpr(row_length,rows,model_levels,2)  &! InOut Randomly sampled horizontal
                                        !  velocity (m/s)
, wbar(row_length,rows,ntrop)          &! Mean vertical velocity at daycount
                                        !  days from winter solstice (m/s)
, wbara(row_length,rows,ntrop)         &! Amplitude and mean of seasonal
, wbarb(row_length,rows,ntrop)         &!  variation of SD of vert. vel.
                                        !  (mb or HPa/s)
, wr(row_length,rows,ntrop,2)          &! InOut Randomly sampled vertical
                                        !  velocity (mb/s)
, wsd(row_length,rows,ntrop)           &! SD vertical velocity at daycount
                                        !  days from winter solstice (m/s)
, wsda(row_length,rows,ntrop)          &! Amplitude and mean of seasonal
, wsdb(row_length,rows,ntrop)           !  variation of SD of vert. vel.
                                        !  (mb/s) roughness length (m)

!-----------------------------------------------------------------------------
! Large scale observational forcing
!-----------------------------------------------------------------------------
! Variables for diagnostic output for observational forcing

REAL ::                                        &
  dap1(row_length,rows,36,model_levels)        &! Instantaneous profiles
, dap2(row_length,rows,36,model_levels)        &! Mean profiles
, dap3(row_length,rows,36,nfor-1,model_levels) &! Mean profiles - timeseries
, dab1(row_length,rows,44)                     &! Instantaneous budgets
, dab2(row_length,rows,44)                     &! Mean budgets
, dab3(row_length,rows,44,nfor-1)              &! Mean budgets - timeseries
, deltap(row_length,rows,model_levels)         &! Layer thickness (Pa)
, factor_rhokh(row_length,rows)                 ! Used to specify surface
                                                ! flux from observation



!-----------------------------------------------------------------------------
! Variables enabling diagnostic output
!-----------------------------------------------------------------------------

INTEGER, PARAMETER ::   &
  nSCMdpkgs = 15         ! No of diags packages

! Note that if nSCMDpkgs is changed it should also be changed
! in PC2_PRES and ATMSTEP2

LOGICAL ::              &
  l_scmdiags(nscmdpkgs)  ! Logicals for diagnostics packages

LOGICAL ::              &
  kill_scmdiags(nscmdpkgs)
                         ! Disable diagnostics for geostrophic
                         ! initialisation call

INTEGER ::              &
  site(row_length*rows)  ! SSFM site WMO number

! Arrays to store initial fields for calculating total increments and
! total increment fields
REAL ::                                      &
  t_start(row_length,rows,model_levels)      &
, u_start(row_length,rows,model_levels)      &
, v_start(row_length,rows,model_levels)      &
, w_start(row_length,rows,0:model_levels)    &
, q_start(row_length,rows,model_levels)      &
, qcl_start(row_length,rows,model_levels)    &
, qcf_start(row_length,rows,model_levels)    &
, cf_start(row_length,rows,model_levels)     &
, cfl_start(row_length,rows,model_levels)    &
, cff_start(row_length,rows,model_levels)

REAL ::                                      &
  t_totalinc(row_length,rows,model_levels)   &
, u_totalinc(row_length,rows,model_levels)   &
, v_totalinc(row_length,rows,model_levels)   &
, w_totalinc(row_length,rows,0:model_levels) &
, q_totalinc(row_length,rows,model_levels)   &
, qcl_totalinc(row_length,rows,model_levels) &
, qcf_totalinc(row_length,rows,model_levels) &
, cf_totalinc(row_length,rows,model_levels)  &
, cfl_totalinc(row_length,rows,model_levels) &
, cff_totalinc(row_length,rows,model_levels)


! Qpos arrays
!---------------------------------------------------------------------
! Hold incs to q to prevent -ve values.
REAL ::                                      &
  dq_qpos_for (row_length,rows,model_levels) &
, dq_qpos     (row_length,rows,model_levels) &
, dqcl_qpos   (row_length,rows,model_levels) &
, dqcf_qpos   (row_length,rows,model_levels)


!-----------------------------------------------------------------------------
! Extra diagnostics boundary layer code
!-----------------------------------------------------------------------------

REAL ::                             &
  resp_w_ft(row_length, rows,npft)   ! Out Wood maintenance respiration
                                     !     (Kg C m^-2 s^-1).

!-----------------------------------------------------------------------------
! Clouds
!-----------------------------------------------------------------------------

REAL ::                                    &
  cf(row_length,rows,tdims%k_start:tdims%k_end)    &! Layer cloud amount
                                            ! (decimal fraction)
, cfl(row_length,rows,tdims%k_start:tdims%k_end)   &! liquid layer cloud amount
, cff(row_length,rows,tdims%k_start:tdims%k_end)   &! frozen layer cloud amount
, area_cloud_fraction(row_length,rows,model_levels)&
                                            ! area cloud fraction
, ccwpin(row_length,rows)                  &! Condensed water path (Kg/m^2)
, rhts(row_length,rows,model_levels)       &! Initial relative humidity
                                            ! wrt TL
, qtts(row_length,rows,model_levels)       &! Initial q_T
, tlts(row_length,rows,model_levels)       &! Initial TL
, ptts(row_length,rows,model_levels)       &! Initial p on theta levs
, tl(row_length,rows)                      &! Liquid water temperature (K)
, qsl_tl(row_length,rows)                   ! Qsat wrt liquid water
                                            ! at temp TL


!-----------------------------------------------------------------------------
! Radiation
!-----------------------------------------------------------------------------
INTEGER ::         &
  daynumber        &! Day in the year (default=1)
, year              ! Year


!-----------------------------------------------------------------------------
! Electric
!-----------------------------------------------------------------------------

REAL :: flash_pot(row_length, rows, tdims%k_start:tdims%k_end)

!-----------------------------------------------------------------------------
! Used in calculation of min/max tropopause
!-----------------------------------------------------------------------------
REAL :: r_ref_rho(model_levels)  ! height of model rho levels above surface.


!-----------------------------------------------------------------------------
! Loop Counters and limits
!-----------------------------------------------------------------------------
INTEGER ::       &
  daysteps       &! Number of timestep in a day
, full_daysteps  &! Number of timesteps in a full day
, nstepsin       &! Number of steps in final day
, total_nsteps   &! Total number of steps in run
, i, j, j2, k    &! General loop counters, array dumpmean
, itype          &
, m1              ! No. of dumps

!-----------------------------------------------------------------------------
! Miscellaneous
!-----------------------------------------------------------------------------
CHARACTER(LEN=errormessagelength) :: &
  cmessage, iomessage         ! Error message if Icode >0

CHARACTER(LEN=*), PARAMETER :: RoutineName = 'SCM_MAIN'

INTEGER ::        &
  error_code      &
, day             &
, yearno          &
, time_sec         ! Actual time of day in secs.

REAL ::                        &
  modug(row_length,rows)       &! Magnitude of geostrophic wind (m/s)
, f_coriolis(row_length,rows)  &! 2*omega*sin(latitude)
, rccb(row_length,rows)        &! Real val. cloud base geoint only
, rcct(row_length,rows)        &! Real cloud top geoint only
, maxinc                        ! Maximum wind increment from geoinit

REAL ::                                     &
  sw_incs(row_length,rows,0:model_levels+1) &
, lw_incs(row_length,rows,0:model_levels)   &
, dirpar_incs(row_length,rows)              &! PAR flux variable
, t1_sd(row_length,rows)                    &! Set to zero initially
, q1_sd(row_length,rows)                    &! Set to zero initially
, sw_tile_rts(row_length*rows,ntiles)       &! Surface net SW on land tiles
, rhokh(row_length,rows,bl_levels)
  ! Exchange coeffs for moisture.
  !       surface:out of SF_EXCH
  !               contains only RHOKH.
  ! above surface:out of KMKH
  !               contains GAMMA(1)*RHOKH(,1)*RDZ(,1)

REAL ::                           &
  gridbox_area_m(row_length,rows)  ! Gridbox area in m^2

INTEGER ::                        &
  land_ice_points                 &! No of ice land points
, soil_points                     &! No of soil points
, land_index(row_length*rows)     &! Should be defined by land_points
, land_ice_index(row_length*rows) &! but land_points has not yet been
, soil_index(row_length*rows)      ! given a value

! Land point specific surface soil parameters
REAL ::                &
  b_exp  (land_points) &! Clapp-Hornberger exponent
, hcap   (land_points) &! Soil heat capacity
, hcon   (land_points) &! Soil thermal conductivity
, satcon (land_points) &! Saturated hydrological conductivity
, sathh  (land_points) &! Saturated soil water suction
, v_sat  (land_points) &! Volumetric SMC at saturation
, v_wilt (land_points) &! Volumetric SMC at wilting point
, v_crit (land_points) &! Volumetric SMC at critical point
, z0m_soil(land_points) ! Bare soil surface roughness length


REAL ::                                               &
  energy_correction                                   &
, tracer_ukca(row_length,rows,tr_levels,tr_ukca)

REAL ::                                               &
  biogenic (row_length,rows,model_levels)

! Local arrays for using the aerosol climatology for NWP

! Internal model switches
LOGICAL :: &
  l_use_arcl(npd_arcl_species)

! Internal array of mass-mixing ratios
REAL, ALLOCATABLE :: arcl(:,:,:,:)


INTEGER ::        &
  n_arcl_species  &! Number of requested species within the climatology
, n_arcl_compnts   ! Corresponding number of requested components

INTEGER :: &
  i_arcl_compnts(npd_arcl_compnts) ! Array index of each component

REAL :: &
  unscl_dry_rho(row_length,rows,model_levels) ! unscaled dry density

REAL ::                                       &
  rhokm   (row_length,rows,0:bl_levels-1)     &
, ch_term (row_length,rows,model_levels-1)

REAL ::                                       &
  photosynth_act_rad(row_length, rows)        &! Net downward shortwave
                                               ! radiation in band 1 (w/m2).
, rad_hr(row_length, rows, 2, bl_levels)      &! BL radiative heating rates
, micro_tends(row_length, rows, 2, bl_levels) &! BL Microphys tendencies
, dOLR(row_length, rows)                      &! TOA - surface upward LW
, sw_tile(row_length*rows, ntiles)            &! Surface net SW on land tiles
, cos_zenith_angle(row_length, rows)

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
REAL :: zlcl_mixed(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!---------------------------------------------------------------------
! Rate of change in forcing increments (1/s)
!---------------------------------------------------------------------
REAL ::                                            &
  ch_t_inc (row_length,rows,nfor-1,model_levels)   &! Temperature
, ch_u_inc (row_length,rows,nfor-1,model_levels)   &! Zonal wind
, ch_v_inc (row_length,rows,nfor-1,model_levels)   &! Meridional wind
, ch_w_inc (row_length,rows,nfor-1,0:model_levels) &! Vertical wind
, ch_q_star(row_length,rows,nfor-1,model_levels)   &! Specific Humidity
, ch_flux_e(row_length,rows,nfor-1)                &! Latent heat flux
, ch_flux_h(row_length,rows,nfor-1)                &! Sensible heat flux
, ch_tstar_forcing(row_length,rows,nfor-1)         &! Surface Temperature
, ch_ustar_forcing(row_length,rows,nfor-1)          ! Friction velocity

REAL ::                                            &
  ch_ug(row_length,rows,nfor-1,model_levels)       &! Geostrophic u-wind
, ch_vg(row_length,rows,nfor-1,model_levels)        ! Geostrophic v-wind

REAL ::                                            &
  ch_q_bg (row_length,rows,nfor-1,model_levels)    &! Observed q
, ch_t_bg (row_length,rows,nfor-1,model_levels)    &! Observed T
, ch_u_bg (row_length,rows,nfor-1,model_levels)    &! Observed u-wind
, ch_v_bg (row_length,rows,nfor-1,model_levels)    &! Observed v-wind
, ch_w_bg (row_length,rows,nfor-1,0:model_levels)   ! Observed w-velocity

!---------------------------------------------------------------------
! Forcings passed to physics routines
!---------------------------------------------------------------------

REAL ::                          &
  flux_e_scm (row_length, rows)  &! Surface latent heat fluxes
, flux_h_scm (row_length, rows)  &! Surface sensible heat fluxes
, ustar_in   (row_length, rows)   ! Surface friction velocity

REAL ::                                       &
  ug_scm (row_length,rows,model_levels)       &! Geostrophic u-wind (m/s)
, vg_scm (row_length,rows,model_levels)        ! Geostrophic v-wind (m/s)

! Increments/timestep due to large-scale horizontal and vertical advection
REAL ::                                       &
  t_inc_scm  (row_length,rows,model_levels)   &! Temperature inc (K)
, u_inc_scm  (row_length,rows,model_levels)   &! Zonal wind inc (m/s)
, v_inc_scm  (row_length,rows,model_levels)   &! Merid wind inc (m/s)
, w_inc_scm  (row_length,rows,0:model_levels) &! Vert  wind inc (m/s)
, q_star_scm (row_length,rows,tdims%k_start:tdims%k_end) 
                                               ! Spec. humid. inc (kg/kg)


! Background field profiles (for relaxation)
REAL ::                                    &
  t_bg_scm(row_length,rows,model_levels)   &! Temperature (K)
, u_bg_scm(row_length,rows,model_levels)   &! Zonal wind (m/s)
, v_bg_scm(row_length,rows,model_levels)   &! Merid wind (m/s)
, w_bg_scm(row_length,rows,0:model_levels) &! Vert  wind (m/s)
, q_bg_scm(row_length,rows,model_levels)    ! Spec. humid (kg/kg)

LOGICAL ::         &
  l_mr_pc2         &! True if using mixing ratios in PC2 section
, l_calc_exner     &
, l_calc_rho

LOGICAL ::         &
  rad_mask(row_length,rows)

! Land based variables
REAL ::                    &
  fexp       (land_points) &! Decay factor in Sat. conductivity in water
                            ! table layer.
, gamtot     (land_points) &! Integrated complete Gamma function.
, ti_mean    (land_points) &! Mean topographic index.
, ti_sig     (land_points) &! St dev. of topographic index. in water
                            ! table layer
, fsat       (land_points) &! Surface saturation fraction.
, fwetl      (land_points) &! Wetland fraction.
, zw         (land_points) &! Water table depth (m).
, sthzw      (land_points) &! Soil moist fract. in deep-zw layer.
, a_fsat     (land_points) &! Fitting parameter for Fsat in LSH model
, c_fsat     (land_points) &! Fitting parameter for Fsat in LSH model
, a_fwet     (land_points) &! Fitting parameter for Fwet in LSH model
, c_fwet     (land_points) &! Fitting parameter for Fwet in LSH model
, catch_snow (land_points, ntiles) &! Snow capacity for NLT tile (kg/m2)
, snow_grnd  (land_points, ntiles)  ! Snow below canopy (kg/m2).


CHARACTER(LEN=300) :: sdum0, sdum1, sdum2

CHARACTER(LEN=10)  :: nfor_str, mx_nfor_str


INTEGER :: asteps_since_triffid
INTEGER :: previous_time(7)      ! Time information for current timestep

! Variables for COSP, for consistency in calls to atmos_physics1/2
! Convective rainfall
REAL :: cosp_crain_3d(row_length,rows,model_levels)
! Convective snowfall
REAL :: cosp_csnow_3d(row_length,rows,model_levels)

! Ancillary fields and fields needed to be kept from timestep to timestep
! Assuming here that nice=nice_use
REAL ::                                &
  ice_fract_n(row_length, rows, nice)  &! Ice categories
, ice_thick_n(row_length, rows, nice)  &
, ti_n       (row_length, rows)        & ! Ice temperature as grid box mean
, ice_k_n(row_length, rows, nice)      & ! Effective cond of sea ice
, tstar_sice_n(row_length, rows, nice_use) & ! Surface ice T
, snodep_sice_n(row_length, rows, nice_use)  ! Snow depth on ice

! Add new variables from extra args in calls to physics routines
REAL :: &
  rhcpt(row_length, rows, model_levels)

! Work array for scmoutput for ug and vg
REAL :: &
  geo_diag(row_length, rows, model_levels)

! TKE based turbulence scheme
REAL ::                                            &
  e_trb(row_length, rows, model_levels)            &
!                   ! TKE defined on theta levels K-1
  , tsq_trb(row_length, rows, model_levels)          &
!                   ! Self covariance of liquid potential temperature
!                   ! (thetal'**2) defined on theta levels K-1
  , qsq_trb(row_length, rows, model_levels)          &
!                   ! Self covariance of total water
!                   ! (qw'**2) defined on theta levels K-1
  , cov_trb(row_length, rows, model_levels)          &
!                   ! Correlation between thetal and qw
!                   ! (thetal'qw') defined on theta levels K-1
  , zhpar_shcu(row_length, rows)
!                   ! Height of mixed layer used to evaluate
!                   ! the non-gradient buoyancy flux


INTEGER :: &
  ichgf    &! No. of timesteps between change in observational forcing
, ilscnt    ! Count for observational forcing

INTEGER :: &
  nout(19)  ! Output units to receive printout of initial
            ! data by PRINT_INITDATA

REAL ::    &
  ti (row_length,rows,model_levels) ! Initial temperature profile, (K)

! Dr Hook
!==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Assign variables normally for MPP configs
CALL decompose(1,1,0,0,model_levels,noSmallHalos=.TRUE.)
CALL Change_Decomposition(decomp_smexe)
at_extremity      = .FALSE.

! Assign Pointers to dummy target arrays
g_p_field               => idum0
g_r_field               => idum0
a_steps_since_riv       => idum0
river_row_length        => idum0
river_rows              => idum0
global_river_row_length => idum0
global_river_rows       => idum0
lambda_p                => rdum0
phi_p                   => rdum3
lambda_u                => rdum0
phi_v                   => rdum3
dlambda_p               => rdum0
dphi_p                  => rdum3
wt_lambda_p             => rdum0
wt_lambda_u             => rdum0
wt_phi_p                => rdum3
wt_phi_v                => rdum3
dust_div1               => rdum4
dust_div2               => rdum4
dust_div3               => rdum4
dust_div4               => rdum4
dust_div5               => rdum4
dust_div6               => rdum4
acc_lake_evap           => rdum3
r_area                  => rdum3
slope                   => rdum3
flowobs1                => rdum3
r_inext                 => rdum3
r_jnext                 => rdum3
r_land                  => rdum3
substore                => rdum3
surfstore               => rdum3
flowin                  => rdum3
bflowin                 => rdum3
trivdir                 => rdum3
trivseq                 => rdum3
twatstor                => rdum3
soil_clay               => rdum3
soil_silt               => rdum3
soil_sand               => rdum3
dust_mrel1              => rdum3
dust_mrel2              => rdum3
dust_mrel3              => rdum3
dust_mrel4              => rdum3
dust_mrel5              => rdum3
dust_mrel6              => rdum3
inlandout_atm           => dum_land
tot_surf_runoff         => dum_land
tot_sub_runoff          => dum_land
xpa                     => rdum2
xua                     => rdum2
xva                     => rdum2
ypa                     => rdum1
yua                     => rdum1
yva                     => rdum2

! Initialise Dummy/Hardwired values

ice_fract_n   = 0.0    ! Was uninitialised before
ice_thick_n   = 0.0    ! Was uninitialised before
ti_n          = 273.0  ! Ice temperature as grid box mean
ice_k_n       = 41.8   ! 2*kappai/de
tstar_sice_n  = 273.0
snodep_sice_n = 0.0
model_type    = mt_single_column  ! Specifies model type
n_rims_to_do  = 1                 ! Size of LAM rim zone
l_regular     = .TRUE.            ! True if NOT variable resolution
l_lbc_old     = .TRUE.            ! True if old lbcs (not active scm)
numcycles     = 1                 ! cycling dynamics-physics to
                                  ! be compatible with UM

i_river_vn         = 1
fexp               = 1.0
gamtot             = 1.0
ti_mean            = 1.0
ti_sig             = 1.0
fsat               = 1.0
fwetl              = 1.0
zw                 = 1.0
sthzw              = 1.0
a_fsat             = 1.0
c_fsat             = 1.0
a_fwet             = 1.0
c_fwet             = 1.0
catch_snow(:,:)    = 1.0
snow_grnd(:,:)     = 1.0

cdnc_ukca_dummy(:,:,:) = 0.0

! Initialise dry_rho; not currently required
DO k = 1, model_levels
  DO j=1, rows
    DO i=1, row_length
      unscl_dry_rho(i,j,k) = 0.0
    END DO
  END DO
END DO

! Array for holding dumps/primary variables/convection sub-step data
ALLOCATE ( resdump(row_length, rows, nprimvars) )

! Allocate Smagorinsky variables
ALLOCATE ( visc_h(1,1,1) )
ALLOCATE ( visc_m(1,1,1) )
ALLOCATE ( rneutml_sq(1,1,1) )
ALLOCATE ( shear(1,1,1) )
ALLOCATE ( max_diff(1,1) )
ALLOCATE ( delta_smag(1,1) )

! Hardwired model switches for SCM

l_mr_pc2     = .FALSE. ! Use mixing ratios in PC2 section
l_murk_bdry  = .FALSE. ! Atmos_physics1 call
lflux_reset  = .FALSE. ! Atmos_physics1 call
l_inland     = .FALSE. ! Re-route inland basin flow to soil moisture

l_phenol  = .FALSE.
rad_mask  = .TRUE.


! Set dimensions of aerosol climatologies - aerosols are not used
! in the SCM so this are hardwired to be unused.
!
! If aerosols are enabled in the future then these should be correctly
! dimensioned using set_arcl_dimensions and set_arcl_clim as in
! atm_step

n_arcl_species = 0
n_arcl_compnts = 1  ! Set to one to prevent out of bounds error in radiation

l_use_arcl(ip_arcl_biom) = .FALSE.
l_use_arcl(ip_arcl_blck) = .FALSE.
l_use_arcl(ip_arcl_sslt) = .FALSE.
l_use_arcl(ip_arcl_sulp) = .FALSE.
l_use_arcl(ip_arcl_dust) = .FALSE.
l_use_arcl(ip_arcl_ocff) = .FALSE.
l_use_arcl(ip_arcl_dlta) = .FALSE.

ALLOCATE ( arcl(1,1,1,1) )

! Mid Latitude values, (ref J-C-Thelen)
o3_trop_level  = 22
o3_trop_height = 11000
t_trop_level   = 22
t_trop_height  = 11000

! Initialise dummy arguments
rdum0    = 1.0
rdum1    = 1.0
rdum2    = 1.0
rdum3    = 1.0

rdum4    = 1.0
rdum5    = 1.0
rdum6    = 1.0

idum0    = 1
dum_land = 1.0

WRITE(nfor_str,'(I5)')    nfor
WRITE(mx_nfor_str,'(I5)') mx_nobs

CALL alloc_common_scm
CALL alloc_forcing
CALL set_forcing_defaults

ti = 0.0
theta_eg = 0.0

!=============================================================================
! Set up vertical level information
!=============================================================================

OPEN(10, FILE=vert_lev, IOSTAT=istatus, ACTION='READ', STATUS='old', &
     IOMSG=iomessage)

IF (istatus /= 0) THEN

  icode=500
  WRITE(cmessage,'(A)')                                              newline//&
    "Error opening " // TRIM(ADJUSTL(vert_lev)) //                   newline//&
    "IoMsg: "//TRIM(iomessage)
  CALL ereport (routinename, icode, cmessage)

END IF

READ(10,vertlevs)

CLOSE(10)

! Initialise error code
error_code = 0

!-----------------------------------------------------------------------------
! Calculate values of grid radius/heights
!-----------------------------------------------------------------------------
  ! Note full model does this in subroutine setcona
ALLOCATE( r_rho_levels      (row_length,rows,  model_levels)   )
ALLOCATE( r_theta_levels    (row_length,rows,0:model_levels)   )
ALLOCATE( r_layer_centres   (row_length,rows,  model_levels+1) )
ALLOCATE( r_layer_boundaries(row_length,rows,0:model_levels+1) )
ALLOCATE( d_layer           (row_length,rows,  model_levels+1) )

IF (.NOT. ALLOCATED(eta_theta_levels)) THEN
  ALLOCATE (eta_theta_levels(0:model_levels))
END IF
IF (.NOT. ALLOCATED(eta_rho_levels)) THEN
  ALLOCATE (eta_rho_levels(model_levels))
END IF

eta_theta_levels(0:model_levels) = eta_theta (1:model_levels+1)
eta_rho_levels(1:model_levels) = eta_rho (1:model_levels)
z_top_theta = z_top_of_model

! DEPENDS ON: calc_levels
CALL calc_levels                                                       &
   ! In
   ( orog, height_gen_method                                           &
   , bl_levels, model_levels                                           &
   , rows, row_length )


z_th(1:model_levels) = r_theta_levels(1,1,1:model_levels) - &
                                                    r_theta_levels(1,1,0)
z_rh(1:model_levels) = r_rho_levels(1,1,1:model_levels)   -  &
                                                   r_theta_levels(1,1,0)
z_rh(model_levels+1) = 2*z_th(model_levels) - z_rh(model_levels)


! Set level heights for physics routines
! Bottom boundary is the surface
DO j = 1, rows
  DO i = 1, row_length
    r_layer_boundaries(i,j,0) = r_theta_levels(i,j,0)
  END DO
END DO
! The first rho level is not used, otherwise boundaries lie on rho levels
! and centres lie on theta levels
DO k = 1, model_levels - 1
  DO j = 1, rows
    DO i = 1, row_length
      r_layer_centres(i,j,k) = r_theta_levels(i,j,k)
      r_layer_boundaries(i,j,k) = r_rho_levels(i,j,k+1)
      d_layer(i,j,k) = r_layer_boundaries(i,j,k) - r_layer_boundaries(i,j,k-1)
    END DO
  END DO
END DO
! The top boundary is extrapolated using twice the distance to the layer centre
k = model_levels
DO j = 1, rows
  DO i = 1, row_length
    r_layer_centres(i,j,k) = r_theta_levels(i,j,k)
    r_layer_boundaries(i,j,k) = 2.0*r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
    d_layer(i,j,k) = r_layer_boundaries(i,j,k) - r_layer_boundaries(i,j,k-1)
  END DO
END DO
! A layer is specified above the top of the model which may be used by the
! radiation scheme to contain the remainder of the atmosphere. This is taken
! to be the same height as the layer below.
k = model_levels + 1
DO j = 1, rows
  DO i = 1, row_length
    d_layer(i,j,k) = d_layer(i,j,k-1)
    r_layer_centres(i,j,k) = r_layer_centres(i,j,k-1) + d_layer(i,j,k)
    r_layer_boundaries(i,j,k) = r_layer_boundaries(i,j,k-1) + d_layer(i,j,k)
  END DO
END DO


! Set radiation switches in rad_input_mod
lrad_ccrad         = l_ccrad

! Set microphysics switches used in mphys_bypass_mod
l_crystals      = l_mcr_qcf2 .OR. (.NOT. l_psd)
mphys_mod_top   = z_top_of_model


!---------------------------------------------------------------------
!     Initialise the array giving the unit nos. for output of the
!     initial data and anything else that needs setting
!---------------------------------------------------------------------

asteps_since_triffid = 0
nout(:)              = 0

! Initialise
rad_hr       = 0.0
micro_tends  = 0.0

ccwpin   = 1.0
ls_rain  = 0.0
ls_snow  = 0.0
ls_graup = 0.0

t_bg_scm = 0.0
q_bg_scm = 0.0
u_bg_scm = 0.0
v_bg_scm = 0.0
w_bg_scm = 0.0

ug_scm(:,:,:) = 0.0
vg_scm(:,:,:) = 0.0

q_star_scm = 0.0
t_inc_scm  = 0.0
u_inc_scm  = 0.0
v_inc_scm  = 0.0
w_inc_scm  = 0.0

flux_h_scm = 0.0
flux_e_scm = 0.0
ustar_in   = 0.0
qcl_inc    = 0.0
qcf_inc    = 0.0

! Initialise past history for convection

DO j=1, rows
  DO i=1, row_length
    deep_flag(i,j)    = 0.0      ! No previous deep convection
    past_precip(i,j)  = 0.0      ! No previous convective precip.
    past_conv_ht(i,j) = 0.0      ! No previous convection.
  END DO
END DO

! Initialise the prognostic variables in the TKE schemes
! The missing values are automatically initialised in the first call of
! the TKE schemes.
e_trb   = rmdi
tsq_trb = rmdi
qsq_trb = rmdi
cov_trb = rmdi

! ilscnt hardwired initialisation to zero until further
! development to correct functionality. Also removed from
! &INOBSFOR namelist to prevent incorrect usage.

ilscnt = 0

! No of levels for Convective Cloud Amount.
IF (l_3d_cca) THEN
  n_cca_lev = model_levels
ELSE
  n_cca_lev = 1
END IF

!-----------------------------------------------------------------------------
! Write out user notification of ancillary files
!-----------------------------------------------------------------------------

sw_spec_file = sw_control(1)%spectral_file
lw_spec_file = lw_control(1)%spectral_file

WRITE(umMessage,'(A)')                                          &
  newline                                                     //&
  ' Using Ancillary files:'//                          newline//&
  ' ----------------------------------------'                 //&
  '-----------------------------------------'//        newline//&
  ' Vertical levels  : '//TRIM(ADJUSTL(vert_lev))//    newline//&
  ' Forcing namelist : '//TRIM(ADJUSTL(scm_nml))//     newline//&
  ' ----------------------------------------'                 //&
  '-----------------------------------------'//        newline//&
  newline
CALL umPrint(umMessage,src='scm_main')

IF (l_netcdf_obs) THEN
  CALL get_netcdf_obs
END IF

! DEPENDS ON: read_scm_nml
CALL read_scm_nml

IF (timestep_nml /= rmdi) THEN
  timestep = timestep_nml
END IF

DO i=1, row_length
  DO j=1, rows
    gridbox_area_m(i,j) = gridbox_area(i,j) * 1e+6
    deltan(i,j)         = SQRT(gridbox_area_m(i,j)/pi)
  END DO
END DO

!---------------------------------------------------------------------------
! Use McClatchey atmospheric profiles where there is misssing data
!---------------------------------------------------------------------------
CALL alloc_mcc(model_levels)
CALL get_mcc(lat, month_init)

! UM levels ABOVE available data/obs_top
!---------------------------------------
DO k=1, model_levels
  IF (z_th(k) > obs_top) THEN
    ozone   (:,:,k) = mcc_o3(k)
    theta   (:,:,k) = mcc_th(k)
    t_bg  (:,:,:,k) = mcc_th_t(k)
    t_inc (:,:,:,k) = 0.0
    qi      (:,:,k) = mcc_q(k)
    q_bg  (:,:,:,k) = mcc_q(k)
    q_star(:,:,:,k) = 0.0
  ELSE
    WHERE (ozone  (:,:,k) == rmdi) ozone  (:,:,k) = mcc_o3(k)
    WHERE (theta  (:,:,k) == rmdi) theta  (:,:,k) = mcc_th(k)
    WHERE (t_bg (:,:,:,k) == rmdi) t_bg (:,:,:,k) = mcc_th_t(k)
    WHERE (qi     (:,:,k) == rmdi) qi     (:,:,k) = mcc_q(k)
    WHERE (q_bg (:,:,:,k) == rmdi) q_bg (:,:,:,k) = mcc_q(k)
  END IF

  ! Check variables on rho levels
  IF (z_rh(k) > obs_top) THEN
    ui(:,:,k) = ui(:,:,k-1)
    vi(:,:,k) = vi(:,:,k-1)
    wi(:,:,k) = 0.0

    u_inc(:,:,:,k) = u_inc(:,:,:,k-1)
    v_inc(:,:,:,k) = v_inc(:,:,:,k-1)
    w_inc(:,:,:,k) = 0.0
  END IF
END DO

DO k=1, model_levels+1
  IF (z_rh(k) > obs_top) THEN
    p_in(:,:,k) = mcc_rh_p(k)
  ELSE
    WHERE (p_in(:,:,k) == rmdi) p_in(:,:,k) = mcc_rh_p(k)
  END IF
END DO

! UM levels BELOW available data/obs_top
!---------------------------------------
! NOTE: Wind profiles and forcing data on UM levels below the obs_bot
!       (default lowest level in namelist.scm) are set to values at
!       obs_bot.  Users wishing to prescribe values should add
!       additional levels to their namelist.scm
!---------------------------------------------------------------------------
DO k=model_levels, 1, -1
  IF (z_th(k) < obs_bot) THEN
    t_inc (:,:,:,k) = t_inc (:,:,:,k+1)
    q_star(:,:,:,k) = q_star(:,:,:,k+1)
    t_bg  (:,:,:,k) = t_bg  (:,:,:,k+1)
    q_bg  (:,:,:,k) = q_bg  (:,:,:,k+1)
    w_inc (:,:,:,k) = w_inc (:,:,:,k+1)
    wi    (:,:,k)   = wi    (:,:,k+1)
  END IF

  IF (z_rh(k) < obs_bot) THEN
    ui(:,:,k) = ui(:,:,k+1)
    vi(:,:,k) = vi(:,:,k+1)

    u_inc(:,:,:,k) = u_inc(:,:,:,k+1)
    v_inc(:,:,:,k) = v_inc(:,:,:,k+1)
  END IF
END DO

! Impose zero w wind at top of model
wi(:,:,model_levels)      = 0.0
w_inc(:,:,:,model_levels) = 0.0


CALL dealloc_mcc

!---------------------------------------------------------------------------
! Now all input has been read in from combination of netcdf/scm
! forcing namelist, check settings
!---------------------------------------------------------------------------

IF ((stats    .AND. obs)     .OR. (stats .AND. noforce)      &
  .OR. (stats    .AND. geoforce) .OR. (obs   .AND. geoforce)     &
  .OR. (geoforce .AND. noforce) .OR. (obs   .AND. noforce)) THEN

  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Warning: More than one forcing logical is set to |'//          newline//&
    '| true. This may produce unexpected results.       |'//          newline//&
    '|                                                  |'
  CALL umPrint(umMessage,src='scm_main')

  WRITE(umMessage,'(4(A12,L1,A39,A))')                                         &
    '|  noforce: ',noforce, '|',                                       newline,&
    '| geoforce: ',geoforce,'|',                                       newline,&
    '|    stats: ',stats,   '|',                                       newline,&
    '|      obs: ',obs,     '|'
  CALL umPrint(umMessage,src='scm_main')

  WRITE(umMessage,'(A)')                                                       &
    '|                                                  |'//          newline//&
    '===================================================='
  CALL umPrint(umMessage,src='scm_main')

END IF ! Test for multiple forcing options


IF ((.NOT. stats)   .AND. (.NOT. obs)                      &
  .AND. (.NOT. noforce) .AND. (.NOT. geoforce)) THEN

  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| No forcing logical set.                          |'//          newline//&
    '| Setting noforce to true as default               |'//          newline//&
    '===================================================='
  CALL umPrint(umMessage,src='scm_main')
  noforce = .TRUE.

END IF ! Test for at least one forcing option.

IF (l_use_dust .OR. l_dust .OR. l_dust_diag) THEN
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Mineral dust scheme not yet implemented in SCM.  |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_use_dust                                      |'//          newline//&
    '|  l_dust                                          |'//          newline//&
    '|  l_dust_diag                                     |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='

  CALL umPrint(umMessage,src='scm_main')

  l_use_dust  = .FALSE.
  l_dust      = .FALSE.
  l_dust_diag = .FALSE.

END IF ! Test for dust scheme options


IF (l_use_soot_direct   .OR. l_use_soot_indirect .OR.                          &
    l_use_soot_autoconv .OR. l_soot) THEN
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Soot chemistry not yet implemented in the SCM.   |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_use_soot                                      |'//          newline//&
    '|  l_soot                                          |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='

  CALL umPrint(umMessage,src='scm_main')

  l_soot              = .FALSE.
  l_use_soot_direct   = .FALSE.
  l_use_soot_indirect = .FALSE.
  l_use_soot_autoconv = .FALSE.

END IF ! Test for soot chemistry options


IF (l_biomass            .OR. l_use_bmass_direct .OR.                          &
    l_use_bmass_indirect .OR. l_use_bmass_autoconv) THEN

  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Biomass scheme not available in the SCM.         |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_biomass                                       |'//          newline//&
    '|  l_use_bmass_direct                              |'//          newline//&
    '|  l_use_bmass_indirect                            |'//          newline//&
    '|  l_use_bmass_autoconv                            |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='

  CALL umPrint(umMessage,src='scm_main')

  l_biomass             = .FALSE.
  l_use_bmass_direct    = .FALSE.
  l_use_bmass_indirect  = .FALSE.
  l_use_bmass_autoconv  = .FALSE.

END IF

IF (l_ocff               .OR. l_use_ocff_direct .OR.                           &
    l_use_ocff_indirect  .OR. l_use_ocff_autoconv) THEN

  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Fossil-fuel OC scheme not available in the SCM.  |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_ocff                                          |'//          newline//&
    '|  l_use_ocff_direct                               |'//          newline//&
    '|  l_use_ocff_indirect                             |'//          newline//&
    '|  l_use_ocff_autoconv                             |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='

  CALL umPrint(umMessage,src='scm_main')

  l_ocff                = .FALSE.
  l_use_ocff_direct     = .FALSE.
  l_use_ocff_indirect   = .FALSE.
  l_use_ocff_autoconv   = .FALSE.

END IF

IF (l_nitrate              .OR. l_use_nitrate_direct .OR.                      &
    l_use_nitrate_indirect .OR. l_use_nitrate_autoconv) THEN

  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Nitrate scheme not available in the SCM.         |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_nitrate                                       |'//          newline//&
    '|  l_use_nitrate_direct                            |'//          newline//&
    '|  l_use_nitrate_indirect                          |'//          newline//&
    '|  l_use_nitrate_autoconv                          |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='

  CALL umPrint(umMessage,src='scm_main')

  l_nitrate               = .FALSE.
  l_use_nitrate_direct    = .FALSE.
  l_use_nitrate_indirect  = .FALSE.
  l_use_nitrate_autoconv  = .FALSE.

END IF

IF (l_rad_deg) THEN
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Spatial degradation of radiation calculation not |'//          newline//&
    '| applicable in the SCM.                           |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_rad_deg                                       |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='
  CALL umPrint(umMessage,src='scm_main')
  l_rad_deg = .FALSE.
END IF

IF (l_ukca_radaer .OR. l_ukca_radaer_sustrat) THEN
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| UKCA RADAER scheme not available in the SCM.     |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_ukca_radaer                                   |'//          newline//&
    '|  l_ukca_radaer_sustrat                           |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='
  CALL umPrint(umMessage,src='scm_main')
  
  l_ukca_radaer         = .FALSE.
  l_ukca_radaer_sustrat = .FALSE.
END IF

IF (l_ukca_aie1 .OR. l_ukca_aie2) THEN
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| UKCA ACTIVATE scheme not available in SCM.       |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_ukca_aie1                                     |'//          newline//&
    '|  l_ukca_aie2                                     |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='
  CALL umPrint(umMessage,src='scm_main')
  
  l_ukca_aie1    = .FALSE.
  l_ukca_aie2    = .FALSE.
END IF

IF (l_glomap_clim_radaer .OR. l_glomap_clim_radaer_sustrat) THEN
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| GLOMAP_CLIM RADAER scheme not available in SCM.  |'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_glomap_clim_radaer                            |'//          newline//&
    '|  l_glomap_clim_radaer_sustrat                    |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='
  CALL umPrint(umMessage,src='scm_main')
  
  l_glomap_clim_radaer         = .FALSE.
  l_glomap_clim_radaer_sustrat = .FALSE.
END IF

IF (l_glomap_clim_aie1 .OR. l_glomap_clim_aie2) THEN
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| GLOMAP_CLIM ACTIVATE scheme not available in SCM.|'//          newline//&
    '| Setting following logicals to false:             |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_glomap_clim_aie1                              |'//          newline//&
    '|  l_glomap_clim_aie2                              |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='
  CALL umPrint(umMessage,src='scm_main')
  
  l_glomap_clim_aie1    = .FALSE.
  l_glomap_clim_aie2    = .FALSE.
END IF

IF (l_leonard_term) THEN

  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Leonard terms not available in the SCM.          |'//          newline//&
    '| Setting following logical to false:              |'//          newline//&
    '|                                                  |'//          newline//&
    '|  l_leonard_term                                  |'//          newline//&
    '|                                                  |'//          newline//&
    '===================================================='

  CALL umPrint(umMessage,src='scm_main')

  l_leonard_term           = .FALSE.

END IF


!-----------------------------------------------------------------------------
! Set l_cosp to false because it is not currently available
! for use in the SCM
!-----------------------------------------------------------------------------
IF (l_cosp) THEN
  l_cosp = .FALSE.
  icode = -1
  cmessage =                                                          newline//&
    'COSP is not currently available in the SCM.'//                   newline//&
    'Setting l_cosp to false.'
  CALL ereport (routinename, icode, cmessage)
END IF

! Check consistency on land masks/land_points
!============================================
land_ice_points = 0
soil_points     = 0
k = 0
l_lice_point(:)=.FALSE.
l_soil_point(:)=.FALSE.

DO j=1, rows
  DO i=1, row_length
    IF (land_sea_mask(i,j)) THEN
      k = k + 1
      land_index(k) = (j-1)*row_length + i

      IF (land_ice_mask(i,j)) THEN
        land_ice_points = land_ice_points + 1
        land_ice_index(land_ice_points) = (j-1)*row_length+i
        l_lice_point(k) = .TRUE.
      END IF

      IF (soil_mask(i,j)) THEN
        soil_points = soil_points + 1
        soil_index(soil_points) = (j-1)*row_length+i
        l_soil_point(k) = .TRUE.
      END IF

    ELSE
      IF (land_ice_mask(i,j) .OR. soil_mask(i,j)) THEN

        WRITE(umMessage,'(A)')                                                 &
          '===================================================='//    newline//&
          '| LAND_ICE_MASK/SOIL_MASK in inconsistent with sea |'//    newline//&
          '| point (land_sea_mask=.FALSE.)                    |'//    newline//&
          '| Setting following logicals to false:             |'//    newline//&
          '|                                                  |'//    newline//&
          '|  land_ice_mask                                   |'//    newline//&
          '|  soil_mask                                       |'//    newline//&
          '|                                                  |'//    newline//&
          '===================================================='
        CALL umPrint(umMessage,src='scm_main')

        land_ice_mask(i,j) = .FALSE.
        soil_mask(i,j)     = .FALSE.

      END IF
    END IF
  END DO
END DO


IF (land_points /= k) THEN
  Icode = 507
  WRITE(umMessage,'(A)')                                                       &
    '===================================================='//          newline//&
    '| Specified total number of land points and        |'//          newline//&
    '| land_sea_mask are inconsisent                    |'//          newline//&
    '===================================================='//          newline//&
    newline                                                                  //&
    'Run ABORTED'
  CALL umPrint('',src='scm_main')

  WRITE(cmessage,'(A)')                                               newline//&
    'Check that land_points and land_sea_mask in'//                   newline//&
    'namelist are consistent.'
  CALL ereport (routinename, icode, cmessage)
END IF

IF (land_points > 0) THEN

  DO j=1, ntype
    DO i=1, land_points
      IF (frac_typ(i,j) == rmdi) THEN
        WRITE(umMessage,'(A,I2,A)')                                            &
          '===================================================='//    newline//&
          '| Land tile fractions (frac_typ) has not been      |'//    newline//&
          '| fully set. ', ntype,                                              &
                           ' fractions must be specified for  |'//    newline//&
          '| each land point in namelist.                     |'//    newline//&
          '===================================================='
        CALL umPrint(umMessage,src='scm_main')

        cmessage =                                                    newline//&
          'Land tile fractions improperly set'
        icode = 1
        CALL ereport (routinename, icode, cmessage)
      END IF
    END DO
  END DO

  DO i=1, land_points
    IF (SUM(frac_typ(i,:)) /= 1.0) THEN
      WRITE(umMessage,'(A)')                                                   &
        '===================================================='//      newline//&
        '| Total tile_fractions (frac_typ) must sum to 1.0  |'//      newline//&
        '| for each land point.                             |'//      newline//&
        '===================================================='
      CALL umPrint(umMessage,src='scm_main')

      cmessage =                                                      newline//&
        'Land tile fractions do not total 1.0'
      icode = 1
      CALL ereport (routinename, icode, cmessage)
    END IF

    IF (soil_type(i) == imdi) THEN
      WRITE(umMessage,'(A)')                                                   &
        '===================================================='//      newline//&
        '| Soil types (soil_type) must be specified by user |'//      newline//&
        '| for each land point.                             |'//      newline//&
        '===================================================='
      CALL umPrint(umMessage,src='scm_main')

      cmessage =                                                      newline//&
        'land point has unspecified soil type'
      icode = 1
      CALL ereport (routinename, icode, cmessage)
    END IF
  END DO

END IF ! test on land_points > 0

IF (obs .AND. .NOT. l_netcdf_obs) THEN
  IF (nfor == imdi) THEN
    icode=503
    WRITE(umMessage,'(A)')                                                     &
      newline                                                                //&
      "=============================================="//              newline//&
      "| Number of observational forcings (nfor) in |"//              newline//&
      "| SCM namelist (&CNTLSCM) has not been set.  |"//              newline//&
      "=============================================="
    CALL umPrint('',src='scm_main')

    CLOSE(10)

    WRITE(cmessage,'(A)') 'Variable NFOR has not been set'
    CALL ereport(routinename, icode, cmessage)

  ELSE IF (nfor > mx_nobs) THEN
    icode=504
    WRITE(umMessage,'(A)')                                                     &
      newline                                                                //&
      "============================================="//               newline//&
      " Specified nfor(" // TRIM(ADJUSTL(nfor_str))                          //&
      ") > mx_nobs("     // TRIM(ADJUSTL(mx_nfor_str)) // ")."//      newline//&
      " This WILL produce incorrect forcings."//                      newline//&
      " Please check your forcing file:"//                            newline//&
      "   " // TRIM(ADJUSTL(scm_nml)) //                              newline//&
      "============================================="//               newline//&
      newline
    CALL umPrint(umMessage,src='scm_main')

    CLOSE(10)

    WRITE(cmessage,'(A)')                                             newline//&
      "Specified nfor(" // TRIM(ADJUSTL(nfor_str))                           //&
      ") > mx_nobs("     // TRIM(ADJUSTL(mx_nfor_str)) // ")."
    CALL ereport(routinename, icode, cmessage)

  END IF
END IF !  obs .and. .not. l_netcdf_obs


IF ( (ndayin == imdi) .OR.                                        &
     (nminin == imdi) .OR.                                        &
     (nsecin == imdi) ) THEN

  ! Run full length of forcing provided using
  ! observed profiles and observation period.
  ndayin =  INT((nfor-1)*obs_pd/rsec_per_day)
  nminin = (INT((nfor-1)*obs_pd - ndayin*rsec_per_day))/60.0
  nsecin = (nfor-1)*obs_pd - ndayin*rsec_per_day - nminin*60.0

END IF

IF (INT(sec_day/timestep) /= REAL(sec_day)/timestep) THEN

  ! Print a warning if we're going to run for more than a day.
  IF (ndayin + REAL(nminin)/(60.0*24.0)                           &
             + REAL(nsecin)/(3600.0*24.0) >  1.0) THEN
    WRITE(umMessage,'(A)')                                                     &
      '================================================'//            newline//&
      '| Warning: Your timestep is such that you have |'//            newline//&
      '| have requested a non-integer number of steps |'//            newline//&
      '| per day.                                     |'//            newline//&
      '================================================'
    CALL umPrint(umMessage,src='scm_main')
  END IF
END IF

full_daysteps = INT(sec_day/timestep)
nstepsin      = INT((nminin*60 + nsecin)/timestep)
total_nsteps  = ndayin * full_daysteps + nstepsin
scm_timestep  = timestep
recip_timestep=1.0/timestep

IF (MOD(obs_pd,timestep) == 0.0) THEN
  ichgf = INT(obs_pd/timestep)
ELSE
  WRITE(umMessage,'(A)')                                                       &
    'Specified run timestep is not'//                                 newline//&
    'a factor of the obs period.'
  CALL umPrint(umMessage,src='scm_main')

  IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

WRITE(umMessage,'(A)')                                            &
  ' ======================================='                    //&
  '========================================'
CALL umPrint(umMessage,src='scm_main')

WRITE(umMessage,'(A11,T22,F7.1)')                                 &
  ' Timestep: ', timestep
CALL umPrint(umMessage,src='scm_main')
WRITE(umMessage,'(A21,T22,F7.1)')                                 &
  ' Observation period: ', obs_pd
CALL umPrint(umMessage,src='scm_main')
WRITE(umMessage,'(A13,T22,I2,A1,I2,A1,I4,A3,I2,A1,I2)')           &
  ' Start date: ', day_init,'/', month_init,'/', year_init, ',  ' &
                 , hour_init, ':', min_init
CALL umPrint(umMessage,src='scm_main')
WRITE(umMessage,'(A13,T22,I2,A6,TR1,I2,A5,I2,A5)')                &
  ' Run length: ', ndayin, ' days,', INT(nminin/60.0), ' hrs,'    &
                 , nminin - INT(nminin/60.0)*60, ' mins'
CALL umPrint(umMessage,src='scm_main')
WRITE(umMessage,'(A)')                                            &
  ' ======================================='                   // &
  '========================================'
CALL umPrint(umMessage,src='scm_main')


!-----------------------------------------------------------------------------
! Code taken from SETCONA
! Calculation of min_trop_level and max_trop_level
! NOTE: min and max_trop_level are used if climatological aerosols
! are chosen.
!-----------------------------------------------------------------------------
!   The tropopause diagnosed for radiative purposes divides theta-levels
!   considered to lie in the stratosphere from those considered to lie
!   in the troposphere: the tropopause is therefore taken as a
!   rho-level. This level is constrained to lie between heights of
!   z_min_trop and z_max_trop. The level is used in the setting of the
!   background aerosol climatology and of ozone profiles, subject to the
!   choice of appropriate options; additionally it is used in the
!   calculation of diagnostics defined at the tropopause.
!
!   Start at the second rho-level because the first is omitted from the
!   grid seen in the physics.
!-----------------------------------------------------------------------------
!   This code is only used if min_trop_level and max_trop_level in
!   &RUNDATA namelist are set both set to 0.  Otherwise values in scm
!   namelist are used.
!-----------------------------------------------------------------------------

IF ((min_trop_level == 0 .AND. max_trop_level == 0) .OR.          &
     min_trop_level <  0  .OR. max_trop_level <  0) THEN

  min_trop_level=2

  DO k=1, model_levels
    r_ref_rho(k) = eta_rho(k)*z_top_of_model
  END DO

  DO ; IF ((r_ref_rho(min_trop_level) >= z_min_trop) .OR.         &
          (min_trop_level == model_levels)) EXIT
    min_trop_level = min_trop_level+1
  END DO

  max_trop_level=min_trop_level
  DO ; IF ( (r_ref_rho(max_trop_level) > z_max_trop) .OR.         &
          (max_trop_level == model_levels) ) EXIT
    max_trop_level = max_trop_level+1
  END DO

  max_trop_level = max_trop_level-1
END IF

!-----------------------------------------------------------------------------

! Derive the initial daynumber in the year and the initial time
! from the UM data structure supplied

! DEPENDS ON: inittime_scm
CALL inittime_scm( dayno_init, time_initi, lcal360 )

time_init = time_initi    ! type real expected elsewhere.

! Initial tape daynumber in the year is the same as the initial
! daynumber

tapeday_init = dayno_init

previous_time(1) = year_init
previous_time(2) = month_init
previous_time(3) = day_init
previous_time(4) = hour_init
previous_time(5) = min_init
previous_time(6) = sec_init
previous_time(7) = dayno_init

! Set time in scm_utils
time_info = previous_time


!-----------------------------------------------------------------------------
! Set longitude to zero so that diagnostics apply to local time
! rather than GMT and convert to radians
!-----------------------------------------------------------------------------

ALLOCATE( true_latitude(row_length,rows))
ALLOCATE( true_longitude(row_length,rows))

DO j=1, rows
  DO i=1, row_length
    true_latitude(i,j) = pi_over_180 * lat(i,j)
    IF (local_time) THEN
      true_longitude(i,j) = 0.0
    ELSE
      true_longitude(i,j) = pi_over_180 * long(i,j)
    END IF
  END DO
END DO


! Allocate and Calculate trig arrays - stored in module
! Full model does this in routine setcona
ALLOCATE( cos_theta_latitude (row_length,rows) )
ALLOCATE( sec_theta_latitude (row_length,rows) )
ALLOCATE( sin_theta_latitude (row_length,rows) )

ALLOCATE( cos_theta_longitude (row_length,rows) )
ALLOCATE( sin_theta_longitude (row_length,rows) )
ALLOCATE( fv_cos_theta_latitude (row_length,rows) )

DO j=1, rows
  DO i=1, row_length
    cos_theta_latitude    (i,j) = COS(true_latitude(i,j))
    sec_theta_latitude    (i,j) = 1.0/cos_theta_latitude(i,j)
    cos_theta_longitude   (i,j) = COS(true_longitude(i,j))
    sin_theta_latitude    (i,j) = SIN(true_latitude(i,j))
    sin_theta_longitude   (i,j) = SIN(true_longitude(i,j))
    FV_cos_theta_latitude (i,j) = cos_theta_latitude(i,j)
  END DO
END DO

! Initialise various arrays
SW_incs     = 0.0
LW_incs     = 0.0
dirpar_incs = 0.0
t1_sd       = 0.0
q1_sd       = 0.0

! If using radiation...
IF (l_radiation) THEN

  ! Set timestep number of first radiation timestep (needed each timestep
  ! in set_l_rad_step to set whether to call radiation)
  it_rad1 = ntrad1

  ! Determine number of model timesteps between the radiation
  ! diagnostics/prognostics timesteps, and calculate the radiation
  ! timestep length.
  CALL set_a_radstep()

END IF

! In SCM A_LW_SEGMENTS should be set to 1 and segment sizes for
! LW and SW should be unset (i.e. -99)
a_lw_segments = 1
a_lw_seg_size = -99
a_sw_seg_size = -99

!-----------------------------------------------------------------------------
! Set l_triffid to false because it is not currently available
! for use in the SCM (UM vn6.2). See Chris D Jones of the
! Terrestrial Carbon Cycle Group, Hadley Centre if you wish to
! use dynamical vegetation.
!-----------------------------------------------------------------------------
IF (l_triffid) THEN
  WRITE(umMessage,'(A)')                                                       &
    '================================================='//             newline//&
    '| Warning: you have requested 19_2A interactive |'//             newline//&
    '| vegetation distribution, which is not         |'//             newline//&
    '| currently available in the SCM.               |'//             newline//&
    '| Setting l_triffid to false                    |'//             newline//&
    '================================================='
  CALL umPrint(umMessage,src='scm_main')
  l_triffid = .FALSE.
END IF

!-----------------------------------------------------------------------------
! Calculate values of delta_lambda and delta_phi
!-----------------------------------------------------------------------------
DO j=1, rows
  DO i=1, row_length
    delta_lambda = SQRT ( gridbox_area_m(i,j) /                            &
                 ( r_theta_levels(i,j,0) * r_theta_levels(i,j,0)           &
                  * fv_cos_theta_latitude(i,j) ) )
    delta_phi = delta_lambda
  END DO
END DO

!-----------------------------------------------------------------------------
! Calculate values of exner
!-----------------------------------------------------------------------------
l_calc_exner = .TRUE.
l_calc_rho   = .TRUE.

! DEPENDS ON: calc_press
CALL calc_press                                                             &
  ! (In)
  ( rows, row_length, p_in, theta, qi                                       &
  , l_calc_exner, l_calc_rho                                                &
  ! (InOut)
  , rho                                                                     &
  ! (Out)
  , exner_theta_levels, exner_rho_levels, p_theta_levels, rp, rp_theta      &
  , p_star )


!-----------------------------------------------------------------------------
! Convert thetai to ti
!-----------------------------------------------------------------------------
DO k=1,model_levels
  DO j=1,rows
    DO i=1,row_length
      ti(i,j,k) = theta(i,j,k)*exner_theta_levels(i,j,k)
    END DO
  END DO
END DO

!-----------------------------------------------------------------------------
! Set up the unit nos. for output.
!-----------------------------------------------------------------------------

nout(1)=6
IF (test)          nout(2) = 22
IF (prindump_step) nout(3) = 30
IF (prindump_day)  nout(4) = 31

IF (prindump_days) THEN
  IF (dump_days(1) > 1) nout(5) = 32
  IF (dump_days(2) > 1) nout(6) = 33
  IF (dump_days(3) > 1) nout(7) = 34
  IF (dump_days(4) > 1) nout(8) = 35
END IF

IF (grafdump_step) nout(9)  = 37
IF (grafdump_day)  nout(10) = 38

IF (grafdump_days) THEN
  IF (dump_days(1) > 1) nout(11) = 39
  IF (dump_days(2) > 1) nout(12) = 40
  IF (dump_days(3) > 1) nout(13) = 41
  IF (dump_days(4) > 1) nout(14) = 42
END IF

IF (obs .AND. prindump_obs) THEN
  DO i=1, 5
    nout(i+14) = 42 + i
  END DO
END IF

!-----------------------------------------------------------------------------
! Write out initial data for run to standard output and to all
! the units to which diagnostics will be written.
!-----------------------------------------------------------------------------

! DEPENDS ON: print_initdata
CALL print_initdata                                                          &
  ( row_length, rows, land_points                                            &
  , ozone_levels, nfor, bl_levels, st_levels, sm_levels                      &
  , ntrop, dayno_init, a_sw_radstep_diag, a_sw_radstep_prog, t_inc           &
  , q_star, u_inc, v_inc, w_inc, ichgf, ilscnt, ti, time_init, nout, 19 )

!-----------------------------------------------------------------------------
! Large scale cloud
!-----------------------------------------------------------------------------
! qcl/qcf values initialised in Run_Init, except for when
! geostropic forcing is used.

cf                    = 0.0
area_cloud_fraction   = 0.0
cfl                   = 0.0
cff                   = 0.0
qcl                   = 0.0
qcf                   = 0.0

! Initialise extra microphysics variables to zero as they are not
! set in Run_Init. Forcing options are currently not available for
! these variables.

qcf2        = 0.0
qrain       = 0.0
qgraup      = 0.0
qcf2_star   = 0.0
qrain_star  = 0.0
qgraup_star = 0.0

! Initialise electric variables. Set to zero for now
flash_pot(:,:,:) = 0.0

! Initialise convection prognostics to zero, if they are used
IF ( l_conv_prog_group_1 ) conv_prog_1(:,:,:)      = 0.0
IF ( l_conv_prog_group_2 ) conv_prog_2(:,:,:)      = 0.0
IF ( l_conv_prog_group_3 ) conv_prog_3(:,:,:)      = 0.0
IF ( l_conv_prog_precip )  conv_prog_precip(:,:,:) = 0.0

! convective cold pool (CCP) prognostics
! not yet recommended for scm use

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ux_ccp(i,j)   = rmdi
    uy_ccp(i,j)   = rmdi
    um_ccp(i,j)   = rmdi
    g_ccp(i,j)    = rmdi
    h_ccp(i,j)    = rmdi
    riso_ccp(i,j) = rmdi
    rdir_ccp(i,j) = rmdi
  END DO
END DO

! ----------------------------------------------------------------------+-------
! Initialise total precipitation prognostic to zero, if used
! More rigorous to define l_totppn_softs but (l_totppn_softs == l_use_ussp) for
! code as introduced, hence l_use_ussp is sufficient unless requirement moves
! beyond simple use for non-orographic gravity scheme source calculation only.
! ----------------------------------------------------------------------+-------
IF ( l_use_ussp ) totalppn(:,:)  = 0.0

daynumber = dayno_init
year = 1

!=============================================================================
! Options to set initial profiles
!=============================================================================
!     (i)   Observational large scale forcing (OBS=TRUE of
!           Namelist LOGIC)
!           Initial data is then from namelist INPROF
!     (ii)  Statistical large scale forcing (STATS=TRUE of
!           Namelist LOGIC)
!           Initial data can either be derived from climate datasets
!           using subroutine INITSTAT or set from namelist
!           INPROF (set ALTDAT=TRUE in namelist LOGIC)
!     (iii) No large-scale forcing initial data is set from namelist
!           INPROF
!     (iv)  Continuation from previous run stored on tape
!     (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!     is overwritten
!=============================================================================

! DEPENDS ON: run_init
CALL run_init                                                               &
  ! (In)
  ( row_length, rows, land_points, land_sea_mask                            &
  , nfor, bl_levels, st_levels, sm_levels                                   &
  , ntiles, nice, nice_use, ntrop, n_cca_lev                                &
  , dayno_init, ichgf                                                       &
  , sec_day, rhcrit(1:model_levels), lcal360, l_snow_albedo   &

  ! (InOut)
  , rho, ti, smcl                                                   &

  ! (Out)
  , iseed, dayno_wint, iccb, icct, cca                                      &
  , cf(:,:,1:), cfl(:,:,1:), cff(:,:,1:), t, q(:,:,1:)                      &
  , qcl(:,:,1:), qcf(:,:,1:), theta_star(:,:,1:)                            &
  , p_star, tstar, u, v, w, w_adv, z0msea, zh                               &
  , t_deep_soil, canopy_gb, tsi, smc, sthf, sthu, snodep, catch, infil_tile &
  , z0_tile, z0h_tile, catch_snow                                           &
  , exner_rho_levels, exner_theta_levels, deltap, p                         &
  , p_theta_levels, rp, rp_theta, ch_t_inc, ch_q_star, ch_u_inc, ch_v_inc   &
  , ch_w_inc, ch_flux_e, ch_flux_h, ch_tstar_forcing, ch_ustar_forcing      &
  , ch_t_bg, ch_q_bg, ch_u_bg, ch_v_bg, ch_w_bg, ch_ug, ch_vg               &
  , flux_e_scm, flux_h_scm, ustar_in, t_inc_scm                             &
  , u_inc_scm, v_inc_scm, w_inc_scm, q_star_scm(:,:,1:), ug_scm, vg_scm     &
  , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                        &
  , tls, qls, uls, vls, wls                                                 &
  , dap1, dap2, dap3, dab1, dab2, dab3                                      &
  , b_exp, hcap, hcon, satcon, sathh, v_sat, v_wilt, v_crit, z0m_soil       &
  , atime, btime, alfada, dbara, dgrada, pa, tbara, tgrada, tsda            &
  , vnbara, vnsda, vpbara, wbara, wsda, alfadb, dbarb, dgradb, pb, tbarb    &
  , tgradb, tsdb, vnbarb, vnsdb, vpbarb, wbarb, wsdb                        &
  , cca_dp, cca_md, cca_sh, bl_w_var )


! Should be set in JULES module
rho_snow_grnd = i_rho_snow_grnd
snowdepth = i_snowdepth
IF (nsmax > 0) THEN
  !  Variables required only with the multilayer scheme.
  !  nsnow is not included, since it does not have a
  !  corresponding variable prefixed with "i_".
  ds = i_ds
  sice = i_sice
  sliq = i_sliq
  tsnowlayer = i_tsnowlayer
  !  rho_snow not required.
  rgrainl = i_rgrainl
END IF

!==========================================================
! cf is initialised in run_init,
! area_cloud_fraction initialised to the same value as cf (bulk fraction)
!==========================================================
DO i=1, row_length
  DO j=1, rows
    DO  k=1, model_levels
      area_cloud_fraction(i,j,k) = cf(i,j,k)
    END DO  ! k
  END DO  ! j
END DO  ! i


! Initialise all surf temperature variables to TSTARI from
! namelist &INPROF
DO i=1, row_length
  DO j=1, rows
    IF (tstar_sea(i,j) == rmdi) THEN
      tstar_sea(i,j)  = tstar(i,j)
    END IF

    IF (tstar_sice(i,j) == rmdi) THEN
      tstar_sice(i,j) = tstar(i,j)
    END IF

    IF (tstar_land(i,j) == rmdi) THEN
      tstar_land(i,j) = tstar(i,j)
    END IF
  END DO
END DO

DO i=1, nlnd
  DO itype=1, ntype
    IF (tstar_tile(i,itype) == rmdi) THEN
      tstar_tile(i,itype) = tstar(i,j)
    END IF
  END DO
END DO

!-----------------------------------------------------------------------------
! For geostrophic forcing
!-----------------------------------------------------------------------------

ALLOCATE( f3_at_u(row_length,rows) )
DO j=1, rows
  DO i=1, row_length
    f_coriolis(i,j) = 2.0 * omega * SIN(lat(i,j) * pi_over_180)
    f3_at_u(i,j) = two_omega * sin_theta_latitude(i,j)
  END DO
END DO

IF (geoinit .AND. geoforce) THEN

  kill_scmdiags(:) = .FALSE.

  DO i=1, row_length
    DO j=1, rows

      modug(i,j) = 0.0
      DO k=1, model_levels
        modug(i,j) = MAX( modug(i,j)                                   &
                        , SQRT(  ug_scm(i,j,k)*ug_scm(i,j,k)           &
                               + vg_scm(i,j,k)*vg_scm(i,j,k)) )
        u(i,j,k) = ug_scm(i,j,k)
        v(i,j,k) = vg_scm(i,j,k)
      END DO

    END DO
  END DO

  !-----------------------------------------------------------------------------
  ! Form restart dump
  !-----------------------------------------------------------------------------

  resdump(:,:,:) = 0.0

  CALL restart_dump                                                         &
    ( row_length, rows, land_points, bl_levels                              &
    , st_levels, sm_levels, n_cca_lev, land_sea_mask, u, v, w               &
    , t, theta, q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:), cf(:,:,1:)             &
    , p, rho, t_deep_soil, smc, canopy_gb                                   &
    , snodep, tstar, zh, z0msea, cca, cca_dp, cca_md, cca_sh                &
    , iccb, icct, smcl, bl_w_var )
  ! Note: if you add more fields to the dump used here (array resdump),
  ! you need to increase the declared array size for resdump, set in the
  ! module scm_utils.
  ! You also need to modify the routine which reads the dump consistently;
  ! dumpinit (called from subroutine runinit and from scm_main).


  DO i=1, row_length
    DO j=1, rows
      rccb(i,j) = REAL(iccb(i,j))
      rcct(i,j) = REAL(icct(i,j))
    END DO
  END DO

  timestep_number = 0
  maxinc = 9999.0

  DO WHILE (maxinc  >   0.1                                                 &
    .AND.  timestep_number  <   full_daysteps)

    timestep_number = timestep_number + 1
    daycount  = 1

    ! DEPENDS ON: timecalc
    CALL timecalc                                                           &
      ( dayno_init, time_init, time_string, lcal360, yearno, day, time_sec  &
      , previous_time, ihour, imin )

    time_info = previous_time

    ! Call the pre-physics routine to set up variables
    ! DEPENDS ON: pre_physics
    CALL pre_physics                                                        &
      ! (In)
      ( row_length, rows, nfor, ichgf, qcl(:,:,1:), qcf(:,:,1:)             &
      , ch_ug, ch_vg, ilscnt, f_coriolis, lcal360, daycount                 &
      , timestep_number, a_sw_radstep_diag, a_sw_radstep_prog               &
      , l_triffid, npft                                                     &
      ! (InOut)
      , u, v, ug_scm, vg_scm, npft_trif                                     &
      ! (Out)
      , co2_mmr )

    ! Reset the radiation timestep flags to false during the initialisation
    ! (they are set in pre_physics)
    l_rad_step_prog = .FALSE.
    l_rad_step_diag = .FALSE.

    !         Now call the physics routines with flags set to
    !         only enable the boundary layer call.

    ! ensure this is always 0 before ap2 call
    zlcl_mixed = 0.0
    l_emcorr_opt = .FALSE.

    
    DO  k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          theta_eg(i,j,k)=theta(i,j,k)
        END DO
      END DO
    END DO

    ! DEPENDS ON: atmos_physics2
    CALL atmos_physics2                                                      &
    ! Parallel variables
       ( row_length, rows, n_proc, n_procx, n_procy, rows, row_length        &
       , numcycles, cycleno                                                  &

    ! Model dimensions.
       , row_length, rows, rows, land_points                                 &
       , bl_levels, st_levels, sm_levels, cloud_levels, land_ice_points      &
       , soil_points, n_cca_lev, ntiles, tr_levels                           &
       , first_constant_r_rho_level, dim_cs1, dim_cs2                        &

    ! Model switches
       , l_dry, l_lbc_old                                                    &
       , lcal360                                                             &

    ! Model Parameters
       , rhcrit(1:model_levels), co2_mmr                                     &
       , tr_vars, tr_ukca                                                    &

    ! In: coordinate information
       , unscl_dry_rho                                                       &
       , delta_lambda, delta_phi, dlambda_p, dphi_p                          &
       , wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                        &
       , lat_rot_np, long_rot_np, f3_at_u                                    &

    ! In: time stepping information.
       , yearno, day, ihour, imin, time_sec                                  &

    ! River routing
       , row_length, row_length, xpa, xua, xva, ypa, yua, yva                &
       , g_p_field, g_r_field, a_steps_since_riv, river_row_length           &
       , river_rows, global_river_row_length, global_river_rows              &
       , river_vel, river_mcoef, i_river_vn                                  &
       , trivdir, trivseq, twatstor, inlandout_atm                           &

    ! Lake evaporation:
       , acc_lake_evap                                                       &

    ! Grid-to-grid river routing
       , r_area, slope, flowobs1, r_inext, r_jnext, r_land                   &
       , substore, surfstore, flowin, bflowin,                               &

    ! diagnostic info
         stashwork3, stashwork5, stashwork8, stashwork9, stashwork19         &
       , stashwork26                                                         &

    ! SCM diagnostics
       , nscmdpkgs, kill_scmdiags, conv_mode, l_emcorr_opt                   &

    ! In: data fields.
       , theta_eg, q, qcl, qcf, qrain, qgraup, qcf2                          &
       , rho, u, v, w, w_adv, p, p_star                                      &
       , exner_rho_levels, exner_theta_levels, land_sea_mask, p_theta_levels &
    ! Mixing ratio prognostics; passed in separately to the q fields so
    ! that any convection scheme which works with mixing ratios
    ! can just use them instead of having to convert the q-fields back to
    ! mixing ratios, which would seem perverse!
       , mix_v, mix_cl, mix_cf, mix_cf2, mix_rain, mix_graup                 &
    ! Note: at this stage, the mixing ratios aren't set, but convection
    ! shouldn't be called in this geostrophic spinup call to physics2,
    ! so these shouldn't be used.

    ! Ancillary fields and fields needed to be kept from timestep to
    ! timestep
       , land_index, land_ice_index, soil_index, canopy_gb, snodep, hcon     &
       , hcap, v_crit, v_wilt, v_sat, z0m_soil, sthf, sthu, sil_orog_land    &
       , ho2r2_orog, sd_orog_land, di, ice_fract, u_0, v_0, u_0, v_0         &
       , cca_dp, cca_md, cca_sh                                              &
       , cca, iccb, icct, cclwp, ccw, lcbase                                 &
       , t_deep_soil, tsi, ti_n, ice_k_n, tstar                              &
       , z0msea, ice_fract_n, ice_thick_n, satcon, sathh, b_exp, smcl        &
       , t1_sd, q1_sd, zh, ddmfx, area_cloud_fraction, cf                    &
       , cfl, cff, ls_rain, ls_rainfrac, ls_snow, ls_graup                   &
       , micro_tends, totalppn, photosynth_act_rad, rad_hr                   &
       , soil_clay, soil_silt, soil_sand, dust_mrel1, dust_mrel2, dust_mrel3 &
       , dust_mrel4, dust_mrel5, dust_mrel6                                  &
       , so2_high_level, so2_em, nh3_em, dms_em, soot_hilem                  &
       , soot_em, ocff_hilem, ocff_em, co2_emits, co2flux                    &
       , deep_flag, past_precip, past_conv_ht                                &

    ! In/Out
       , theta_star                                                          &
       , q_star_scm, qcl_star, qcf_star, qrain_star, qgraup_star, qcf2_star  &
       , cf_star, cfl_star, cff_star                                         &
       , u_inc_scm, v_inc_scm, w_inc_scm(:,:,1:model_levels)                 &
       , sum_eng_fluxes, sum_moist_flux                                      &

    ! Tracer fields(InOut)
       , aerosol, free_tracers, tracer_ukca                                  &
       , dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6    &
       , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new             &
       , soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld               &
       , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss, co2             &

    ! RIVERS (InOut)
       , tot_surf_runoff, tot_sub_runoff                                     &

    ! Out: fields
       , ntml, cumulus, nbdsc, ntdsc                                         &
       , rhcpt,row_length, rows                                              &
       , zlcl_mixed                                                          &
    ! Additional variables for Land Surface Scheme
       , frac_typ, frac_disturb, canht, lai                                  &
       , canopy(1:land_points,1:ntiles)                                      &
       , catch(1:land_points,1:ntiles), catch_snow, snow_grnd, snow_tile     &
       , z0_tile, z0h_tile, tstar_tile, tsurf_elev_surft                     &
       , infil_tile(1:land_points,1:ntiles)                                  &
       , rgrain(1:land_points,1:ntiles)                                      &
       , cs, gs, co2_dim_row, co2_dim_len                                    &
       , asteps_since_triffid, timestep_number                               &
       , g_leaf_acc                                                          &
       , g_leaf_phen_acc, npp_ft_acc, resp_w_ft, resp_s_acc                  &
       , land_pts_trif, npft_trif, dolr, lw_down, sw_tile, fland_ctile       &
       , tstar_land, tstar_sea, tstar_sice_n, tstar_sice                     &
       , albsoil, cos_zenith_angle                                           &

    ! INOUT variables for TKE based turbulence schemes
       , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                        &
    ! OUT variable to store BL w-variance diagnostic
       , bl_w_var                                                            &
    ! IN/OUT convection prognostics
       , conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip             &
       , ux_ccp, uy_ccp, um_ccp, g_ccp, h_ccp, riso_ccp, rdir_ccp            &
    ! Additional variables required for large-scale hydrology:
       , fexp, gamtot, ti_mean, ti_sig, fsat, fwetl, zw                      &
       , sthzw, a_fsat, c_fsat, a_fwet, c_fwet                               &

    ! JULES 2 prognostics (InOut)
       , snowdepth, rho_snow_grnd                                            &
       , nsnow, ds, sice, sliq, tsnowlayer, rho_snow, rgrainl                &

    ! FLake lake scheme prognostics (InOut)
       , lake_depth, lake_fetch, lake_t_mean, lake_t_mxl                     &
       , lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape                     &
       , lake_g_dt                                                           &

    ! Cariolle ozone
       , ozone_tracer                                                        &

    ! Additional screen-level variables
       , tscrndcl_ssi, tscrndcl_tile, tstbtrans                              &

    ! Variables required for COSP (out)
       , cosp_crain_3d, cosp_csnow_3d                                        &

    ! Prescribed surface forcing and roughness lengths
       , flux_e_scm, flux_h_scm, ustar_in, l_spec_z0, z0m_scm, z0h_scm       &

    ! Error information
       , error_code )


    ! Calculate max increment and copy winds for safe-keeping
    maxinc = 0.0
    DO i=1, row_length
      DO j=1, rows
        DO k=1, model_levels
          maxinc = MAX(maxinc,(u_inc_scm(i,j,k)**2 +                        &
                               v_inc_scm(i,j,k)**2))
          ui(i,j,k) = u(i,j,k)
          vi(i,j,k) = v(i,j,k)
        END DO
        maxinc = SQRT(maxinc)/(f_coriolis(i,j) * timestep*modug(i,j))
      END DO
    END DO


    !--------------------------------------------------------------------------
    ! Copy initial data back from DUMP.
    !--------------------------------------------------------------------------

    CALL dumpinit                                                           &
      ( row_length, rows, land_points, bl_levels                            &
      , st_levels, sm_levels, n_cca_lev, land_sea_mask, u, v, w             &
      , t, theta, q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:), cf(:,:,1:)           &
      , p, rho, t_deep_soil, smc, canopy_gb                                 &
      , snodep, tstar, zh, z0msea, cca, cca_dp, cca_md, cca_sh              &
      , rccb, rcct, smcl, bl_w_var )
! Note: if you add more fields to the dump used here (array resdump),
! you need to increase the declared array size for resdump, set in the
! module scm_utils.
! You also need to modify the routine which saves the dump consistently;
! restart_dump (called from scm_main).
! You also need to modify the other call to dumpinit, from run_init.

    DO i=1, row_length
      DO j=1, rows
        iccb(i,j) = INT(rccb(i,j))
        icct(i,j) = INT(rcct(i,j))

        ! Copy saved U,V back

        DO k=1, model_levels
          u(i,j,k) = ui(i,j,k)
          v(i,j,k) = vi(i,j,k)
        END DO
      END DO
    END DO
  END DO                   ! maxinc and timestep_number < daysteps

  resdump(:,:,:) = 0.0

  WRITE(umMessage,'(A,E13.5,A,I5,A)')                                          &
    "Geostrophic wind initialised."//                                 newline//&
    "Max relative wind change at end = ", maxinc, " after ",                   &
    timestep_number," steps"
  CALL umPrint(umMessage,src='scm_main')

  daynumber = dayno_init
  year = 1

END IF                     ! Geostrophic forcing initialising.

!-----------------------------------------------------------------------------
! Initialise the output diagnostic system
!-----------------------------------------------------------------------------
! DEPENDS ON: setup_diags
CALL setup_diags                                               &
  ( row_length, rows, bl_levels, sm_levels                     &! In
  , st_levels, land_points, ntiles, nsmax, n_vis_thresh        &! In
  , cloud_levels                                               &! In
  , total_nsteps, timestep, full_daysteps                      &! In
  , a_sw_radstep_prog, a_sw_radstep_diag, ntrad1               &! In
  , daycount, timestep_number, nscmdpkgs                       &! In
  , l_scmdiags, scmop)                ! Out

! If PC2 is off then l_SCMDiags(SCMDiag_pc2) must be false
IF (i_cld_vn /= i_cld_pc2) THEN
  IF (l_scmdiags(scmdiag_pc2)) THEN
    WRITE(umMessage,'(A)')                                                     &
      ' ==================================='//                        newline//&
      ' | Warning: you have requested PC2 |'//                        newline//&
      ' | diagnostics but PC2 is off,     |'//                        newline//&
      ' | resetting l_SCMDiag_PC2 logical |'//                        newline//&
      ' ==================================='
    CALL umPrint(umMessage,src='scm_main')
  END IF
  l_scmdiags(scmdiag_pc2) = .FALSE.
END IF

! The availability of Surface based diagnostics packages is determined
! by the surface type.
!
! Surface Package     (SCMDiags_surf) - always available
! Land Points Package (SCMDiags_land) - Only if land_sea_mask TRUE
! Sea Points Package  (SCMDiags_sea)  - Only if land_sea_mask FALSE
!
! This only works in the SCM because rows, row_length are size 1 here

DO j=1, rows
  DO i=1, row_length
    IF (land_sea_mask(i,j)) THEN

      ! Land point
      IF (l_scmdiags(scmdiag_sea)) THEN
        WRITE(umMessage,'(A)')                                                 &
          ' ============================================'//           newline//&
          ' | Warning: You have requested sea          |'//           newline//&
          ' | diagnostics but this is a land point,    |'//           newline//&
          ' | resetting l_SCMDiag_sea logical to FALSE |'//           newline//&
          ' ============================================'
        CALL umPrint(umMessage,src='scm_main')
      END IF

      l_scmdiags(scmdiag_sea) = .FALSE.

    ELSE

      ! Sea point
      IF (l_scmdiags(scmdiag_land)) THEN
        WRITE(umMessage,'(A)')                                                 &
          ' ============================================='//          newline//&
          ' | Warning: You have requested land          |'//          newline//&
          ' | diagnostics but this is a sea point,      |'//          newline//&
          ' | resetting l_SCMDiag_land logical to FALSE |'//          newline//&
          ' ============================================='
        CALL umPrint(umMessage,src='scm_main')
      END IF

      l_scmdiags(scmdiag_land) = .FALSE.

      IF (l_scmdiags(scmdiag_mlsnow)) THEN
        WRITE(umMessage,'(A)')                                &
          ' ===============================================' // newline //&
          ' | Warning: You have requested multilayer      |' // newline //&
          ' | snow diagnostics but this is a sea point,   |' // newline //&
          ' | resetting l_SCMDiag_mlsnow logical to FALSE |' // newline //&
          ' ==============================================='
        CALL umPrint(umMessage,src='scm_main')
      END IF

      l_scmdiags(scmdiag_mlsnow) = .FALSE.

    END IF ! land_sea_mask
  END DO
END DO

WRITE(umMessage,'(A)')                                                         &
  'Initial profiles of water species are considered as '                     //&
  'SPECIFIC values' //                                                newline//&
  'Forcing profiles of water species are considered as '                     //&
  'SPECIFIC values'
CALL umPrint(umMessage,src='scm_main')

IF (l_mr_physics) THEN
  WRITE(umMessage,'(A)')                                                     &
    ' Atmos_Physics 1 and 2 will consider water species as MIXING RATIOS values'
ELSE
  WRITE(umMessage,'(A)')                                                     &
    ' Atmos_Physics 1 and 2 will consider water species as SPECIFIC values'
END IF
CALL umPrint(umMessage,src='scm_main')

!-----------------------------------------------------------------------------
! Loop over days
!-----------------------------------------------------------------------------

! Timestepping proper is about to begin, switch the
! diagnostic system on (as long as at least one stream is open)

IF (main_diag_switch /= 0 .AND. ANY(scmop%strm%switch /= 0)) THEN
      ! any() is a Fortran90 intrinsic funtion
  scmop%on = .TRUE.
END IF


scm_timestep_count = 0

DO daycount=1, ndayin+1

  !----------------------------------------------------------------------------
  ! Reset sinusoidal distribution every 10 days if climate stats required
  !----------------------------------------------------------------------------

  IF (stats .AND. (daycount  ==  1                                          &
    .OR. (ancyc .AND. MOD(daycount, change_clim)  ==  0))) THEN
    ! DEPENDS ON: statday
    CALL StatDay                                                            &
      ! (In)
      ( row_length, rows, ntrop                                             &
      , atime, btime, dayno_wint, deltan, daycount                          &
      , tbara, tbarb, tsda, tsdb, dbara, dbarb, vnbara, vnbarb              &
      , vnsda, vnsdb, vpbara, vpbarb, wbara, wbarb, wsda, wsdb              &
      , alfada, alfadb, pa, pb, p, tgrada                                   &
      , tgradb, dgrada, dgradb, cort, cord, corvn, corw                     &
      , tdash, ddash, ctbar, ctsd, at, cdbar, cdsd, ad                      &
      , cvnbar, cvnsd, avn, cwbar, cwsd, aw                                 &
      , tbar, tsd, dbar, dsd                                                &
      , vnbar, vnsd, vpbar, wbar, wsd )

    !--------------------------------------------------------------------------
    !     Calculate values of exner
    !--------------------------------------------------------------------------

    l_calc_exner = .TRUE.
    l_calc_rho   = .FALSE.

    ! DEPENDS ON: calc_press
    CALL calc_press                                                         &
      ! (In)
      ( rows, row_length, p                                                 &
      , theta, q(:,:,1:), l_calc_exner, l_calc_rho                          &
      ! (InOut)
      , rho                                                                 &
      ! (Out)
      , exner_theta_levels, exner_rho_levels, p_theta_levels                &
      , rp, rp_theta, p_star)

    !--------------------------------------------------------------------------
    ! Initialise PX and PY arrays for calculation of vertical fluxes later
    !--------------------------------------------------------------------------

    DO i=1, row_length
      DO j=1, rows
        DO k=1, ntrop
          px(i,j,k) = 1.0 / LOG(p(i,j,k+1)/ p(i,j,k))
        END DO

        DO k=1, ntrop-1
          py(i,j,k) = 1.0 / LOG(p(i,j,k+2)/ p(i,j,k))
        END DO
      END DO
    END DO
  END IF                   ! stats

  !============================================================================
  ! Options for forcing
  !----------------------------------------------------------------------------
  !
  !       Observational forcing : use observed values of T,Q,U,V
  !       to calculate their increments due to large scale effects.
  ! OR
  !       Statistical forcing : take random sample from Normal
  !       (Gaussian) distribution with mean and sd climlogical
  !       average to calculate increments to T and Q due to large scale
  !       effects.
  !=============================================================================
  !
  !       Loop over timesteps
  !
  !
  !       If it is the last day in the run and a full day is not
  !       required, loop over the number of timesteps required
  !       otherwise Do the full number of timesteps in a day.

  IF (daycount  ==  ndayin+1 .AND. nstepsin  /=  full_daysteps) THEN
    daysteps = nstepsin
  ELSE
    daysteps = full_daysteps
  END IF

  DO timestep_number=1, daysteps

    ! Update local_timestep_count
    scm_timestep_count = scm_timestep_count + 1

    !--------------------------------------------------------------------------
    ! VARIABLE VALUES: q     is vapour         at timelevel n
    !                  qcl   is liquid         at timelevel n
    !                  t     is temperature    at timelevel n
    !                  theta is potential temp at timelevel n
    !--------------------------------------------------------------------------

    IF (main_diag_switch /= 0) THEN
      !
      !------------------------------------------------------------------------
      ! SCM PC2 Diagnostics Package
      !------------------------------------------------------------------------
      IF (l_scmdiags(scmdiag_pc2)) THEN

        CALL scmoutput(q(:,:,1:), 'q_timen'                                 &
          , 'Vapour at start of timestep', 'kg/kg'                          &
          , t_inst, d_wet, default_streams, '', routinename)

        CALL scmoutput(qcl(:,:,1:), 'qcl_timen'                             &
          , 'Liquid at start of timestep', 'kg/kg'                          &
          , t_inst, d_wet, default_streams, '', routinename)

        CALL scmoutput(t, 't_timen'                                         &
          , 'Temperature at start of timestep', 'K'                         &
          , t_inst, d_all, default_streams, '', routinename)

        CALL scmoutput(theta, 'th_timen'                                    &
          , 'Potential temperature at start of timestep', 'K'               &
          , t_inst, d_all, default_streams, '', routinename)

      END IF ! l_SCMDiags(SCMDiag_PC2)

    END IF ! main_diag_switch /= 0

    ! Reset the wind increments to zero -
    ! The other increments are reset in physics1

    u_inc_scm(:,:,:) = 0.0
    v_inc_scm(:,:,:) = 0.0
    w_inc_scm(:,:,:) = 0.0

    !-------------------------------------------------------------------------
    ! Calculate the year (in run) and actual time and day number
    ! for labelling  of the diagnostics only.
    !-------------------------------------------------------------------------

    ! DEPENDS ON: timecalc
    CALL timecalc                                                           &
      ( dayno_init, time_init, time_string, lcal360, yearno, day, time_sec  &
      , previous_time, ihour, imin)

    time_info = previous_time

    ! If there is no annual cycle, the year and day numbers
    ! will be the init ones.
    IF (.NOT. ancyc) THEN
      year = 1
      day = dayno_init
    END IF

    !-------------------------------------------------------------
    ! Convert temperature to potential temperature
    !         t_inc_scm   to theta_star
    !-------------------------------------------------------------

    DO  k=1, model_levels   
      DO j=1, rows
        DO i=1, row_length

          theta(i,j,k)      = t(i,j,k)         / exner_theta_levels(i,j,k)
          theta_star(i,j,k) = t_inc_scm(i,j,k) / exner_theta_levels(i,j,k)

        END DO
      END DO
    END DO


    IF (i_cld_vn == i_cld_pc2) THEN

      !-----------------------------------------------------------------------
      ! Calculate initial relative humidity w.r.t. Liquid temperature TL
      ! Is used by PC2 routines at end of timestep. Jeremy Price Feb 2005.
      !-----------------------------------------------------------------------
      DO k=1, model_levels

        DO j=1, rows
          DO i=1, row_length
            tl(i,j) = t(i,j,k) - (lc/cp)*qcl(i,j,k)
          END DO ! i
        END DO ! j

        !Change this to SCM once #3386 is on and we merge to HoT
        IF ( l_new_qsat_cntl ) THEN
          CALL qsat_wat_new(qsl_tl,tl,p_theta_levels(:,:,k),row_length,rows)
        ELSE
          ! DEPENDS ON: qsat_wat
          CALL qsat_wat (qsl_tl, tl, p_theta_levels(1,1,k), row_length*rows )
        END IF

        DO j=1, rows
          DO i=1, row_length
            rhts(i,j,k) = (q(i,j,k) + qcl(i,j,k))/qsl_tl(i,j)
            tlts(i,j,k) = tl(i,j)
            ptts(i,j,k) = p_theta_levels(i,j,k)
            qtts(i,j,k) = q(i,j,k) + qcl(i,j,k)
          END DO ! i
        END DO ! j

      END DO ! k

    END IF  ! i_cld_pc2


    !-------------------------------------------------------------------------
    ! Save initial fields for later to calculate total increments
    !-------------------------------------------------------------------------
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          t_start(i,j,k) = t(i,j,k)
          u_start(i,j,k) = u(i,j,k)
          v_start(i,j,k) = v(i,j,k)
        END DO
      END DO
    END DO

    DO k=0, model_levels
      DO j=1, rows
        DO i=1, row_length
          w_start(i,j,k) = w(i,j,k)
        END DO
      END DO
    END DO

    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          q_start(i,j,k)   = q(i,j,k)
          qcl_start(i,j,k) = qcl(i,j,k)
          qcf_start(i,j,k) = qcf(i,j,k)
          cf_start(i,j,k)  = cf(i,j,k)
          cfl_start(i,j,k) = cfl(i,j,k)
          cff_start(i,j,k) = cff(i,j,k)
        END DO
      END DO
    END DO

    !-------------------------------------------------------------------------
    ! Call physics one before advection step
    !-------------------------------------------------------------------------

    ! DEPENDS ON: pre_physics
    CALL pre_physics                                                        &
      ! (In)
      ( row_length, rows, nfor, ichgf, qcl(:,:,1:), qcf(:,:,1:)             &
      , ch_ug, ch_vg, ilscnt, f_coriolis, lcal360, daycount                 &
      , timestep_number, a_sw_radstep_diag, a_sw_radstep_prog               &
      , l_triffid, npft                                                     &
      ! (InOut)
      , u, v, ug_scm, vg_scm, npft_trif                                     &
      ! (Out)
      , co2_mmr )


    !------------------------------------------------------------------|
    ! UM5.x timestepping stores all wind increments in u_inc_scm and   |
    ! v_inc_scm. These are then added to u,v in atmos_physics2.        |
    !------------------------------------------------------------------|

    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length

          !---------------------------------------------------|
          ! Calculate increments from geostrophic forcing     |
          ! These are then added to u_inc_scm and v_inc_scm   |
          ! before atmos_physics2 as atmos_physics1 sets      |
          ! u_inc_scm and v_inc_scm to zero.                  |
          !---------------------------------------------------|

          uinc_geo(i,j,k) = u(i,j,k)-u_start(i,j,k)
          vinc_geo(i,j,k) = v(i,j,k)-v_start(i,j,k)

          !--------------------------------|
          ! Reset u,v back to time-level n |
          !--------------------------------|

          u(i,j,k) = u_start(i,j,k)
          v(i,j,k) = v_start(i,j,k)

        END DO  ! i
      END DO  ! j
    END DO  ! k

    !-----------------------------------
    ! Convert to mixing ratios if needed
    !-----------------------------------
    IF (l_mr_physics .OR. l_mr_conv) THEN

      ! Convert to mixing ratios
      ! DEPENDS ON: q_to_mix
      CALL q_to_mix                                                         &
        ( row_length, rows, halo_i, halo_j                                  &
        , q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:)                               &
        , qcf2(:,:,1:), qrain(:,:,1:), qgraup(:,:,1:)                       &
        , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                             &
        , mix_v(:,:,1:), mix_cl(:,:,1:), mix_cf(:,:,1:)                     &
        , mix_cf2(:,:,1:), mix_rain(:,:,1:), mix_graup(:,:,1:) )

    END IF

    IF (l_mr_physics) THEN

      ! Now place mixing ratio values back into d1
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length

            q(i,j,k)     = mix_v(i,j,k)     ! Vapour
            qcl(i,j,k)   = mix_cl(i,j,k)    ! Liquid
            qcf(i,j,k)   = mix_cf(i,j,k)    ! Ice
            qcf2(i,j,k)  = mix_cf2(i,j,k)   ! Ice 2
            qrain(i,j,k) = mix_rain(i,j,k)  ! Rain
            qgraup(i,j,k)= mix_graup(i,j,k) ! Graupel

          END DO
        END DO
      END DO

    END IF  ! l_mr_physics

    DO  k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          theta_eg(i,j,k)=theta(i,j,k)
        END DO
      END DO
    END DO

    !==========================================================================

    ! DEPENDS ON: atmos_physics1
    CALL atmos_physics1                                                       &
    ! Parallel variables
       ( row_length, rows, n_proc, n_procx, n_procy, rows, row_length         &

    ! model dimensions.
       , row_length, rows, rows, land_points                                  &
       , bl_levels, st_levels, sm_levels                                      &
       , ozone_levels, cloud_levels, land_ice_points, soil_points             &
       , n_cca_lev, ntiles, salt_dim1, salt_dim2, salt_dim3                   &
       , tr_levels, tr_ukca, cdnc_dim1, cdnc_dim2, cdnc_dim3                  &
       , co2_dim_len, co2_dim_row, co2_dim_lev                                &
       , n_arcl_species, n_arcl_compnts, i_arcl_compnts                       &

    ! Model switches
       , l_lbc_old                                                            &

    ! Setting lflux_reset to false but, if using energy correction, will
    ! need to be worked out every timestep
       , l_ukca_chem, l_ukca_set_trace_gases                                  &
       , l_ukca_strat, l_ukca_strattrop                                       &
       , l_ukca_prescribech4, l_use_arcl                                      &

    ! Model Parameters
       , rhcrit(1:model_levels)                                               &
       , min_trop_level, max_trop_level                                       &
       , ngrgas, grgas_addr                                                   &

    ! Parameter for stochastic physics random parameters2
       , m_ci                                                                 &

    ! In: coordinate information
       , delta_lambda, delta_phi, lat_rot_NP, long_rot_NP                     &

    ! In: time stepping information.
       , yearno, day, ihour, imin, isec, previous_time,                       &

    ! Diagnostic info
         stashwork1, stashwork2, stashwork4, stashwork6, stashwork14          &
       , stashwork21                                                          &
    ! SCM diagnostics
       , nscmdpkgs, l_scmdiags                                                &

    ! In data fields.
       , theta_eg, q, qcl, qcf, qcf2, qrain, qgraup, rho,  u, v, w, p         &
       , p_star, exner_rho_levels, exner_theta_levels, land_sea_mask          &
       , p_theta_levels, fland_ctile, frac_control, cdnc_ukca_dummy           &

    ! ancillary fields and fields needed to be kept from timestep to
    ! timestep
       , land_index, rgrain(1:land_points,1:ntiles), soot, canht              &
       , ntml, cumulus                                                        &
       , ice_fract_n, ice_fract_n, ice_thick_n                                &
       , cca_dp, cca_md, cca_sh                                               &
       , cca, iccb, icct, cclwp, ccw, lcbase, totalppn                        &
       , tstar, tstar_land, tstar_sea, tstar_sice_n                           &
       , sice_alb, land_alb, snodep, snodep_sice_n                            &
       , ozone, sw_incs, lw_incs                                              &
       , dirpar_incs, o3_trop_level, o3_trop_height, t_trop_level             &
       , t_trop_height, zh, sd_orog_land, orog_grad_xx_land                   &
       , orog_grad_xy_land, orog_grad_yy_land, area_cloud_fraction            &
       , cf, cfl, cff, aerosol_em, arcl, albsoil, albobs_sw, albobs_vis       &
       , albobs_nir, lai, snow_tile, frac_typ, tstar_tile, z0_tile            &
       , dolr_rts, lw_down, sw_tile_rts, es_space_interp, rad_mask            &
       , cos_zenith_angle                                                     &
       , easyaerosol_sw, easyaerosol_lw, easyaerosol_cdnc                     &

    ! Variables for COSP
       , cosp_crain_3d,cosp_csnow_3d                                          &

    ! JULES 2 prognostics (In)
       , snowdepth,lake_h_ice, z0msea, chloro_sea                             &
    ! In variable storing BL w-variance diagnostic
       , bl_w_var                                                             &
    ! InOut
       , theta_star                                                           &
       , q_star_scm, qcl_star, qcf_star, qcf2_star, qrain_star, qgraup_star   &
       , cf_star, cfl_star, cff_star, u_inc_scm, v_inc_scm                    &
       , energy_correction, sum_eng_fluxes, sum_moist_flux, aerosol           &
       , flash_pot                                                            &
       , dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6     &
       , so2, so4_aitken, so4_accu, so4_diss, nh3                             &
       , soot_new, soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld      &
       , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss                   &
       , co2, free_tracers, tracer_ukca, biogenic, asteps_since_triffid       &
       , ukca_radaer                                                          &

    ! Out
       , ls_rain, ls_rainfrac, ls_snow, ls_graup, micro_tends, unscl_dry_rho  &
       , photosynth_act_rad, rad_hr, dolr, sw_tile                            &

    ! Error information
       , error_code )


          !-------------------------------
          ! Convert to specific humidities
          !-------------------------------

    IF (l_mr_physics) THEN

      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length

            ! Copy q and qstar variables to mix variables
            mix_v_star(i,j,k)     = q_star_scm(i,j,k)  ! Vapour
            mix_cl_star(i,j,k)    = qcl_star(i,j,k)    ! Liquid
            mix_cf_star(i,j,k)    = qcf_star(i,j,k)    ! Ice
            mix_cf2_star(i,j,k)   = qcf2_star(i,j,k)   ! Ice2
            mix_rain_star(i,j,k)  = qrain_star(i,j,k)  ! Rain
            mix_graup_star(i,j,k) = qgraup_star(i,j,k) ! Graupel

          END DO
        END DO
      END DO


      ! Convert mixing ratios back to specific humidities
      ! DEPENDS ON: mix_to_q
      CALL mix_to_q                                                         &
        ( row_length, rows, halo_i, halo_j                                  &
        , mix_v(:,:,1:), mix_cl(:,:,1:), mix_cf(:,:,1:)                     &
        , mix_cf2(:,:,1:), mix_rain(:,:,1:), mix_graup(:,:,1:)              &
        , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                             &
        , q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:)                               &
        , qcf2(:,:,1:), qrain(:,:,1:), qgraup(:,:,1:) )

      ! Convert mixing ratio increments (mix_star) back
      ! to specific humidity increments (q_star_scm)
      ! DEPENDS ON: calc_q_star
      CALL calc_q_star                                                      &
        ( row_length, rows                                                  &
        , halo_i, halo_j, offx, offy                                        &
        , mix_v(:,:,1:), mix_cl(:,:,1:), mix_cf(:,:,1:)                     &
        , mix_cf2(:,:,1:), mix_rain(:,:,1:), mix_graup(:,:,1:)              &
        , mix_v_star, mix_cl_star, mix_cf_star                              &
        , mix_cf2_star, mix_rain_star, mix_graup_star                       &
        , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                             &
        , q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:)                               &
        , qcf2(:,:,1:), qrain(:,:,1:), qgraup (:,:,1:)                      &
        , q_star_scm(:,:,1:), qcl_star(:,:,1:), qcf_star (:,:,1:)           &
        , qcf2_star(:,:,1:), qrain_star(:,:,1:), qgraup_star(:,:,1:) )

    END IF  ! l_mr_physics


    !------------------------------------------------------------------|
    ! Add on increments from geostrophic wind forcing.                 |
    !------------------------------------------------------------------|

    IF (geoforce) THEN
      u_inc_scm(:,:,:) = u_inc_scm(:,:,:) + uinc_geo(:,:,:)
      v_inc_scm(:,:,:) = v_inc_scm(:,:,:) + vinc_geo(:,:,:)
    END IF

    !-------------------------------------------------------------------------
    ! Convert theta_star to t_inc_scm for call to forcing
    !-------------------------------------------------------------------------

    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length

          t_inc_scm(i,j,k) = theta_star(i,j,k) * exner_theta_levels(i,j,k)

          ! set qcl_inc and qcf_inc to increments
          qcl_inc(i,j,k) = qcl_star(i,j,k)
          qcf_inc(i,j,k) = qcf_star(i,j,k)
        END DO
      END DO
    END DO

    !-------------------------------------------------------------------------
    ! VARIABLE VALUES:
    !   qcl_inc and qcl_star = liq. water incs      (atmos_physics1)
    !   q_star_scm           = vapour incs          (atmos_physics1)
    !   t_inc_scm            = temp. incs           (atmos_physics1)
    !   theta_star           = pot. temp.incs       (atmos_physics1)
    !-------------------------------------------------------------------------

    !-------------------------------------------------------------------------
    ! PC2 section.
    ! Store increments from atmos_physics1. We can then work out the forcing
    ! by subtraction of these _earliest values from the values returned after
    ! forcing.
    !-------------------------------------------------------------------------

    IF (i_cld_vn == i_cld_pc2) THEN
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            q_earliest   (i,j,k) = q_star_scm (i,j,k)
            qcl_earliest (i,j,k) = qcl_inc    (i,j,k)
            t_earliest   (i,j,k) = t_inc_scm  (i,j,k)
          END DO
        END DO
      END DO
    END IF  ! i_cld_pc2

    IF (main_diag_switch /= 0) THEN

      !-----------------------------------------------------------------------
      ! SCM PC2 Diagnostics Package
      !-----------------------------------------------------------------------
      IF (l_scmdiags(scmdiag_pc2)) THEN

        CALL scmoutput(q_earliest, 'dq_earliest'                            &
          , 'q_star vapour incs from atmos_phys1', 'kg/kg'                  &
          , t_inst, d_wet, default_streams, '', routinename)

        CALL scmoutput(qcl_earliest, 'dqcl_earliest'                        &
          , 'qcl_inc liq water incs atmos_phys1', 'kg/kg'                   &
          , t_inst, d_wet, default_streams, '', routinename)

        CALL scmoutput(t_earliest, 'dt_earliest'                            &
          , 't_inc temp incs from atmos_phys1', 'K'                         &
          , t_inst, d_all, default_streams, '', routinename)

      END IF ! l_SCMDiags(SCMDiag_PC2)

    END IF ! main_diag_switch /= 0

    !-------------------------------------------------------------------------
    ! Include any forcing required
    !-------------------------------------------------------------------------

    IF (stats .OR. obs) THEN
      ! DEPENDS ON: forcing
      CALL forcing                                                          &
        ! (In)
        ( row_length, rows, nfor, bl_levels                                 &
        , st_levels, sm_levels, ntrop, sec_day, timestep_number, daycount   &
        , dayno_wint, daysteps, nscmdpkgs, ichgf, t                         &
        , q(:,:,1:), qcl(:,:,1:)                                            &
        , qcf(:,:,1:), u, v, w, l_scmdiags, p                               &
        , exner_theta_levels(:,:,1:), rp, r_theta_levels                    &
        , ch_tstar_forcing, ch_ustar_forcing, ch_flux_h, ch_flux_e          &
        , ch_t_inc, ch_q_star, ch_u_inc, ch_v_inc, ch_w_inc                 &
        , ch_t_bg,  ch_q_bg,   ch_u_bg,  ch_v_bg,  ch_w_bg                  &
        , ad, at, avn, aw, cdbar, ctbar, cvnbar, cwbar, dbar                &
        , tbar, vnbar, wbar, vpbar, cdsd, ctsd, cvnsd, cwsd, dsd, tsd, vnsd &
        , wsd, tdash, ddash, deltan, px, py                                 &

        ! (InOut)
        , ilscnt, flux_h_scm, flux_e_scm, ustar_in, tstar, ti, qi, t_inc_scm&
        , q_star_scm(:,:,1:), qcl_inc, qcf_inc                              &
        , u_inc_scm, v_inc_scm, w_inc_scm                                   &
        , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                  &
        , tls, qls, uls, vls, wls, tr, qr, vnr, vpr, wr                     &

        ! (Out)
        , factor_rhokh, rhokh, dab1, dap1 )
    END IF

    !-------------------------------------------------------------------------
    ! VARIABLE VALUES:
    !   qcl_inc    = liquid water increments        (atmos_physics1)
    !              + forcing increments.
    !   t_inc_scm  = temperature increments         (atmos_physics1)
    !              + forcing increments
    !   q_star_scm = vapour increments              (atmos_physics1)
    !              + forcing increments
    !   wls        = The current vertical velocity forcing tendency (m/s)/day
    !-------------------------------------------------------------------------
    ! Update model vertical velocities to be consistent with
    ! subsidence forcing
    !-------------------------------------------------------------------------

          ! Want latest w whether vertical advection on or off.
          ! w = Current w from start of timestep
    w(:,:,:) = w(:,:,:) + w_inc_scm(:,:,:)
    w_adv(:,:,:) = w(:,:,:)


    !-------------------------------------------------------------------------
    ! Update SST to be consistent with tstar_forcing (Coastal tiling)
    !-------------------------------------------------------------------------
    DO j=1, rows
      DO i=1, row_length
        tstar_sea(i,j) = tstar(i,j)
      END DO ! i
    END DO ! j

    !-------------------------------------------------------------------------
    ! This is a PC2 section.
    ! Calculate the forcing of vapour, liquid, temperature and pressure.
    !-------------------------------------------------------------------------

    IF (i_cld_vn == i_cld_pc2) THEN
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            q_forcing(i,j,k)   = q_star_scm(i,j,k) - q_earliest(i,j,k)
            qcl_forcing(i,j,k) = qcl_inc(i,j,k)    - qcl_earliest(i,j,k)
            t_forcing(i,j,k)   = t_inc_scm(i,j,k)  - t_earliest(i,j,k)
            p_forcing(i,j,k)   = 0.0
          END DO
        END DO
      END DO
    END IF  ! i_cld_pc2

    !-------------------------------------------------------------------------
    ! Increments need to be converted to increment+value for t, q, qcl and qcf
    !-------------------------------------------------------------------------

    !  Before going into physics2, t_inc_scm, q_inc, qcl_inc and qcf_inc are
    !  converted to increment plus value whilst the winds stay as increments



    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length

          t_inc_scm(i,j,k) = t_inc_scm(i,j,k) + t(i,j,k)
          q_star_scm(i,j,k) = q_star_scm(i,j,k) + q(i,j,k)
          qcl_star(i,j,k)   = qcl_inc(i,j,k)    + qcl(i,j,k)
          qcf_star(i,j,k)   = qcf_inc(i,j,k)    + qcf(i,j,k)

          ! At present qcf2_inc etc are not used as there is
          ! no forcing option.

          qcf2_star(i,j,k)   = qcf2_star(i,j,k)   + qcf2(i,j,k)
          qrain_star(i,j,k)  = qrain_star(i,j,k)  + qrain(i,j,k)
          qgraup_star(i,j,k) = qgraup_star(i,j,k) + qgraup(i,j,k)

          cf_star (i,j,k) = cf_star (i,j,k) + cf(i,j,k)
          cfl_star(i,j,k) = cfl_star(i,j,k) + cfl(i,j,k)
          cff_star(i,j,k) = cff_star(i,j,k) + cff(i,j,k)
        END DO
      END DO
    END DO


    IF (l_qpos_for) THEN
      ! Check that forcing has not caused q to < qlimit

      ! NOTE: This has no equivalent in the full UM. It is justified
      !       as the forcing is decoupled from the model and so in the
      !       SCM there is no other mechanism to prevent a q < qlimit

      WRITE(sdum1,'(ES9.2)') qlimit
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length

            IF (q_star_scm(i,j,k) < qlimit) THEN
              ! The increment to q from the forcing routine will cause
              ! q to below qlimit. Don't allow this and output the
              ! required q to stop it going below qlimit
              dq_qpos_for(i,j,k) = qlimit - q_star_scm(i,j,k)
              q_star_scm(i,j,k)  = qlimit

              WRITE(sdum0,'(I3)') k

              CALL scm_message                                        &
                 ( 'q < '//TRIM(ADJUSTL(sdum1))//                     &
                   ': q reset on level '//TRIM(ADJUSTL(sdum0)) )
            ELSE
              dq_qpos_for(i,j,k) = 0.0
            END IF
          END DO      ! i
        END DO      ! j
      END DO      ! k

      IF (L_SCMDiags(SCMDiag_forc) .OR.                               &
          L_SCMDiags(SCMDiag_incs)) THEN

        CALL scmoutput                                                &
           ( dq_qpos_for, 'dq_qpos_for'                               &
           , 'Specific humidity adjustment to LS forcing to '//       &
             'maintain qlimit','kg/kg'                                &
           , t_avg, d_wet, default_streams, '', routinename )

      END IF ! l_scmdiags

    END IF ! l_qpos_for


    ! Convert t_inc_scm back to theta_star for call to Atmos_Physics2.
    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length

          theta_star(i,j,k) = t_inc_scm(i,j,k) / exner_theta_levels(i,j,k)

       END DO
      END DO
    END DO

    !-----------------------------------------------------------------------
    ! theta_star, q_star_scm, qcl_star and qcf_star now all store
    ! Increment + Whole Values
    !-----------------------------------------------------------------------



      ! Convert X_star to mixing ratios
    ! DEPENDS ON: q_to_mix
    CALL q_to_mix                                                           &
      ( row_length, rows                                                    &
      , halo_i, halo_j                                                      &
      , q_star(:,:,:,1:), qcl_star(:,:,1:), qcf_star(:,:,1:)                &
      , qcf2_star(:,:,1:), qrain_star(:,:,1:), qgraup_star(:,:,1:)          &
      , l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup                               &
      , mix_v_star, mix_cl_star, mix_cf_star                                &
      , mix_cf2_star, mix_rain_star, mix_graup_star )

    !-----------------------------------------------------------------------------
    ! VARIABLE VALUES:
    !   qcl_star   = liq. water       (at timelevel n)
    !              + liq. water incs  (atmos_physics1 & forcing)
    !   t_inc_scm  = temp.            (at timelevel n)
    !              + temp. incs       (atmos_physics1 & forcing)
    !   q_star_scm = vapour           (at timelevel n)
    !              + vapour incs      (atmos_physics1 & forcing)
    !   theta_star = pot. temp.       (at timelevel n)
    !              + pot. temp. incs  (atmos_physics1 & forcing)
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    ! Call physics2
    !-----------------------------------------------------------------------------
    ! ensure this is always 0 before ap2 call
    zlcl_mixed = 0.0
    l_emcorr_opt = .TRUE.


    DO  k=1, model_levels
      DO j=1, rows
        DO i=1, row_length 
          theta_eg(i,j,k)=theta(i,j,k)
        END DO
      END DO
    END DO

    ! DEPENDS ON: atmos_physics2
    CALL atmos_physics2                                                      &
    ! Parallel variables
       ( row_length, rows, n_proc, n_procx, n_procy,rows, row_length         &
       , numcycles, cycleno                                                  &

    ! Model dimensions.
       , row_length, rows, rows, land_points, bl_levels                      &
       , st_levels, sm_levels, cloud_levels, land_ice_points                 &
       , soil_points, n_cca_lev, ntiles, tr_levels                           &
       , first_constant_r_rho_level, dim_cs1, dim_cs2                        &

    ! Model switches
       , l_dry, l_lbc_old                                                    &
       , lcal360                                                             &

    ! Model parameters
       , rhcrit(1:model_levels), co2_mmr                                     &
       , tr_vars, tr_ukca                                                    &

    ! In: coordinate information
       , unscl_dry_rho                                                       &
       , delta_lambda, delta_phi, dlambda_p, dphi_p                          &
       , wt_lambda_p, wt_lambda_u, wt_phi_p, wt_phi_v                        &
       , lat_rot_np, long_rot_np, f3_at_u                                    &

    ! In: time stepping information.
       , yearno, day, ihour, imin, time_sec                                  &

    ! River routing
       , row_length, row_length, xpa, xua, xva, ypa, yua, yva                &
       , g_p_field, g_r_field, a_steps_since_riv, river_row_length           &
       , river_rows, global_river_row_length, global_river_rows              &
       , river_vel, river_mcoef, i_river_vn                                  &
       , trivdir, trivseq, twatstor, inlandout_atm                           &

    ! Lake evaporation:
       , acc_lake_evap                                                       &

    ! Grid-to-grid river routing
       , r_area, slope, flowobs1, r_inext, r_jnext, r_land                   &
       , substore, surfstore, flowin, bflowin,                               &

    ! Diagnostics info
         stashwork3, stashwork5, stashwork8, stashwork9, stashwork19         &
       , stashwork26                                                         &

    ! SCM diagnostics
       , nSCMdpkgs, l_SCMdiags, conv_mode, l_emcorr_opt                      &

    ! In: data fields.
       , theta_eg, q, qcl, qcf, qrain, qgraup, qcf2                          &
       , rho, u, v, w, w_adv, p, p_star                                      &
       , exner_rho_levels, exner_theta_levels                                &
       , land_sea_mask, p_theta_levels                                       &
    ! Mixing ratio prognostics; passed in separately to the q fields so
    ! that any convection scheme which works with mixing ratios
    ! can just use them instead of having to convert the q-fields back to
    ! mixing ratios, which would seem perverse!
       , mix_v, mix_cl, mix_cf, mix_cf2, mix_rain, mix_graup                 &

    ! Ancillary fields and fields needed to be kept from timestep to
    ! timestep
       , land_index, land_ice_index, soil_index, canopy_gb                   &
       , snodep, hcon, hcap, v_crit, v_wilt, v_sat, z0m_soil, sthf           &
       , sthu, sil_orog_land                                                 &
       , ho2r2_orog, sd_orog_land, di, ice_fract, u_0, v_0, u_0, v_0         &
       , cca_dp, cca_md, cca_sh                                              &
       , cca, iccb, icct, cclwp, ccw, lcbase                                 &
       , t_deep_soil, tsi, ti_n, ice_k_n, tstar                              &
       , z0msea, ice_fract_n, ice_thick_n, satcon, sathh, b_exp, smcl        &
       , t1_sd, q1_sd, zh, ddmfx, area_cloud_fraction, cf                    &
       , cfl, cff, ls_rain, ls_rainfrac, ls_snow, ls_graup                   &
       , micro_tends, totalppn, photosynth_act_rad, rad_hr                   &
       , soil_clay, soil_silt, soil_sand, dust_mrel1 ,dust_mrel2             &
       , dust_mrel3, dust_mrel4, dust_mrel5, dust_mrel6                      &
       , so2_high_level, so2_em, nh3_em, dms_em, soot_hilem                  &
       , soot_em, ocff_hilem, ocff_em, co2_emits, co2flux                    &
       , deep_flag, past_precip, past_conv_ht                                &

    ! InOut
       , theta_star, q_star_scm                                              &
       , qcl_star, qcf_star, qrain_star, qgraup_star, qcf2_star              &
       , cf_star, cfl_star, cff_star                                         &
       , u_inc_scm, v_inc_scm, w_inc_scm(:,:,1:model_levels)                 &
       , sum_eng_fluxes, sum_moist_flux                                      &

    ! InOut: tracer fields
       , aerosol, free_tracers, tracer_ukca                                  &
       , dust_div1, dust_div2, dust_div3, dust_div4, dust_div5, dust_div6    &
       , so2, dms, so4_aitken, so4_accu, so4_diss, nh3, soot_new             &
       , soot_aged, soot_cld, bmass_new, bmass_aged, bmass_cld               &
       , ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss, co2             &

    ! InOut: RIVERS
       , tot_surf_runoff, tot_sub_runoff                                     &

    ! Out: fields
       , ntml, cumulus, nbdsc, ntdsc                                         &
       , rhcpt,row_length, rows                                              &
       , zlcl_mixed                                                          &
    ! Additional variables for Land Surface Scheme
       , frac_typ, frac_disturb, canht, lai                                  &
       , canopy(1:land_points,1:ntiles)                                      &
       , catch(1:land_points,1:ntiles), catch_snow, snow_grnd, snow_tile     &
       , z0_tile, z0h_tile, tstar_tile, tsurf_elev_surft                     &
       , infil_tile(1:land_points,1:ntiles)                                  &
       , rgrain(1:land_points,1:ntiles)                                      &
       , cs, gs, co2_dim_row, co2_dim_len                                    &
       , asteps_since_triffid, timestep_number                               &
       , g_leaf_acc                                                          &
       , g_leaf_phen_acc, npp_ft_acc, resp_w_ft, resp_s_acc                  &
       , land_pts_trif, npft_trif, dolr, lw_down, sw_tile, fland_ctile       &
       , tstar_land, tstar_sea, tstar_sice_n, tstar_sice                     &
       , albsoil, cos_zenith_angle                                           &

    ! INOUT variables for TKE based turbulence schemes
       , e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu                        &
    ! OUT variable to store BL w-variance diagnostic
       , bl_w_var                                                            &
    ! IN/OUT convection prognostics
       , conv_prog_1, conv_prog_2, conv_prog_3, conv_prog_precip             &
       , ux_ccp, uy_ccp, um_ccp, g_ccp, h_ccp, riso_ccp, rdir_ccp            &
    ! Additional variables required for large-scale hydrology:
       , fexp, gamtot, ti_mean, ti_sig, fsat, fwetl, zw                      &
       , sthzw, a_fsat, c_fsat, a_fwet, c_fwet                               &

    ! JULES 2 prognostics (InOut)
       , snowdepth, rho_snow_grnd                                            &
       , nsnow, ds, sice, sliq, tsnowlayer, rho_snow, rgrainl                &

    ! FLake lake scheme prognostics (InOut)
       , lake_depth, lake_fetch, lake_t_mean, lake_t_mxl                     &
       , lake_t_ice, lake_h_mxl, lake_h_ice,  lake_shape                     &
       , lake_g_dt                                                           &

    ! Cariolle ozone
       , ozone_tracer                                                        &

    ! Additional screen-level variables
       , tscrndcl_ssi, tscrndcl_tile, tstbtrans                              &

    ! Variables required for COSP (out)
       , cosp_crain_3d, cosp_csnow_3d                                        &

    ! Prescribed surface forcing and roughness lengths
       , flux_e_scm, flux_h_scm, ustar_in, l_spec_z0, z0m_scm, z0h_scm       &

    ! Error information
       , error_code                                                          &
       )

    !--------------------------------------------------------------------------
    ! VARIABLE VALUES
    !   qcl_star   = liq. water      (at timelevel n)
    !              + liq. water incs (atmos_physics1, forcing & atmos_physics2)
    !   q_star_scm = vapour          (at timelevel n)
    !              + vapour incs     (atmos_physics1, forcing & atmos_physics2)
    !   theta_star = pot. temp.      (at timelevel n)
    !              + pot. temp. incs (atmos_physics1, forcing & atmos_physics2)
    !--------------------------------------------------------------------------


    IF (stats .OR. obs .OR. noforce) THEN

      ! rho and pressure out of sync with T, but can't update with IdealGL
      ! settings i.e. P=rho.R.T as this will cause a crash in the physics
      ! routines. I assume that Pressure and density need to be updated
      ! using exner_prime


              ! Update theta, q and winds
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length

            theta(i,j,k) = theta_star(i,j,k)

            !  Convert theta back to t
            t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)

            ! add wind increments to winds
            u(i,j,k) = u(i,j,k) + u_inc_scm(i,j,k)
            v(i,j,k) = v(i,j,k) + v_inc_scm(i,j,k)

            ! Vertical wind either prescribed for vertical advection
            ! forcing and updated after s_forcing or left at initial
            ! value.

          END DO
        END DO
      END DO

      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            q      (i,j,k) = q_star_scm  (i,j,k)
            qcl    (i,j,k) = qcl_star    (i,j,k)
            qcf    (i,j,k) = qcf_star    (i,j,k)
            qcf2   (i,j,k) = qcf2_star   (i,j,k)
            qrain  (i,j,k) = qrain_star  (i,j,k)
            qgraup (i,j,k) = qgraup_star (i,j,k)

            IF (i_cld_vn == i_cld_pc2) THEN
              cf (i,j,k) = cf_star (i,j,k)
              cfl(i,j,k) = cfl_star(i,j,k)
              cff(i,j,k) = cff_star(i,j,k)
            END IF  ! i_cld_pc2

          END DO
        END DO
      END DO

    END IF   ! stats .OR. obs .OR. noforce

    !--------------------------------------------------------------------------
    ! VARIABLE VALUES FOR THE NON-PC2 SITUATION.
    !   (FOR PC2 SEE NEXT SECTION OF CODE FOR ADDITIONAL TERMS)
    !
    !   qcl   = liq. water  (at timelevel n+1)
    !   q     = vapour      (at timelevel n+1)
    !   theta = pot. temp.  (at timelevel n+1)
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! PC2 section
    !--------------------------------------------------------------------------

    IF (i_cld_vn == i_cld_pc2) THEN

      ! The PC2 section needs to
      ! 1) increment condensation and cloud fractions due to the forcing
      ! 2) initiate cloud.
      ! 3) set area_cloud_fraction to bulk_cloud_fraction, cf, and correct
      !    theta.
      !
      ! 1. Calculate condensation increments which result from the forcing
      ! by using the homogeneous forcing approach

      !------------------------------------------------------------------------
      ! SCM PC2 Diagnostics Package
      !------------------------------------------------------------------------
      IF (main_diag_switch /= 0) THEN
        IF (l_scmdiags(scmdiag_pc2)) THEN

          ! save fields before incrementing with PC2 response to forcing
       DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            q_earliest   (i,j,k) = q  (i,j,k)
            qcl_earliest (i,j,k) = qcl(i,j,k)
            t_earliest   (i,j,k) = t  (i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

        END IF ! l_scmdiags(scmdiag_pc2)
      END IF ! main_diag_switch /= 0

      CALL pc2_homog_plus_turb                                              &
        ( p_theta_levels(1,1,1:model_levels), model_levels                  &
        , timestep, t, cf(:,:,1:), cfl(:,:,1:), cff(:,:,1:)                 &
        , q(:,:,1:), qcl(:,:,1:), t_forcing                                 &
        , q_forcing, qcl_forcing, p_forcing                                 &
        , 0.0, 0.0, l_mr_pc2)

      ! 1a. We have already applied the forcing to q, qcl and t in the
      !     forcing section. PC2_homog_plus_turb has now done this again
      !     so we need to subtract off these increments. Just the
      !     condensation will therefore remain.

      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            q(i,j,k)   = q(i,j,k)   - q_forcing(i,j,k)
            qcl(i,j,k) = qcl(i,j,k) - qcl_forcing(i,j,k)
            t(i,j,k)   = t(i,j,k)   - t_forcing(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k

      IF (main_diag_switch /= 0) THEN
        IF (l_scmdiags(scmdiag_pc2)) THEN
      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length
            qcl_inc (i,j,k)   = qcl(i,j,k) - qcl_earliest (i,j,k)
            q_star_scm(i,j,k) = q(i,j,k)   - q_earliest   (i,j,k)
            t_inc_scm(i,j,k)  = t(i,j,k)   - t_earliest   (i,j,k)
            cf_work (i,j,k)   = cf (i,j,k) - cf_star (i,j,k)
            cfl_work(i,j,k)   = cfl(i,j,k) - cfl_star(i,j,k)
            cff_work(i,j,k)   = cff(i,j,k) - cff_star(i,j,k)
          END DO ! i
        END DO ! j
      END DO ! k
          CALL scmoutput(qcl_inc, 'dqcl_pc2forc'                              &
            , 'PC2 qcl increment response to forcing', 'kg/kg'                &
            , t_inst, d_wet, default_streams, '', routinename)

          CALL scmoutput(q_star_scm(:,:,1:), 'dq_pc2forc'                     &
            , 'PC2 q increment response to forcing', 'kg/kg'                  &
            , t_inst, d_wet, default_streams, '', routinename)

          CALL scmoutput(t_inc_scm, 'dt_pc2forc'                              &
            , 'PC2 T increment response to forcing', 'K'                      &
            , t_inst, d_all, default_streams, '', routinename)

          CALL scmoutput(cf_work, 'dbcf_pc2forc'                              &
            , 'PC2 bulk cf increment response to forcing', ''                 &
            , t_inst, d_wet, default_streams, '', routinename)

          CALL scmoutput(cfl_work, 'dcfl_pc2forc'                             &
            , 'PC2 cfl increment response to forcing', ''                     &
            , t_inst, d_wet, default_streams, '', routinename)

          CALL scmoutput(cff_work, 'dcff_pc2forc'                             &
            , 'PC2 cff increment response to forcing', ''                     &
            , t_inst, d_wet, default_streams, '', routinename)

        END IF ! l_scmdiags(scmdiag_pc2)
      END IF ! main_diag_switch /= 0

      ! 2. Initiate cloud if necessary and check for sensible values.
      !    Use a dummy argument to receive the increment information.
      !    NOTE: SCM increments are output within pc2_initiation_ctl

      CALL pc2_initiation_ctl                                               &
        ( row_length, rows,   zlcl_mixed                                    &
        , l_mr_pc2                                                          &
        , nscmdpkgs, l_scmdiags                                             &
        , t, q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:)                            &
        , cf(:,:,1:), cfl(:,:,1:), cff(:,:,1:)                              &
        , rhts, tlts, qtts, ptts                                            &
        , area_cloud_fraction, p, p_star, p_theta_levels(1,1,1)             &
        , iccb, cumulus                                                     &
        , rhcpt                                                             &
        , t_work, q_work, qcl_work, qcf_work, cf_work, cfl_work, cff_work)


      IF (i_cld_area == acf_off) THEN

        ! 3. Set area_cloud_fraction to bulk_cloud_fraction, cf, and update
        !    theta from new temperature

        DO k=1, model_levels
          DO j=1, rows
            DO i=1, row_length
              area_cloud_fraction(i,j,k) = cf(i,j,k)
              theta(i,j,k) = t(i,j,k)/exner_theta_levels(i,j,k)
            END DO ! i
          END DO ! j
        END DO ! k

      ELSE IF (i_cld_area == acf_brooks) THEN

        CALL ls_acf_brooks(                                               &
            fv_cos_theta_latitude                                         &
          , cf(:,:,1:), cfl(:,:,1:), cff(:,:,1:)                          &
          , cumulus, area_cloud_fraction )

      END IF ! i_cld_area

    END IF  ! i_cld_pc2

    !--------------------------------------------------------------------------
    ! End of PC2 section.
    !
    ! Variables q, qcl, t, theta and cloud fractions
    ! are now all at timelevel n+1
    !--------------------------------------------------------------------------

    !-----------------------------------------------------------------------
    ! Check for negative q, qcl and qcf
    ! This would be done by q_pos in the full UM, l_qpos is specified
    ! in run_diffusion namelist
    !-----------------------------------------------------------------------
    IF (l_qpos) THEN
      ! Check q, qcl, and qcf. are not -ve
      dq_qpos   (:,:,:) = 0.0
      dqcl_qpos (:,:,:) = 0.0
      dqcf_qpos (:,:,:) = 0.0

      DO k=1, model_levels
        DO j=1, rows
          DO i=1, row_length

            IF (q(i,j,k) < 0.0) THEN
              WRITE(sdum0,'(I3)') k
              CALL scm_message                                       &
                 ( 'q < qlimit: q added on level '//                 &
                   TRIM(ADJUSTL(sdum0)))
              dq_qpos(i,j,k) = qlimit - q(i,j,k)
              q(i,j,k) = qlimit
            END IF

            IF (qcl(i,j,k) < 0.0) THEN
              WRITE(sdum0,'(I3)') k
              CALL scm_message                                       &
                 ( 'qcl < 0.0: qcl added on level '//                &
                   TRIM(ADJUSTL(sdum0)))
              dqcl_qpos(i,j,k) = - qcl(i,j,k)
              qcl(i,j,k) = 0.0
            END IF

            IF (qcf(i,j,k) < 0.0) THEN
              WRITE(sdum0,'(I3)') k
              CALL scm_message                                       &
                 ( 'qcf < 0.0: qcl added on level '//                &
                   TRIM(ADJUSTL(sdum0)))
              dqcf_qpos(i,j,k) = - qcf(i,j,k)
              qcf(i,j,k) = 0.0
            END IF

          END DO      ! i
        END DO      ! j
      END DO      ! k

      IF (main_diag_switch /= 0) THEN
        IF (l_scmdiags(scmdiag_forc) .OR.                            &
            l_scmdiags(scmdiag_incs)) THEN

          CALL scmoutput(dq_qpos,'dq_qpos'                           &
             , 'Q inc to prevent q < qlimit','kg/kg'                 &
             ,  t_avg,d_wet,default_streams,'',routinename)

          CALL scmoutput(dqcl_qpos,'dqcl_qpos'                       &
             , 'Qcl inc to prevent -ve qcl','kg/kg'                  &
             ,  t_avg,d_wet,default_streams,'',routinename)

          CALL scmoutput(dqcf_qpos,'dqcf_qpos'                       &
             , 'Qcf inc to prevent -ve qcf','kg/kg'                  &
             ,  t_avg,d_wet,default_streams,'',routinename)

        END IF        ! l_scmdiags
      END IF        ! main_diag_switch /= 0

    END IF        ! l_qpos


    !--------------------------------------------------------------------------
    ! Calculate some final diagnostics, and write the lot out if needs be.
    !--------------------------------------------------------------------------

    IF (main_diag_switch /= 0) THEN
      !
      !------------------------------------------------------------------------
      ! SCM PC2 Diagnostics Package
      !------------------------------------------------------------------------
      IF (l_scmdiags(scmdiag_pc2)) THEN

        CALL scmoutput(qcl(:,:,1:), 'qcl_n1_afterpc2'                       &
          , 'qcl at end of timestep after pc2', 'kg/kg'                     &
          , t_inst, d_wet, default_streams, '', routinename)

        CALL scmoutput(q(:,:,1:), 'q_n1_afterpc2'                           &
          , 'q at end of timestep after pc2', 'kg/kg'                       &
          , t_inst, d_wet, default_streams, '', routinename)

        CALL scmoutput(theta, 'th_n1_afterpc2'                              &
          , 'theta at end of timestep after pc2', 'K'                       &
          , t_inst, d_all, default_streams, '', routinename)

      END IF ! l_scmdiags(scmdiag_pc2)

      !------------------------------------------------------------------------
      ! SCM Increments Diagnostics Package
      !------------------------------------------------------------------------
      IF (l_scmdiags(scmdiag_incs)) THEN
        DO k=1, model_levels
          DO j=1, rows
            DO i=1, row_length
              t_totalinc(i,j,k) = t(i,j,k) - t_start(i,j,k)
              u_totalinc(i,j,k) = u(i,j,k) - u_start(i,j,k)
              v_totalinc(i,j,k) = v(i,j,k) - v_start(i,j,k)
              w_totalinc(i,j,k) = w(i,j,k) - w_start(i,j,k)

              q_totalinc   (i,j,k) = q  (i,j,k) - q_start  (i,j,k)
              qcl_totalinc (i,j,k) = qcl(i,j,k) - qcl_start(i,j,k)
              qcf_totalinc (i,j,k) = qcf(i,j,k) - qcf_start(i,j,k)
              cf_totalinc  (i,j,k) = cf (i,j,k) - cf_start (i,j,k)
              cfl_totalinc (i,j,k) = cfl(i,j,k) - cfl_start(i,j,k)
              cff_totalinc (i,j,k) = cff(i,j,k) - cff_start(i,j,k)

            END DO
          END DO
        END DO

        CALL scmoutput(t_totalinc, 'dt_total'                               &
          , 'Total increment to T', 'K'                                     &
          , t_avg, d_all, default_streams, '', routinename)

        CALL scmoutput(u_totalinc, 'du_total'                               &
          , 'Total increment to u', 'm/s'                                   &
          , t_avg, d_all, default_streams, '', routinename)

        CALL scmoutput(v_totalinc, 'dv_total'                               &
          , 'Total increment to v', 'm/s'                                   &
          , t_avg, d_all, default_streams, '', routinename)

        DO k=1, model_levels
          a2out(:,:,k) = w_totalinc(:,:,k)
        END DO

        CALL scmoutput(a2out, 'dw_total'                                    &
          , 'Total increment to w', 'm/s'                                   &
          , t_avg, d_all, default_streams, '', routinename)

        CALL scmoutput(q_totalinc, 'dq_total'                               &
          , 'Total increment to q', 'kg/kg'                                 &
          , t_avg, d_wet, default_streams, '', routinename)

        CALL scmoutput(qcl_totalinc, 'dqcl_total'                           &
          , 'Total increment to qcl', 'kg/kg'                               &
          , t_avg, d_wet, default_streams, '', routinename)

        CALL scmoutput(qcf_totalinc, 'dqcf_total'                           &
          , 'Total increment to qcf', 'kg/kg'                               &
          , t_avg, d_wet, default_streams, '', routinename)

        CALL scmoutput(cf_totalinc, 'dbcf_total'                            &
          , 'Total increment to bulk cf', ' '                               &
          , t_avg, d_wet, default_streams, '', routinename)

        CALL scmoutput(cfl_totalinc, 'dcfl_total'                           &
          , 'Total increment to cfl', ' '                                   &
          , t_avg, d_wet, default_streams, '', routinename)

        CALL scmoutput(cff_totalinc, 'dcff_total'                           &
          , 'Total increment to cff', ' '                                   &
          , t_avg, d_wet, default_streams, '', routinename)

        IF (geoforce) THEN
          CALL scmoutput(uinc_geo, 'du_geo'                                 &
            , 'Geostrophic forcing increment to u', 'm/s'                   &
            , t_avg, d_all, default_streams, '', routinename)

          CALL scmoutput(vinc_geo, 'dv_geo'                                 &
            , 'Geostrophic forcing increment to v', 'm/s'                   &
            , t_avg, d_all, default_streams, '', routinename)

          IF (ug_opt >=1) THEN
            geo_diag(:,:,:) = ug_scm(:,:,:)

            CALL scmoutput(geo_diag, 'u_g'                                  &
               , 'Zonal geostrophic wind', 'm/s'                            &
               , t_avg, d_all, default_streams, '', routinename)
          END IF

          IF (vg_opt >= 1) THEN
            geo_diag(:,:,:) = vg_scm(:,:,:)

            CALL scmoutput(geo_diag, 'v_g'                                  &
               , 'Meridional geostrophic wind', 'm/s'                       &
               , t_avg, d_all, default_streams, '', routinename)
          END IF

        END IF

      END IF ! l_SCMDiags(SCMDiag_incs)


       ! Store some diagnostic
      ! DEPENDS ON: dgnstcs_scm_main
      CALL dgnstcs_scm_main                                                 &
        ( row_length, rows, land_points, tr_levels, tr_vars, sm_levels      &
        , nsmax, st_levels, ntype                                           &
        , rho, timestep, u, v, t                                            &
        , theta, q(:,:,1:), qcl(:,:,1:), qcf(:,:,1:), cf(:,:,1:), cca, ccw  &
        , t_deep_soil, p_star, tstar, smc, canopy_gb, snodep, zh            &
        , z0msea, smcl, sthu, sthf, gs                                      &
        , nsnow, ds, sice, sliq, tsnowlayer, rho_snow, rgrainl,rgrain       &
        , lw_incs, photosynth_act_rad                                       &
        , tstar_tile, aerosol(:,:,1:), free_tracers                         &
        , p_theta_levels(:,:,1:), p, iccb, icct                             &
        , w, w_adv, area_cloud_fraction                                     &
        , cf(:,:,1:), cfl(:,:,1:), cff(:,:,1:), cclwp, nscmdpkgs, l_scmdiags)

      ! Initialise the diagnostic output files if this is the
      ! end of the first timestep for which the system was on
      IF (scmop%first_pass) THEN
        scmop%first_pass = .FALSE.

        ! The list of diagnostics should be finalised now, so
        ! dump_streams_init can be called.
        ! DEPENDS ON: dump_streams_init
        CALL dump_streams_init                                              &
          ! (InOut)
          ( scmop                                                           &
          ! (In)
          , row_length, rows                                                &
          , bl_levels, cloud_levels                                         &
          , ozone_levels, st_levels, sm_levels, ntiles                      &
          , year_init, month_init, day_init                                 &
          , hour_init, min_init, sec_init                                   &
          , timestep, ndayin,nminin, nsecin, sec_day                        &
          , ndayin*full_daysteps+nstepsin                                   &
          , a_sw_radstep_prog, a_sw_radstep_diag                            &
          , z_top_of_model, first_constant_r_rho_level                      &
          , eta_theta, eta_rho, orog, r_theta_levels                        &
          , r_rho_levels, netcdf_chunksize )

        !=======================================================
        ! ADD REMINDER TO CHECK NFOR IS CORRECT
        ! Added here it doesn't get push off screen
        ! and missed by user
        !=======================================================

        WRITE(umMessage,'(A)')                                        newline//&
          " |========================================"//              newline//&
          " | Namelist NFOR = " // TRIM(ADJUSTL(nfor_str)) //         newline//&
          " | NOTE: INCORRECT SPECIFICATION OF NFOR"//                newline//&
          " |       WILL PRODUCE UNINTENDED RESULTS."//               newline//&
          " |========================================"//              newline//&
          newline
        CALL umPrint(umMessage,src='scm_main')


      END IF

      ! Write output diagnostics to file(s)
      ! DEPENDS ON: dump_streams
      CALL dump_streams                                                     &
        ( scmop, day, time_sec, row_length, rows, model_levels, dayno_init  &
        , INT(timestep_number*timestep), site, lat, long                    &
        , time_initi, year, lcal360 )

    END IF ! main_diag_switch /= 0

    IF (l_ts_log) THEN
      CALL scm_message('Complete')
    END IF

  END DO   ! timestep_number

END DO   ! daycount

! Close the output files
! DEPENDS ON: dump_streams_end
CALL dump_streams_end(scmop)


DEALLOCATE( arcl )                  ! Aerosol climatology array for NWP
DEALLOCATE( true_latitude )
DEALLOCATE( true_longitude )
DEALLOCATE( cos_theta_latitude )
DEALLOCATE( sec_theta_latitude )
DEALLOCATE( sin_theta_latitude )
DEALLOCATE( cos_theta_longitude )
DEALLOCATE( sin_theta_longitude )
DEALLOCATE( FV_cos_theta_latitude )
DEALLOCATE( r_rho_levels )
DEALLOCATE( r_theta_levels )
DEALLOCATE( r_layer_centres )
DEALLOCATE( r_layer_boundaries )
DEALLOCATE( d_layer )
DEALLOCATE( eta_rho_levels )
DEALLOCATE( eta_theta_levels )
DEALLOCATE( f3_at_u )

DEALLOCATE( resdump )

CALL dealloc_common_scm
CALL dealloc_forcing

WRITE(umMessage,'(A)')                        &
  ' ------------------------'//      newline//&
  ' Run completed'//                 newline//&
  ' ------------------------'
CALL umPrint(umMessage,src='scm_main')

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE scm_main
