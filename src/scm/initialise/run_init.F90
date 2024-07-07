! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! SUBROUTINE Run_Init
!
! Purpose:
!   Called by scm_main (Single Column Model main routine) to do
!   initialisations.
!
! Code Description:
!    Language - FORTRAN 90
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Single Column Model
!
!=====================================================================
!     OPTIONS TO SET INITIAL PROFILES
!=====================================================================
! (i)   Observational large scale forcing (OBS=TRUE of namelist LOGIC)
!         Initial data is then from namelist INPROF
! (ii)  Statistical large scale forcing (STATS=TRUE of namelist LOGIC)
!         Initial data can either be derived from climate datasets
!         using subroutine INITSTAT or set from namelist
!         INPROF (set ALTDAT=TRUE in namelist LOGIC)
! (iii) No large-scale forcing initial data is set fron namelist
!         INPROF
! (iv)  Continuation from previous run stored on tape
!         (Set TAPEIN=TRUE in namelist LOGIC).  All other initial data
!         is overwritten
!=====================================================================

SUBROUTINE run_init                                                     &
  ! (In)
  ( row_length, rows, land_pts, land_sea_mask, nfor, nbl_levs           &
  , nsoilt_levs, nsoilm_levs, ntiles, nice, nice_use, ntrop, n_cca_lev  &
  , dayno_init, ichgf                                                   &
  , sec_day, rhcrit                                                     &
  , lcal360, l_snow_albedo                                              &
  ! (InOut)
  , rho, ti, smcl                                                       &
  ! (Out)
  , iseed, dayno_wint, iccb, icct, cca, cf, cfl, cff, t, q, qcl         &
  , qcf, theta_star, pstar, tstar, u, v, w, w_adv, z0msea, zh, t_deep_soil&
  , canopy_gb, tsi, smc, sthf, sthu, snodep, catch, infil_tile, z0_tile &
  , z0h_tile_bare, catch_snow, exner_rho_levels, exner_theta_levels, deltap, p&
  , p_theta_levels, rp, rp_theta, ch_tls, ch_qls, ch_uls, ch_vls, ch_wls&
  , ch_flux_e, ch_flux_h, ch_tstar_forcing, ch_ustar_forcing            &
  , ch_t_bg, ch_q_bg, ch_u_bg, ch_v_bg, ch_w_bg, ch_ug, ch_vg           &
  , flux_e_scm, flux_h_scm, ustar_in, t_inc_scm                         &
  , u_inc_scm, v_inc_scm, w_inc_scm, q_star_scm, ug_scm, vg_scm         &
  , t_bg_scm, q_bg_scm, u_bg_scm, v_bg_scm, w_bg_scm                    &
  , tls, qls, uls, vls, wls, dap1, dap2, dap3, dab1, dab2, dab3         &
  , b_exp, hcap, hcon, satcon, sathh, v_sat, v_wilt, v_crit, z0m_soil, atime&
  , btime, alfada, dbara, dgrada, pa, tbara, tgrada, tsda, vnbara, vnsda&
  , vpbara, wbara, wsda, alfadb, dbarb, dgradb, pb, tbarb, tgradb, tsdb &
  , vnbarb, vnsdb, vpbarb, wbarb, wsdb                                  &
  , cca_dp, cca_md, cca_sh, bl_w_var )


!USE relevant JULES routines
USE jules_init_mod,         ONLY: jules_init
USE freeze_soil_mod,        ONLY: freeze_soil
USE sparm_mod,              ONLY: sparm
USE infiltration_rate_mod,  ONLY: infiltration_rate
USE tilepts_mod,            ONLY: tilepts

USE jules_surface_types_mod
USE water_constants_mod, ONLY: rho_water
USE planet_constants_mod, ONLY: r, cp, kappa, pref, g
USE missing_data_mod,    ONLY: rmdi, imdi
USE level_heights_mod,   ONLY: r_theta_levels
USE nlsizes_namelist_mod, ONLY: model_levels
USE gen_phys_inputs_mod, ONLY: l_mr_physics

USE s_main_force, ONLY:                                                 &
   year_init, tapeday_init, soil_type, iccbi, iccti, nml_inprof_thetal  &
 , smi_opt, ndayin, resdump_days, runno_in, runno_out, timestep         &
 , ug, vg, canopy_gbi, ccai, smci, snodepi, t_deep_soili, tstari, ui, vi&
 , wi, z0mseai, frac_typ, canht, lai, fsmc, smcli, sth, altdat, geoforce&
 , geoinit, noforce, obs, obs_surf, stats, tapein, tapeout, prindump_obs&
 , exname_in, exname_out, t_bg, q_bg, u_bg, v_bg, w_bg                  &
 , theta, p_in, flux_e, flux_h, tstar_forcing, ustar_forcing, qi, t_inc &
 , q_star, u_inc, v_inc, w_inc, lat, long                               &
 , ccai_dp, ccai_md, ccai_sh                                            &
 , bl_w_var_i                                                           &
 , scm_clapp_levs      => clapp_levs                                    &
 , scm_sathh_levs      => sathh_levs                                    &
 , scm_hcap_levs       => hcap_levs                                     &
 , scm_hcon_levs       => hcon_levs                                     &
 , scm_satcon_levs     => satcon_levs                                   &
 , scm_smvccl_levs     => smvccl_levs                                   &
 , scm_smvcwt_levs     => smvcwt_levs                                   &
 , scm_smvcst_levs     => smvcst_levs

USE jules_surface_mod,   ONLY: l_aggregate

USE jules_soil_mod, ONLY: dzsoil_jules => dzsoil

! These parameters are dimensioned as (land_field,1,sm_levels or 0:sm_levels)
! in jules_init.
! Note the singleton 2nd dimension is in preparation for soil
! tiling. Due to soil tiling being a very large change, adding the soiltile
! dimension as a hard-coded singleton is being done as an intermediate step
USE p_s_parms,  ONLY:                                                   &
    jules_clapp_levs  => bexp_soilt                                     &
  , jules_sathh_levs  => sathh_soilt                                    &
  , jules_hcap_levs   => hcap_soilt                                     &
  , jules_hcon_levs   => hcon_soilt                                     &
  , jules_satcon_levs => satcon_soilt                                   &
  , jules_smvccl_levs => smvccl_soilt                                   &
  , jules_smvcwt_levs => smvcwt_soilt                                   &
  , jules_smvcst_levs => smvcst_soilt

USE scm_utils, ONLY:                                                    &
    old_nml, resdump
USE dumpinit_mod, ONLY: dumpinit

USE init_soil_mod, ONLY: nsoilp, b_exp_typ, v_crit_typ, v_sat_typ,      &
                         v_wilt_typ, hcap_typ, hcon_typ, sathh_typ,     &
                         satcon_typ

USE ereport_mod, ONLY: ereport
USE Control_Max_Sizes
USE umPrintMgr, ONLY:                                                   &
    umPrint,                                                            &
    umMessage
USE atm_fields_bounds_mod, ONLY: tdims
USE random_num_gen

USE errormessagelength_mod, ONLY: errormessagelength

USE parkind1, ONLY: jpim, jprb
USE yomhook,  ONLY: lhook, dr_hook

IMPLICIT NONE


!-----------------------------------------------------------------------------
! Arguments with INTENT(In)
!-----------------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                  &
  row_length           &! Leading x dimension of SCM arrays.
, rows                 &! Leading y dimension of SCM arrays.
, land_pts             &! Number of land points to be processed.
, nfor                 &! Number terms for observational forcing
, nbl_levs             &! Number of Boundary layer levels
, nsoilt_levs          &! Number of soil temperature levels
, nsoilm_levs          &! Number of soil moisture levels
, ntiles               &! Number of surface tiles
, nice                 &! Number of sea-ice categories
, nice_use             &! Number of sea-ice categories used fully
, ntrop                &! Max number of levels in the troposphere
, n_cca_lev             ! No of levels for cca

LOGICAL, INTENT(IN) :: land_sea_mask(row_length,rows)  ! True if land point


!-----------------------------------------------------------------------------

INTEGER, INTENT(IN) ::                                                  &
  dayno_init               &! Initial day in year
, ichgf                     ! No. of timesteps between change
                            ! in observational forcing

INTEGER, INTENT(IN) ::                                                  &
  sec_day

REAL, INTENT(IN) ::                                                     &
  rhcrit(model_levels)      ! Critical humidity for cloud formation.

! COMMENTS REFER TO TRUE STATUS
!------------------------------
LOGICAL, INTENT(IN) ::                                                  &
  lcal360                  &! Use 360 day year
, l_snow_albedo

INTEGER ::                                                              &
  land_index      (row_length*rows)

REAL, INTENT(OUT) ::                                                    &
  t_bg_scm(row_length,rows,model_levels)                                &
, u_bg_scm(row_length,rows,model_levels)                                &
, v_bg_scm(row_length,rows,model_levels)                                &
, w_bg_scm(row_length,rows,0:model_levels)                              &
, q_bg_scm(row_length,rows,model_levels)

!-----------------------------------------------------------------------------
! Arguments with INTENT(InOut)
!-----------------------------------------------------------------------------

REAL, INTENT(INOUT) ::                                                  &
  rho(row_length,rows,model_levels)    &!
, smcl(land_pts,nsoilm_levs)            ! Soil moisture content
                                        ! in layers (Kg/m2)

REAL, INTENT(INOUT) ::                                                  &
  ti(row_length,rows,model_levels)      ! Initial temperature profile (K)

!-----------------------------------------------------------------------------
! Arguments with INTENT(Out)
!-----------------------------------------------------------------------------
  ! Random generator variables
INTEGER, INTENT(OUT) ::                                                 &
  iseed                  ! Seed for random number generator

INTEGER, INTENT(OUT) ::                                                 &
  dayno_wint             ! Day number relative to winter solstice

INTEGER, INTENT(OUT) ::                                                 &
  iccb(row_length,rows) &! Convective cloud base
, icct(row_length,rows)  ! Convective cloud top

REAL, INTENT(OUT) ::                                                    &
  cca(row_length,rows,n_cca_lev)            &! Convective cloud amount
, cca_dp(row_length,rows,n_cca_lev)         &! Deep Conv cloud amount
, cca_md(row_length,rows,n_cca_lev)         &! Mid-level Conv cloud amount
, cca_sh(row_length,rows,n_cca_lev)         &! Shallow Conv cloud amount
, cf(row_length,rows,model_levels)          &! layer cloud amount
                                             ! (decimal fraction)
, cfl(row_length,rows,model_levels)         &! liquid layer cloud amount
, cff(row_length,rows,model_levels)          ! frozen layer cloud amount

REAL, INTENT(OUT) ::                                                    &
  bl_w_var(row_length,rows,model_levels)     ! w-variance for mixed phase

REAL, INTENT(OUT) ::                                                    &
  t(row_length,rows,model_levels)           &! Temperature(K)
, q(row_length,rows,model_levels)           &! Specific humidity (kg/kg)
, qcl(row_length,rows,model_levels)         &! Cloud water content(kg/kg)
, qcf(row_length,rows,model_levels)         &! Cloud ice content (kg/kg)
, theta_star(row_length,rows,model_levels)  &! Potential temp. increments(K)
, pstar(row_length,rows)                    &! Surface pressure (Pa)
, tstar(row_length,rows)                     ! Surface temperature (K)

REAL, INTENT(OUT) ::                                                    &
  u(row_length,rows,model_levels)           &! Zonal wind (m/s)
, v(row_length,rows,model_levels)           &! Meridional wind (m/s)
, w(row_length,rows,0:model_levels)

REAL, INTENT(OUT) ::                                                    &
  w_adv(row_length,rows,0:model_levels)                                 &
, z0msea(row_length,rows)                   &! Sea surface roughness length
, zh(row_length, rows)                       ! Height above surface of top
                                             ! of boundary layer (m)

REAL, INTENT(OUT) ::                                                    &
  t_deep_soil(land_pts,nsoilt_levs)   ! Deep soil temperatures (K)
                                      ! top level not included,=surface

REAL, INTENT(OUT) ::                                                    &
  canopy_gb(land_pts)         &! Canopy water content (kg/m2)
, tsi(row_length,rows)        &! Temperature of sea-ice
, smc(land_pts)               &! Soil moisture content(Kg/m^2)
, sthf(land_pts,nsoilm_levs)  &! Frozen soil moisture content of each
                               ! layer as a fraction of saturation (kg/m^2)
, sthu(land_pts,nsoilm_levs)  &! Unfrozen soil moisture content of each
                               ! layer as a fraction of saturation (kg/m^2)
, snodep(row_length,rows)      ! Snow depth (kg/m^2)

REAL, INTENT(OUT) ::                                                    &
  catch(row_length*rows,ntiles)      &! Surface/canopy water capacity of
                                      ! snow-free land tiles (kg/m2)
, infil_tile(row_length*rows,ntiles) &! Maximum surface infiltration rate
                                      ! for each tile (kg/m2/s)
, z0_tile(row_length*rows,ntiles)    &! Roughness length for each tile (m)
, z0h_tile_bare(row_length*rows,ntiles) &! Thermal roughness length
                                      ! for each tile (m)
, catch_snow(row_length*rows)         ! Snow capacity for NLT tile (kg/m^2)

REAL, INTENT(OUT) ::                                                    &
  exner_rho_levels(row_length,rows,model_levels+1)                      &
, exner_theta_levels(row_length,rows,tdims%k_start:tdims%k_end)         &
, deltap(row_length,rows,model_levels)              &! Layer Thickness
, p(row_length,rows,model_levels+1)                 &! Pressure rho levels
, p_theta_levels(row_length,rows,tdims%k_start:tdims%k_end)  &
                                                     ! Pressure theta levels
, rp(row_length,rows,model_levels+1)                &! 1/p on rho levels
, rp_theta(row_length,rows,model_levels)             ! 1/p on theta levels


! Rate of change of increments due to large-scale horizontal and
! vertical advection (per second) of ...
REAL, INTENT(OUT) ::                                                    &
  ch_tls(row_length,rows,nfor-1,model_levels)   &! ... Temp increment
, ch_qls(row_length,rows,nfor-1,model_levels)   &! ... Specific humidity
, ch_uls(row_length,rows,nfor-1,model_levels)   &! ... Zonal  wind
, ch_vls(row_length,rows,nfor-1,model_levels)   &! ... Merid. wind
, ch_wls(row_length,rows,nfor-1,0:model_levels)  ! ... Vert.  wind


! Rate of change of background fields for ...
REAL, INTENT(OUT) ::                                                    &
  ch_t_bg(row_length,rows,nfor-1,model_levels)  &! ... Temperature
, ch_q_bg(row_length,rows,nfor-1,model_levels)  &! ... Specific humidity
, ch_u_bg(row_length,rows,nfor-1,model_levels)  &! ... Zonal  wind
, ch_v_bg(row_length,rows,nfor-1,model_levels)  &! ... Merid. wind
, ch_w_bg(row_length,rows,nfor-1,0:model_levels) ! ... Vert.  wind


! Rate of change of geostrophic winds ...
REAL, INTENT(OUT) ::                                                    &
  ch_ug(row_length,rows,nfor-1,model_levels)    &! ... Zonal  wind
, ch_vg(row_length,rows,nfor-1,model_levels)     ! ... Merid. wind


! Rate of change of forcing of surface ...
REAL, INTENT(OUT) ::                                                    &
  ch_flux_e(row_length,rows,nfor-1)            &! ... flux_e
, ch_flux_h(row_length,rows,nfor-1)            &! ... flux_h
, ch_tstar_forcing(row_length,rows,nfor-1)     &! ... temperature
, ch_ustar_forcing(row_length,rows,nfor-1)      ! ... friction velocity


REAL, INTENT(OUT) ::                                                    &
  flux_e_scm(row_length,rows)                                           &
, flux_h_scm(row_length,rows)                                           &
, ustar_in(row_length,rows)


! Forcing increments (X/timestep)
REAL, INTENT(OUT) ::                                                    &
  q_star_scm (row_length,rows,model_levels)   &! Spec. humid. inc.
                                               ! (Kg/Kg)/timestep
, t_inc_scm  (row_length,rows,model_levels)   &! Temp. inc.  (K/timestep)
, u_inc_scm  (row_length,rows,model_levels)   &! Zonal  wind (m/s)/timestep
, v_inc_scm  (row_length,rows,model_levels)   &! Merid. wind (m/s)/timestep
, w_inc_scm  (row_length,rows,0:model_levels)  ! Vert.  wind (m/s)/timestep


! Geostrophics wind forcing
REAL, INTENT(OUT) ::                                                    &
  ug_scm     (row_length,rows,model_levels)   &! Geo. zonal  wind (m/s)
, vg_scm     (row_length,rows,model_levels)    ! Geo. merid. wind (m/s)


! Current LS forcing tendencies (X/day)
REAL, INTENT(OUT) ::                                                    &
  qls (row_length,rows,model_levels)   &! Spec. humid. inc. (Kg/Kg)/day
, tls (row_length,rows,model_levels)   &! Temp. inc.  (K/day)
, uls (row_length,rows,model_levels)   &! Zonal  wind (m/s)/day
, vls (row_length,rows,model_levels)   &! Merid. wind (m/s)/day
, wls (row_length,rows,0:model_levels)  ! Vert.  wind (m/s)/day



!---------------------------------------------------------------------
! Large scale observational forcing
!---------------------------------------------------------------------

! Variables for diagnostic output for observational forcing

REAL, INTENT(OUT) ::                           &! These don't appear
  dap1(row_length,rows,36,model_levels)        &! to do anything
, dap2(row_length,rows,36,model_levels)        &! except get initialised
, dap3(row_length,rows,36,nfor-1,model_levels)                          &
, dab1(row_length,rows,44)                                              &
, dab2(row_length,rows,44)                                              &
, dab3(row_length,rows,44,nfor-1)

REAL, INTENT(OUT) ::                                                    &
  b_exp  (land_pts)     &! Clapp-Hornberger exponent
, hcap   (land_pts)     &! Soil heat capacity
, hcon   (land_pts)     &! Soil thermal conductivity
, satcon (land_pts)     &! Saturated hydrological conductivity
, sathh  (land_pts)     &! Saturated soil water suction
, v_sat  (land_pts)     &! Vol. soil moisture content at saturation
, v_wilt (land_pts)     &! Vol. soil moisture content at wilting point
, v_crit (land_pts)     &! Vol. soil moisture content at critical point
, z0m_soil(land_pts)

!---------------------------------------------------------------------
!     Large scale statistical forcing
!---------------------------------------------------------------------

! Variable for statistical forcing

REAL, INTENT(OUT) ::                                                    &
  atime                                &! Constants for calculating annual
, btime                                 ! cycle

! Amplitude of seasonal variation of ...
REAL, INTENT(OUT) ::                                                    &
  alfada(row_length,rows)              &! ... Tuning Factor
, dbara(row_length,rows,model_levels)  &! ... Mean Dew pt. dep. (K)
, dgrada(row_length,rows,model_levels) &! ... Dew pt. dep. gradient (K/km)
, pa(row_length,rows, model_levels+1)  &! ... Pressure
, tbara(row_length,rows,model_levels)  &! ... Temperature (K)
, tgrada(row_length,rows,model_levels) &! ... Temperature gradient (K/km)
, tsda(row_length,rows,model_levels)   &! ... SD of temperature (K)
, vnbara(row_length,rows,model_levels) &! ... Velocity VN (m/s)
, vnsda(row_length,rows,model_levels)  &! ... SD of velocity VN (m/s)
, vpbara(row_length,rows,model_levels) &! ... Velocity VP (m/s)
, wbara(row_length,rows,ntrop)         &! ... Vert. Vel. (mb or HPa/s)
, wsda(row_length,rows,ntrop)           ! ... SD of Vert. Vel. (mb/s)

! Mean of seasonal variation of ...
REAL, INTENT(OUT) ::                                                    &
  alfadb(row_length,rows)              &! ... Tuning Factor
, dbarb(row_length,rows,model_levels)  &! ... Mean Dew pt. dep (K)
, dgradb(row_length,rows,model_levels) &! ... Dew pt. dep. gradient (K/km)
, pb(row_length,rows, model_levels+1)  &! ... Pressure
, tbarb(row_length,rows,model_levels)  &! ... Temperature (K)
, tgradb(row_length,rows,model_levels) &! ... Temperature gradient (K/km)
, tsdb(row_length,rows,model_levels)   &! ... SD of temperature (K)
, vnbarb(row_length,rows,model_levels) &! ... Velocity VN (m/s)
, vnsdb(row_length,rows,model_levels)  &! ... SD of velocity VN (m/s)
, vpbarb(row_length,rows,model_levels) &! ... Velocity VP (m/s)
, wbarb(row_length,rows,ntrop)         &! ... Vert. Vel. (mb or HPa/s)
, wsdb(row_length,rows,ntrop)           ! ... SD of Vert. Vel. (mb/s)



!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'RUN_INIT'
CHARACTER(LEN=errormessagelength)  :: cmessage ! Error message if Icode >0
INTEGER                     :: icode

INTEGER, PARAMETER ::                                                   &
  smc_spec = 0               &! Initial smcl specifed via namelist
, smc_fsmc = 1               &! Initial smcl calculated via fsmc
, smc_sth  = 2                ! Initial smcl calculated via sth

INTEGER ::                                                              &
  i, j, k, l, m              &! Loop counters
, runno                      &! Run number of experiment on tape
, nresdump                   &! Number of restart dumps
, andayy1, andayy2           &! To calculate andday using time2sec
, tapeday                    &! Tape year day
, tapedump_no                &! Number of dumps on tape
, dummy

INTEGER ::                                                              &
  ErrorStatus

INTEGER ::                                                              &
  ntml_tmp(row_length,rows)

INTEGER ::                                                              &
  tile_pts(ntype)            &! Number of land points which
                              ! include the nth surface type
, tile_index(land_pts,ntype)  ! Indices of land points which
                              ! include the nth surface type

REAL ::                                                                 &
  andayy                 &! Number of days in 1 year. (for one year effects)
, rccb(row_length,rows)  &! Conv. cloud base (real values) for DUMP purposes
, rcct(row_length,rows)   ! Conv. cloud top  (real values) for DUMP purposes

REAL ::                                                                 &
  fsmc_min(nsoilp)                                                      &
, fsmc_max(nsoilp)

REAL ::                                                                 &
  thl_to_sl              &! Top of mixed layer difference (K)
, sl_bl                   ! Mixed-layer value for s_L (K)


REAL ::                                                                 &
  factor                                                                &
, tstpfd  ! Timestep fraction of day

LOGICAL ::                                                              &
  cumulus_tmp(row_length,rows)

LOGICAL ::         &
  l_calc_exner     &
, l_calc_rho

CHARACTER(LEN=8) ::                                                     &
  exname                  ! Name of experiment on tape

!---------------------------------------------------------------------
!     Define site specific soil parameters and Initialise SWNOCZ
!---------------------------------------------------------------------

  ! Dr Hook
  !==============================
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb) :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

tstpfd = timestep/sec_day


IF (land_pts == 1) land_index(:) = 1

! The routine jules_init is part of the JULES repository. In the
! full um it sets the following variables from atm_fields_mod which
! the SCM doesn't have. So default values based on soil_type of
! the surface of SCM column is used to determine default values
! of the following JULEs parameters

! For the scm the call takes a few extra arguments due to differences in
! module from the main UM
CALL jules_init(land_index,                                             &
                land_pts, ntiles, nsoilm_levs                           &
                )

IF (land_pts >= 1) THEN
  DO i=1, land_pts
    jules_clapp_levs(i,1,:)  = b_exp_typ  (soil_type(i))
    jules_sathh_levs(i,1,:)  = sathh_typ  (soil_type(i))
    jules_hcap_levs(i,1,:)   = hcap_typ   (soil_type(i))
    jules_hcon_levs(i,1,:)   = hcon_typ   (soil_type(i))
    jules_satcon_levs(i,1,:) = satcon_typ (soil_type(i))
    jules_smvccl_levs(i,1,:) = v_crit_typ (soil_type(i))
    jules_smvcwt_levs(i,1,:) = v_wilt_typ (soil_type(i))
    jules_smvcst_levs(i,1,:) = v_sat_typ  (soil_type(i))
  END DO

  ! Jules has set the defaults, now see if SCM user has overrriden them
  ! in the SCM forcing namelists, test on first element, if unset
  ! it should be rmdi
  IF (scm_clapp_levs  (1,1) /= rmdi)                                    &
                      jules_clapp_levs  (:,1,:) = scm_clapp_levs  (:,:)

  IF (scm_sathh_levs  (1,1) /= rmdi)                                    &
                      jules_sathh_levs  (:,1,:) = scm_sathh_levs  (:,:)

  IF (scm_hcap_levs   (1,1) /= rmdi)                                    &
                      jules_hcap_levs   (:,1,:) = scm_hcap_levs   (:,:)

  IF (scm_hcon_levs   (1,1) /= rmdi)                                    &
                      jules_hcon_levs   (:,1,:) = scm_hcon_levs   (:,:)

  IF (scm_satcon_levs (1,1) /= rmdi)                                    &
                      jules_satcon_levs (:,1,:) = scm_satcon_levs (:,:)

  IF (scm_smvccl_levs (1,1) /= rmdi)                                    &
                      jules_smvccl_levs (:,1,:) = scm_smvccl_levs (:,:)

  IF (scm_smvcwt_levs (1,1) /= rmdi)                                    &
                      jules_smvcwt_levs (:,1,:) = scm_smvcwt_levs (:,:)

  IF (scm_smvcst_levs (1,1) /= rmdi)                                    &
                      jules_smvcst_levs (:,1,:) = scm_smvcst_levs (:,:)


  ! Set surface values passed to atmos_phyics2 via argument lists
  DO i=1, land_pts
    b_exp  (i) = jules_clapp_levs  (i,1,1)
    sathh  (i) = jules_sathh_levs  (i,1,1)
    hcap   (i) = jules_hcap_levs   (i,1,1)
    hcon   (i) = jules_hcon_levs   (i,1,0)
    satcon (i) = jules_satcon_levs (i,1,0)
    v_crit (i) = jules_smvccl_levs (i,1,1)
    v_wilt (i) = jules_smvcwt_levs (i,1,1)
    v_sat  (i) = jules_smvcst_levs (i,1,1)
    ! z0m_soil will be passed to sparm() so needs to be initialised here.
    ! In future, z0 for soil should come in from the JULES namelist.
    z0m_soil (i) = 0.0
  END DO

END IF ! Test on land pts

!---------------------------------------------------------------------
!     Calculate number of days in year
!---------------------------------------------------------------------

! DEPENDS ON: time2sec
CALL time2sec (year_init, 1, 1, 0, 0, 0, 0, 0, andayy1, dummy, lcal360)

! DEPENDS ON: time2sec
CALL time2sec (year_init+1, 1, 1, 0, 0, 0, 0, 0, andayy2, dummy, lcal360)
andayy = andayy2 - andayy1


IF (stats) THEN

  !---------------------------------------------------------------------
  !       Calculate day number relative to winter solstice
  !       (ie day 351 assuming 360 day calendar)
  !---------------------------------------------------------------------

  IF (dayno_init  <   351) THEN
    dayno_wint = dayno_init + 9
  ELSE
    dayno_wint = dayno_init - 351
  END IF

  !---------------------------------------------------------------------
  !       Derive initial data from climate datasets
  !---------------------------------------------------------------------
  DO k=1,model_levels+1
    DO j=1, rows
      DO i=1, row_length
        p(i,j,k) = p_in(i,j,k)
      END DO
    END DO
  END DO

  ! DEPENDS ON: initstat
  CALL initstat                                                         &
    ( row_length, rows, ntrop, andayy, dayno_wint, q, t                 &
    , lat, long, p_in, pa, pb, alfada, alfadb, tbara, tbarb, tsda, tsdb &
    , tgrada, tgradb, dbara, dbarb, dgrada, dgradb, vnbara, vnbarb, vnsda&
    , vnsdb, vpbara, vpbarb, wbara, wbarb, wsda, wsdb, atime, btime     &
    , p_theta_levels(:,:,1:) )

  DO k=1,model_levels
    DO j=1, rows
      DO i=1, row_length
        p_in(i,j,k) = p(i,j,k)
      END DO
    END DO
  END DO
  !---------------------------------------------------------------------
  !       Initialise random generator
  !---------------------------------------------------------------------


  CALL random_initial(iseed)
END IF                     ! stats

!---------------------------------------------------------------------
!     Set initial data from &INPROF
!---------------------------------------------------------------------

DO j=1, rows
  DO i=1, row_length

    IF (stats .OR. obs .OR. noforce .OR. geoforce) THEN
      !  Do for stats as well for now as we have no start data otherwise.
      DO k=1, model_levels
        u(i,j,k) = ui(i,j,k)
        v(i,j,k) = vi(i,j,k)
        w(i,j,k) = wi(i,j,k)
        w_adv(i,j,k) = w(i,j,k)
        t(i,j,k) = ti(i,j,k)
      END DO

      w(i,j,0)     = 0.0
      w_adv(i,j,0) = w(i,j,0)

      DO k = 1, model_levels
        q(i,j,k) = qi(i,j,k)
      END DO
    END IF                   ! (obs .or. noforce)

    IF (stats .AND. altdat) THEN
      DO k = 1, model_levels
        t(i,j,k) = ti(i,j,k)
      END DO
      DO k=1, model_levels
        q(i,j,k) = qi(i,j,k)
      END DO
    END IF                   ! (stats .and. altdat)

    IF (geoforce .AND. geoinit) THEN
      DO k=1, model_levels
        u(i,j,k) = ug(i,j,1,k)
        v(i,j,k) = vg(i,j,1,k)
      END DO
    END IF                   ! geoforce and geoinit

  END DO
END DO                     ! i

!---------------------------------------------------------------------
!     Calculate rates of change for large scale observational forcing
!---------------------------------------------------------------------

IF (obs) THEN

  IF (old_nml) THEN
    ! Reproduce old run format style
    ! Old style only had u_inc, v_inc, and w_inc. and treated
    ! them as background fields or forcing tendencies depending on
    ! whether l_windrlx was .TRUE. or .FALSE.

    ! So make wind _bg = wind _inc fields if relaxation is in old format
    ! so it will reproduce same results in restructured forcing code later
    ! on

    u_bg (:,:,:,:) = u_inc (:,:,:,:)
    v_bg (:,:,:,:) = v_inc (:,:,:,:)
    w_bg (:,:,:,:) = w_inc (:,:,:,:)

    ! This is not the preferred method, users with old namelists
    ! generated from earlier documentation wishing to use
    ! these in the restructured code should really rename u_inc, v_inc
    ! and w_inc to u_bg, v_bg and w_bg with units (m/s)

    ! In addition they should source forcing data for u_inc, v_inc and
    ! w_inc in the units of (m/s)/day

  END IF

  factor = ichgf * timestep

  DO l=1, (nfor-1)
    ch_t_bg(:,:,l,:) = (t_bg(:,:,l+1,:) - t_bg(:,:,l,:)) / factor
    ch_q_bg(:,:,l,:) = (q_bg(:,:,l+1,:) - q_bg(:,:,l,:)) / factor
    ch_u_bg(:,:,l,:) = (u_bg(:,:,l+1,:) - u_bg(:,:,l,:)) / factor
    ch_v_bg(:,:,l,:) = (v_bg(:,:,l+1,:) - v_bg(:,:,l,:)) / factor
    ch_w_bg(:,:,l,:) = (w_bg(:,:,l+1,:) - w_bg(:,:,l,:)) / factor

    ch_tls(:,:,l,:)  = (t_inc  (:,:,l+1,:) - t_inc  (:,:,l,:)) / factor
    ch_qls(:,:,l,:)  = (q_star (:,:,l+1,:) - q_star (:,:,l,:)) / factor
    ch_uls(:,:,l,:)  = (u_inc  (:,:,l+1,:) - u_inc  (:,:,l,:)) / factor
    ch_vls(:,:,l,:)  = (v_inc  (:,:,l+1,:) - v_inc  (:,:,l,:)) / factor
    ch_wls(:,:,l,:)  = (w_inc  (:,:,l+1,:) - w_inc  (:,:,l,:)) / factor
  END DO



  ! Initialise background fields
  t_bg_scm(:,:,:) = t_bg(:,:,1,:)
  q_bg_scm(:,:,:) = q_bg(:,:,1,:)
  u_bg_scm(:,:,:) = u_bg(:,:,1,:)
  v_bg_scm(:,:,:) = v_bg(:,:,1,:)
  w_bg_scm(:,:,:) = w_bg(:,:,1,:)

  ! Initialise LS forcing tendencies
  tls (:,:,:) = t_inc (:,:,1,:)
  qls (:,:,:) = q_star(:,:,1,:)
  uls (:,:,:) = u_inc (:,:,1,:)
  vls (:,:,:) = v_inc (:,:,1,:)
  wls (:,:,:) = w_inc (:,:,1,:)

  ! Set initial increments for timestep
  t_inc_scm (:,:,:) = tstpfd * tls(:,:,:)
  q_star_scm(:,:,:) = tstpfd * qls(:,:,:)
  u_inc_scm (:,:,:) = tstpfd * uls(:,:,:)
  v_inc_scm (:,:,:) = tstpfd * vls(:,:,:)
  w_inc_scm (:,:,:) = tstpfd * wls(:,:,:)


  IF (obs_surf) THEN

    DO l=1, (nfor-1)
      ch_flux_h(:,:,l) = (flux_h(:,:,l+1) - flux_h(:,:,l)) / factor
      ch_flux_e(:,:,l) = (flux_e(:,:,l+1) - flux_e(:,:,l)) / factor
      ch_tstar_forcing(:,:,l) = (tstar_forcing(:,:,l+1)                 &
                                   - tstar_forcing(:,:,l)) / factor
      ch_ustar_forcing(:,:,l) = (ustar_forcing(:,:,l+1)                 &
                                   - ustar_forcing(:,:,l)) / factor
    END DO

    flux_h_scm(:,:) = flux_h(:,:,1)
    flux_e_scm(:,:) = flux_e(:,:,1)
    tstar(:,:)      = tstar_forcing(:,:,1)
    ustar_in(:,:)   = ustar_forcing(:,:,1)

  END IF  ! obs_surf

END IF ! obs

IF (geoforce) THEN
  DO k = 1, model_levels
    DO l = 1, nfor - 1
      DO j = 1, rows
        DO i = 1, row_length
          ch_ug(i,j,l,k)  = (ug(i,j,l+1,k) - ug(i,j,l,k))               &
                          /             (ichgf * timestep)
          ch_vg(i,j,l,k)  = (vg(i,j,l+1,k) - vg(i,j,l,k))               &
                          /             (ichgf * timestep)
        END DO
      END DO
    END DO
  END DO

  ug_scm (:,:,:) = ug (:,:,1,:)
  vg_scm (:,:,:) = vg (:,:,1,:)

END IF


IF (tapein) THEN

  !---------------------------------------------------------------------
  !       Read tape data if required to carry on from previous run
  !---------------------------------------------------------------------

  READ (50) exname,runno,tapedump_no
  WRITE(umMessage,'(A17)')                                              &
      ' from tape header'
  CALL umPrint(umMessage,src='run_init')
  WRITE(umMessage,'(A30,A6)')                                           &
      ' expt. name is                ',exname
  CALL umPrint(umMessage,src='run_init')
  WRITE(umMessage,'(A30,i4)')                                           &
      ' run no. is                   ',runno
  CALL umPrint(umMessage,src='run_init')
  WRITE(umMessage,'(A30,i4)')                                           &
      ' no. of dumps on tape are     ',tapedump_no
  CALL umPrint(umMessage,src='run_init')

  !---------------------------------------------------------------------
  !       Check for correct data set
  !---------------------------------------------------------------------

  IF (exname  ==  exname_in .AND. runno  ==  runno_in) THEN

    !---------------------------------------------------------------------
    !         Look for correct day - tapeday_init input in namelist
    !         INDATA.
    !---------------------------------------------------------------------

    DO i=1, tapedump_no
      IF (stats) THEN
        READ (50) tapeday, resdump, iv, iy, drive
      ELSE IF (obs) THEN
        READ (50) tapeday, resdump
      END IF
      WRITE(umMessage,'(A16,I4,A22,I4)')                                &
        ' tape year day= ', tapeday,                                    &
        '      start year day= ', tapeday_init
      CALL umPrint(umMessage,src='run_init')

      IF (tapeday  ==  (tapeday_init-1) ) THEN
        GO TO 9999
      END IF

      !           If the end of the tape is reached and the specified day
      !           tapeday_init not found o/p error message and stop run.

      IF (i  ==  tapedump_no) THEN
        icode = 520
        WRITE(umMessage,'(A39,A6,A4,I3)')                               &
          ' Initial day not found on data set m20.',exname,'.run',runno
        CALL umPrint(umMessage,src='run_init')
        CALL umPrint('',src='run_init')


        CALL ereport(routinename, icode, cmessage)
      END IF
    END DO
    9999  CLOSE (50)

    !---------------------------------------------------------------------
    !         Read initial data from tape in DUMP format.
    !---------------------------------------------------------------------

    CALL dumpinit                                                       &
      ( row_length, rows, land_pts, nbl_levs                            &
      , nsoilt_levs, nsoilm_levs, n_cca_lev, land_sea_mask, u, v, w     &
      , t, theta(:,:,1:), q, qcl, qcf, cf, p, rho, t_deep_soil, smc     &
      , canopy_gb                                                       &
      , snodep, tstar, zh, z0msea, cca, cca_dp, cca_md, cca_sh          &
      , rccb, rcct, smcl, bl_w_var )
    ! Note: if you add more fields to the dump used here (array resdump),
    ! you need to increase the declared array size for resdump, set in the
    ! module scm_utils.
    ! You also need to modify the routine which saves the dump consistently;
    ! restart_dump (called from scm_main).
    ! You also need to modify the other call to dumpinit, from scm_main.


    !         Sometimes (when initial wind in dump is arbitrary) we want
    !         to reset the wind to geostrophic. To do this set geoinit to
    !         true in logic namelist

    DO j=1, rows
      DO i=1, row_length
        IF (geoinit .AND. geoforce) THEN
          DO k=1, model_levels
            u(i,j,k) = ug(i,j,1,k)
            v(i,j,k) = vg(i,j,1,k)
          END DO
        END IF
        iccb(i,j) = INT(rccb(i,j))
        icct(i,j) = INT(rcct(i,j))
        tsi(i,j) = tstar(i,j)
      END DO
    END DO

    !---------------------------------------------------------------------
    !         restore the state of the basic generator
    !---------------------------------------------------------------------

    CALL restore_random_state(drive,iv,iy)
  ELSE
    icode = 521
    WRITE(umMessage,'(A21,A6,A4,I3,A10)')                               &
      'Initial data set m20', exname_in,'.run',runno_in,' not found'
    CALL umPrint(umMessage,src='run_init')
    CALL umPrint('',src='run_init')

    CALL ereport(routinename, icode, cmessage)
  END IF                   ! (exname  ==  exname_in
                        ! .and. runno  ==  runno_in)


ELSE                      ! not tapein
  !---------------------------------------------------------------------
  !       Set initial values if no tape data to be used
  !---------------------------------------------------------------------
  !

  IF (land_pts > 0) THEN
    !-----------------------------------------------------------------------
    ! Call TILEPTS to initialise TILE_PTS and TILE_INDEX
    !-----------------------------------------------------------------------
    CALL tilepts (land_pts, frac_typ, tile_pts, tile_index)
    !-----------------------------------------------------------------------
    ! Initialise tiled and gridbox mean vegetation parameters
    !-----------------------------------------------------------------------
    WRITE(umMessage,*) 'RUN_INIT: CALLING SPARM'
    CALL umPrint(umMessage,src='run_init')

    CALL sparm(land_pts, ntiles, tile_pts, tile_index,                  &
               frac_typ, canht, lai, z0m_soil,                          &
               catch_snow, catch, z0_tile, z0h_tile_bare)

    CALL infiltration_rate(land_pts, ntiles, tile_pts, tile_index,      &
                           satcon, frac_typ, infil_tile)

  END IF

  DO k=1, model_levels+1
    DO j=1, rows
      DO i=1, row_length
        p(i,j,k) = p_in(i,j,k)
      END DO
    END DO
  END DO

  DO k=1, n_cca_lev
    DO j=1, rows
      DO i=1, row_length
        cca(i,j,k) = ccai(i,j)
        cca_dp(i,j,k) = ccai_dp(i,j)
        cca_md(i,j,k) = ccai_md(i,j)
        cca_sh(i,j,k) = ccai_sh(i,j)
      END DO
    END DO
  END DO

  DO j=1, rows
    DO i=1, row_length
      tstar  (i,j) = tstari  (i,j)
      tsi    (i,j) = tstar   (i,j)
      z0msea (i,j) = z0mseai (i,j)
      iccb   (i,j) = iccbi   (i,j)
      icct   (i,j) = iccti   (i,j)
      snodep (i,j) = snodepi (i,j) ! Includes snow depth on land/sea
    END DO
  END DO

  DO k = 1, model_levels
    DO j=1, rows
      DO i=1, row_length
        bl_w_var(i,j,k) = bl_w_var_i(i,j)
      END DO
    END DO
  END DO

  IF (land_pts > 0) THEN

    !--------------------------------------------------------------
    ! Initialise the canopy water
    !--------------------------------------------------------------
    DO i=1, land_pts
      IF (canopy_gbi(i) == rmdi) THEN
        icode = 522
        CALL umPrint('====================================================',&
            src='run_init')
        CALL umPrint('| Initial Gridbox Mean Canopy Water (canopy_gbi)   |',&
            src='run_init')
        CALL umPrint('| must be specified for each land point in         |',&
            src='run_init')
        CALL umPrint('| namelist.                                        |',&
            src='run_init')
        CALL umPrint('====================================================',&
            src='run_init')


        CALL ereport(routinename, icode, cmessage)

      ELSE
        canopy_gb(i) = canopy_gbi(i)
      END IF
    END DO

    !--------------------------------------------------------------
    ! Initialise the soil moisture content (SMC)
    !--------------------------------------------------------------
    DO i=1, land_pts
      IF (smci(i) == rmdi) THEN
        icode =523
        CALL umPrint('====================================================',&
            src='run_init')
        CALL umPrint('| Initial Soil Moisture Content (smci) must be     |',&
            src='run_init')
        CALL umPrint('| specified for each land point in namelist.       |',&
            src='run_init')
        CALL umPrint('====================================================',&
            src='run_init')

        CALL ereport(routinename, icode, cmessage)
      ELSE
        smc(i) = smci(i)
      END IF
    END DO

    !--------------------------------------------------------------
    ! Initialise the deep soil temperature
    !--------------------------------------------------------------
    DO k=1, nsoilm_levs
      DO i=1, land_pts
        IF (t_deep_soili(i,k) == rmdi) THEN
          icode=524
          CALL umPrint('===================================================='&
              ,src='run_init')
          CALL umPrint('| Initial Deep Soil Temperature (t_deep_soili)     |'&
              ,src='run_init')
          WRITE(umMessage,'(A15,I2,A18)')'| must specify ',nsoilt_levs, &
              ' soil temperatures','|'
          CALL umPrint(umMessage,src='run_init')
          CALL umPrint('| for each land point in namelist.                 |'&
              ,src='run_init')
          CALL umPrint('===================================================='&
              ,src='run_init')

          CALL ereport(routinename, icode, cmessage)
        ELSE
          t_deep_soil(i,k) = t_deep_soili(i,k)
        END IF
      END DO
    END DO               ! nsoilm_levs

    !--------------------------------------------------------------
    ! Initialise the soil moisture in root zone layers (SMCL)
    !--------------------------------------------------------------

    SELECT CASE (smi_opt)
    CASE (smc_spec)
      ! Soil moisture (smcli) to be specified by user namelist.

      ! Check for missing data
      DO k=1, nsoilm_levs
        DO i=1, land_pts
          IF (smcli(i,k) == rmdi) THEN
            icode=525
            CALL umPrint(                                               &
                '====================================================', &
                src='run_init')
            CALL umPrint(                                               &
                '| Initial Soil Moisture Content in layers (smcli)  |', &
                src='run_init')
            WRITE(umMessage,'(A25,I2,A21)')'| must specify values on ', &
                nsoilm_levs,' soil moisture levels','|'
            CALL umPrint(umMessage,src='run_init')
            CALL umPrint(                                               &
                '| for each land point in namelist.                 |', &
                src='run_init')
            CALL umPrint(                                               &
                '====================================================', &
                src='run_init')

            CALL ereport(routinename, icode, cmessage)
          ELSE
            smcl(i,k) = smcli(i,k)
          END IF
        END DO
      END DO

    CASE (smc_fsmc)
      ! Calculate max & min limits of FSMC
      DO i=1, nsoilp
        fsmc_min(i) = -  v_wilt_typ(i)                                  &
                      / (v_crit_typ(i) - v_wilt_typ(i))
        fsmc_max(i) =   (v_sat_typ(i)  - v_wilt_typ(i))                 &
                      / (v_crit_typ(i) - v_wilt_typ(i))
      END DO

      ! Soil moisture to be initialised by soil stress factor
      ! (FSMC) from namelist.


      DO i=1, land_pts

        ! Check for missing data
        IF (fsmc(i) == rmdi) THEN
          icode=526
          CALL umPrint(                                                 &
              '====================================================',   &
              src='run_init')
          CALL umPrint(                                                 &
              '| Soil moisture stress factor (fsmc) must be       |',   &
              src='run_init')
          CALL umPrint(                                                 &
              '| specified for each land point in namelist.       |',   &
              src='run_init')
          CALL umPrint(                                                 &
              '====================================================',   &
              src='run_init')
          CALL ereport(routinename, icode, cmessage)
        END IF

        ! Check fsmc is within range for given soil type
        IF ((fsmc(i) < fsmc_min(soil_type(i))) .OR.                     &
            (fsmc(i) > fsmc_max(soil_type(i)))) THEN
          icode=527
          CALL umPrint                                                  &
              ('====================================================',  &
              src='run_init')
          CALL umPrint                                                  &
              ('| Soil moisture stress factor (fsmc) is outside    |',  &
              src='run_init')
          CALL umPrint                                                  &
              ('| acceptable range. For the specified soil types,  |',  &
              src='run_init')
          CALL umPrint                                                  &
              ('| fsmc should be in ranges:                        |',  &
              src='run_init')
          CALL umPrint                                                  &
              ('|                                                  |',  &
              src='run_init')
          WRITE(umMessage,'(A19,F6.3,A10,F6.3)')'|  1) Ice        : ',  &
              fsmc_min(1),' < fsmc < ',                                 &
              fsmc_max(1),                 '|'
          CALL umPrint(umMessage,src='run_init')
          WRITE(umMessage,'(A19,F6.3,A10,F6.3)')'|  2) Clay       : ',  &
              fsmc_min(2),' < fsmc < ',                                 &
              fsmc_max(2),                 '|'
          CALL umPrint(umMessage,src='run_init')
          WRITE(umMessage,'(A19,F6.3,A10,F6.3)')'|  3) Loam       : ',  &
              fsmc_min(3),' < fsmc < ',                                 &
              fsmc_max(3),                 '|'
          CALL umPrint(umMessage,src='run_init')
          WRITE(umMessage,'(A19,F6.3,A10,F6.3)')'|  4) Loamy Sand : ',  &
              fsmc_min(4),' < fsmc < ',                                 &
              fsmc_max(4),                 '|'
          CALL umPrint(umMessage,src='run_init')
          CALL umPrint                                                  &
              ('|                                                  |',  &
              src='run_init')
          CALL umPrint                                                  &
              ('====================================================',  &
              src='run_init')

          CALL ereport(routinename, icode, cmessage)
        END IF
      END DO

      DO k=1, nsoilm_levs
        DO i=1, land_pts
          smcl(i,k) = rho_water*dzsoil_jules(k)                         &
                    * ( fsmc(i)*v_crit(i) + (1-fsmc(i))*v_wilt(i) )
        END DO
      END DO

    CASE (smc_sth)
      ! Soil moisture is to be initialised by total soil
      ! moisture as a fraction of saturation (STH).

      DO k=1, nsoilm_levs
        DO i=1, land_pts

          ! Check for missing data
          IF (sth(i,k) == rmdi) THEN
            icode=528
            CALL umPrint                                                &
                ('====================================================',&
                src='run_init')
            CALL umPrint                                                &
                ('| Total Soil Moisture (sth) as a fraction of       |',&
                src='run_init')
            WRITE(umMessage,'(A37,I2,A5,A1)')                           &
                '| saturation, must specify values on ', nsoilm_levs,   &
                ' soil','|'
            CALL umPrint(umMessage,src='run_init')
            CALL umPrint                                                &
                ('| moisture levels for each land point in namelist. |',&
                src='run_init')
            CALL umPrint                                                &
                ('====================================================',&
                src='run_init')

            CALL ereport(routinename, icode, cmessage)
          END IF

          ! Check sth range
          IF ((sth(i,k) < 0.0) .OR.                                     &
              (sth(i,k) > 1.0)) THEN
            icode=529
            CALL umPrint                                                &
                ('====================================================',&
                src='run_init')
            CALL umPrint                                                &
                ('| Specified values for Total Soil Moisture (sth)   |',&
                src='run_init')
            CALL umPrint                                                &
                ('| as a fraction of saturation, must be within the  |',&
                src='run_init')
            CALL umPrint                                                &
                ('| range: 0.0 ~ 1.0                                 |',&
                src='run_init')
            CALL umPrint                                                &
                ('====================================================',&
                src='run_init')

            CALL ereport(routinename, icode, cmessage)
          END IF

          smcl(i,k) = sth(i,k)*rho_water*v_sat(i) * dzsoil_jules(k)

        END DO
      END DO

    CASE (imdi)
      icode=530
      CALL umPrint                                                      &
          ('====================================================',      &
          src='run_init')
      CALL umPrint                                                      &
          ('| An initialization option for soil moisture       |',      &
          src='run_init')
      CALL umPrint                                                      &
          ('| content in layers (smi_opt) must be specified    |',      &
          src='run_init')
      CALL umPrint                                                      &
          ('| when running over land points                    |',      &
          src='run_init')
      CALL umPrint                                                      &
          ('====================================================',      &
          src='run_init')


      CALL ereport(routinename, icode, cmessage)
    END SELECT
  END IF                    ! land_pts

END IF                     ! tapein

!---------------------------------------------------------------------
!     Initialise the frozen and unfrozen soil moisture.
!---------------------------------------------------------------------
CALL freeze_soil(land_pts, nsoilm_levs, jules_clapp_levs(:,1,:),          &
                 dzsoil_jules, jules_sathh_levs(:,1,:), smcl, t_deep_soil,&
                 jules_smvcst_levs(:,1,:), sthu, sthf )

!---------------------------------------------------------------------
!     Calculate pressure, exner_theta_levels and pstar
!---------------------------------------------------------------------

! DEPENDS ON: calc_press
CALL calc_press                                                         &
! In
    ( rows, row_length, p, theta, q                                     &
    , l_calc_exner, l_calc_rho                                          &
! InOut
    , rho                                                               &
! Out
    , exner_theta_levels, exner_rho_levels, p_theta_levels, rp, rp_theta&
    , pstar )


!---------------------------------------------------------------------
!     Calculate DELTAP
!---------------------------------------------------------------------

DO k=1, model_levels
  DO j=1, rows
    DO i=1, row_length
      deltap(i,j,k) = p(i,j,k+1)-p(i,j,k)
    END DO
  END DO
END DO
!---------------------------------------------------------------------
!     Initialise cloud water (QCL,QCF)
!---------------------------------------------------------------------

IF (.NOT. tapein) THEN
  IF (NML_inprof_thetal == 0) THEN

    ! Standard, to input theta and qi as theta and qv
    DO k=1, model_levels
      ! DEPENDS ON: initqlcf
      CALL initqlcf                                                     &
         ( p_theta_levels(:,:,k), rhcrit, q(:,:,k), t(:,:,k)            &
         , model_levels, row_length, rows, cf(:,:,k)                    &
         , qcf(:,:,k), qcl(:,:,k), nbl_levs, k )
    END DO
    DO i=1, row_length
      DO j=1, rows
        DO  k=1, model_levels
          IF ( qcl(i,j,k)+qcf(i,j,k) > 0.0 ) THEN
            ! use simple water-content-dependent split
            cfl(i,j,k) = cf(i,j,k)*qcl(i,j,k)/(qcl(i,j,k)+qcf(i,j,k))
            cff(i,j,k) = 1.0 - cfl(i,j,k)
          ELSE
            cfl(i,j,k) = 0.0
            cff(i,j,k) = 0.0
          END IF
        END DO  ! k
      END DO  ! j
    END DO  ! i

  ELSE
    !-------------------------------------------------------
    ! Alternative, to input theta and qi as theta_l and qt
    ! where theta_l = (T -L*qcl/cp)/pi = theta - L*qcl/(cp*pi)
    !-------------------------------------------------------
    ! Set up fields ready to pass to ls_cld

    DO j=1, rows
      DO i=1, row_length
        ! Just initialise assuming NOT convection
        cumulus_tmp(i,j) = .FALSE.

        DO k=1, model_levels
          qcf(i,j,k) = 0.0
        END DO

        !-------------------------------------------------------
        ! Impose well-mixed SL(=T -L*qcl/cp+g*z/cp) where input
        ! theta_l is well-mixed in the BL
        ! because SL is the BL's conserved thermodynamic variable
        !-------------------------------------------------------
        IF (NML_inprof_thetal == 2) THEN

          DO k=1, model_levels
            ! Convert "T" back to "theta"=theta_l from namelist "theta"
            t(i,j,k) = t(i,j,k)/exner_theta_levels(i,j,k)
          END DO

          ! Find where t (currently =theta_l from namelist inprof)
          ! stops being well-mixed.
          k=2
          DO WHILE ( k < model_levels .AND.                             &
                     ABS(t(i,j,k)-t(i,j,1)) < 0.01)
            k=k+1
          END DO

          ntml_tmp(i,j) = k-1

          ! mixed layer s_L = mixed_layer_theta_l * pi(surf)
          sl_bl = t(i,j,1)*(pstar(i,j)/pref)**kappa

          DO k=1, ntml_tmp(i,j)
            ! Save theta_l*pi
            thl_to_sl = t(i,j,k)*exner_theta_levels(i,j,k)

            ! Convert s_L to T_L and
            ! set mixed layer T_L = mixed layer s_L - gz/cp
            t(i,j,k) = sl_bl                                            &
                     - g*(r_theta_levels(i,j,k)-r_theta_levels(i,j,0))/cp

            ! Calc change on converting from mixed theta_l to mixed s_L
            thl_to_sl = t(i,j,k) - thl_to_sl
          END DO

          ! For rest of profile, adjust by same amount as top of mixed
          ! layer, in order to preserve inversion jump
          DO k=ntml_tmp(i,j)+1, model_levels
            t(i,j,k) = t(i,j,k)*exner_theta_levels(i,j,k) + thl_to_sl
          END DO

        END IF  ! test on NML_inprof_thetal=2
      END DO
    END DO

    ! Arguments passed to ls_cld:
    !   t = T_L = T-L*qcl/cp
    !   q = q_w = q+qcl
    ! DEPENDS ON: ls_cld
    CALL ls_cld                                                               &
      ( p_theta_levels(1,1,1), rhcrit, model_levels, nbl_levs, row_length     &
      , rows, ntml_tmp, cumulus_tmp, l_mr_physics, t, cf, q, qcf, qcl         &
      , cfl, cff, ErrorStatus )

    DO k=1, model_levels
      DO j=1, rows
        DO i=1, row_length
          ! Reset initial profiles
          ti(i,j,k) = t(i,j,k)
          qi(i,j,k) = q(i,j,k)
        END DO
      END DO
    END DO

  END IF
END IF

!---------------------------------------------------------------------
!       Convert temperature to potential temperature and t_inc to
!       theta_star
!---------------------------------------------------------------------

! Set T to Theta and t_inc to theta_star
DO k=1, model_levels
  DO j=1, rows
    DO i=1, row_length
      theta(i,j,k)      = t(i,j,k)                                      &
                        / exner_theta_levels(i,j,k)

      ! Note: This replicates old code behaviour. This has been identified and
      !       points to a bug that may have always been present.
      !       This bug means that forcings on the 1st timestep are a
      !       factor of 86400/timestep larger that subsequent timesteps.
      !       This is to be fixed at UM8.0, theta_star may be affected in a
      !       similar manner by this bug.

      theta_star(i,j,k) = t_inc(i,j,nfor,k)                             &
                        / exner_theta_levels(i,j,k)

    END DO
  END DO
END DO

!---------------------------------------------------------------------
!     Zero diagnostics for observational forcing
!---------------------------------------------------------------------

IF (obs .AND. prindump_obs) THEN
  DO i=1, row_length
    DO m=1, rows
      DO k=1, model_levels
        DO j=1, 36
          dap1(i,m,j,k) = 0.0
          dap2(i,m,j,k) = 0.0
          DO  l = 1, nfor-1
            dap3(i,m,j,k,l) = 0.0
          END DO
        END DO
      END DO
      DO j=1, 44
        dab1(i,m,j) = 0.0
        dab2(i,m,j) = 0.0
        DO k=1, nfor-1
          dab3(i,m,j,k) = 0.0
        END DO
      END DO
    END DO                   ! i
  END DO                    ! m
END IF                     ! obs


!---------------------------------------------------------------------
!     Radiation: read spectral files
!---------------------------------------------------------------------

ErrorStatus = 0
! DEPENDS ON: initphys
CALL initphys( ErrorStatus, Cmessage )
IF (ErrorStatus  /=  0) THEN
  CALL ereport(routinename, errorstatus, cmessage)
END IF

!---------------------------------------------------------------------
!     Write restart dump information to tape
!---------------------------------------------------------------------

IF (tapeout) THEN
  nresdump = INT( ndayin / resdump_days)
  IF (MOD(ndayin, resdump_days)  /=  0) THEN
    nresdump = nresdump + 1
  END IF
  WRITE (55) exname_out,runno_out,nresdump
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE run_init

