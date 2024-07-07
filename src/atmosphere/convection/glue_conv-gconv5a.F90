! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE glue_conv_5a_mod

USE ereport_mod, ONLY: ereport
USE umPrintMgr

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GLUE_CONV_5A_MOD'

CONTAINS
!
!  Gather-scatter routine for deep and shallow convection points
!
SUBROUTINE glue_conv_5a(np_field,npnts,nlev,nbl,call_number,seg_num  &
,                 th,q,qcl,qcf, qrain, qgraup, qcf2               &
,                 cf_liquid,cf_frozen,bulk_cf,pstar               &
,                 bland,u,v,w,tracer                              &

                  ! Out
,                 dthbydt, dqbydt, dqclbydt, dqcfbydt             &
,                 dcflbydt, dcffbydt, dbcfbydt, dubydt, dvbydt    &
,                 rain, snow, rain_3d, snow_3d                    &
,                 cca0_dp, cca0_md, cca0_sh                       &
,              cca0, iccb0, icct0, cclwp0, ccw0, lcbase0, cca0_2d &
, lctop, lcca, cca,  iccb,  icct,  cclwp,  ccw,  lcbase,  cca_2d  &
,                 freeze_lev,deep_cfl_limited,mid_cfl_limited     &
,                 l_mid_all,kterm_deep,kterm_shall                &
,                 precip_deep,precip_shall,precip_mid,precip_cong &
,                 wstar_dn,wstar_up,mb1,mb2,kterm_congest         &
,                 n_cumulus,uw0,vw0,w_max,zlcl,zlcl_uv,ztop_uv    &
,                 entrain_coef,deep_flag,past_precip,past_conv_ht &
,                 cape_out                                        &
,                 n_dp,n_cg, n_sh, n_md                           &
,                 r_rho,r_theta,rho,rho_theta,rho_dry             &
,                 rho_dry_theta, delta_smag                       &
,                 exner_layer_boundaries                          &
,                 exner_layer_centres                             &
,                 p_layer_boundaries                              &
,                 p_layer_centres                                 &
,                 z_theta, z_rho                                  &
,                 timestep,t1_sd,q1_sd                            &
,                 ntml,ntpar,conv_type,l_shallow_bl               &
,                 l_pc2_diag_sh_pts, l_congestus, l_mid           &
,                 cumulus_bl,wstar,wthvs,delthvu_bl,ql_ad         &
,                 qsat_lcl, ftl, fqt                              &
,                 l_tracer, ntra, trlev, n_cca_lev                &
,                 l_mcr_qrain, l_mcr_qgraup                       &
,                 l_mcr_qcf2, l_calc_dxek, l_q_interact           &
,                 up_flux_half, up_flux, dwn_flux                 &
,                 entrain_up,detrain_up,entrain_dwn,detrain_dwn   &
,                 uw_deep,vw_deep,uw_shall,vw_shall,uw_mid,vw_mid &
,         wqt_flux_sh,wthetal_flux_sh,wthetav_flux_sh,wql_flux_sh &
,          mf_deep,mf_congest,mf_shall,mf_midlev                  &
,          dt_deep,dt_congest,dt_shall,dt_midlev                  &
,          dq_deep,dq_congest,dq_shall,dq_midlev                  &
,          du_deep,du_congest,du_shall,du_midlev                  &
,          dv_deep,dv_congest,dv_shall,dv_midlev                  &
,       ind_cape_reduced, cape_ts_used, ind_deep, ind_shall, w2p  &
,       dt_dd, dq_dd, du_dd, dv_dd, area_ud, area_dd  )

! Purpose:
!  Gather-scatter routine for deep and shallow convection points.
!  Interface to deep, shallow and mid level convection

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3 programming standards v8.3

USE cv_run_mod,  ONLY:                                            &
    l_mom,         adapt,           termconv,                     &
    cca2d_sh_opt,  cca2d_md_opt,    cca2d_dp_opt,                 &
    cca_sh_knob,   cca_md_knob,     cca_dp_knob,                  &
    ccw_sh_knob,   ccw_md_knob,     ccw_dp_knob,                  &
    iconv_shallow, iconv_congestus, iconv_mid,    iconv_deep,     &
    l_anvil, l_ccrad, l_3d_cca, l_conv_hist, l_pc2_diag_sh,       &
    l_new_dd

USE cv_stash_flg_mod, ONLY:                                       &
    flg_up_flx, flg_up_flx_half, flg_dwn_flx,                     &
    flg_entr_up, flg_entr_dwn, flg_detr_up, flg_detr_dwn,         &
    flg_uw_dp, flg_vw_dp, flg_uw_shall, flg_vw_shall,             &
    flg_uw_mid, flg_vw_mid,                                       &
    flg_wqt_flux, flg_wthetal_flux,                               &
    flg_wthetav_flux, flg_wql_flux,                               &
    flg_mf_deep, flg_mf_congest, flg_mf_shall, flg_mf_midlev,     &
    flg_dt_deep, flg_dt_congest, flg_dt_shall, flg_dt_midlev,     &
    flg_dq_deep, flg_dq_congest, flg_dq_shall, flg_dq_midlev,     &
    flg_du_deep, flg_du_congest, flg_du_shall, flg_du_midlev,     &
    flg_dv_deep, flg_dv_congest, flg_dv_shall, flg_dv_midlev,     &
    flg_w_eqn,   flg_area_ud,    flg_area_dd

USE cv_hist_constants_mod, ONLY:                                  &
    decay_period

USE water_constants_mod, ONLY: tm

USE planet_constants_mod, ONLY: g, cp

USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
USE ereport_mod, ONLY: ereport
USE rad_input_mod, ONLY: l_cca_dp_prog, l_cca_md_prog, l_cca_sh_prog

USE model_domain_mod,      ONLY: model_type, mt_single_column
USE atm_fields_bounds_mod, ONLY: ScmRowLen, ScmRow
USE s_scmop_mod,           ONLY: default_streams,                       &
                                 t_inst,d_wet,d_point
USE scmoutput_mod,         ONLY: scmoutput
USE errormessagelength_mod, ONLY: errormessagelength

USE congest_conv_mod, ONLY: congest_conv
USE deep_conv_5a_mod, ONLY: deep_conv_5a
USE deep_turb_conv_mod, ONLY: deep_turb_conv
USE mid_conv_5a_mod, ONLY: mid_conv_5a
USE shallow_conv_5a_mod, ONLY: shallow_conv_5a
USE shallow_turb_conv_mod, ONLY: shallow_turb_conv

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    l_new_qsat_conv !Currently defaults to FALSE

USE gen_phys_inputs_mod, ONLY: l_mr_physics

IMPLICIT NONE

! Model Constants required


!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent IN:


INTEGER, INTENT(IN) :: &
  np_field             & ! No of points in a full field (note all
                         ! multi-dimensional fields passed to routine MUST be
                         ! dimensioned with np_field)
 ,npnts                & ! No. of points in segment
 ,nlev                 & ! No. of model layers
 ,nbl                  & ! No. of boundary layer levels
 ,call_number          & ! convection call number with physics timestep
 ,seg_num                ! segment number

! NOTE - loops over points should be over npnts or less NEVER over np_field.

INTEGER, INTENT(IN) :: ntml(npnts) ! Top level of surface mixed
                                   ! layer defined relative to
                                   ! theta,q grid

INTEGER, INTENT(IN) :: ntpar(npnts) ! Top level of initial parcel
                                    ! ascent in BL scheme defined
                                    ! relative to theta,q grid

INTEGER, INTENT(IN) :: n_cca_lev! No. of convective cloud
                                ! amount levels (1 for 2D,
                                ! nlevs for 3D)

INTEGER, INTENT(IN) :: ntra     ! No. of tracer fields

INTEGER, INTENT(IN) :: trlev    ! No. of model levels on which
                                ! tracers are included

INTEGER, INTENT(IN) :: &
  n_dp                 & ! Number of deep points
 ,n_cg                 & ! Number of congestus points
 ,n_sh                 & ! Number of shallow points
 ,n_md                   ! Number of mid points


REAL, INTENT(IN) :: &
  wstar(npnts)      & ! Convective velocity scale (m/s)
 ,wthvs(npnts)      & ! Surface flux of THV  (Pa m/s2)
 ,delthvu_bl(npnts) & ! Integral of undilute parcel buoyancy over convective
                      ! cloud layer (Kelvin m)
 ,delta_smag(npnts) & ! grid size used in smagorinsky length scale
 ,ql_ad(npnts)      & ! adiabatic liquid water content at inversion (kg/kg)
 ,qsat_lcl(npnts)   & ! qsat at cloud base (kg/kg)
 ,ftl(npnts)        & ! Surface sensible heat flux from BL (W/m2)
                      !  i.e. cp*rho*w'tl'
 ,fqt(npnts)          ! Total water flux from surface (kg/m2/s) i.e. rho*w'qT'


REAL, INTENT(IN) ::               &
  r_rho(np_field,nlev)            & ! radius rho levels(m)
 ,r_theta(np_field,0:nlev)        & ! theta levels (m)
 ,rho(np_field,nlev)              & ! density kg/m3
 ,rho_theta(np_field,nlev)        & ! th lev kg/m3
 ,rho_dry(np_field,nlev)          & ! dry denisty on rho levels (kg/m3)
 ,rho_dry_theta(np_field,nlev)      ! dry denisty on theta levels (kg/m3)

REAL, INTENT(IN) ::                       &
  exner_layer_centres(np_field,0:nlev)    & ! Exner
 ,exner_layer_boundaries(np_field,0:nlev) & ! Exner at half level above
                                            ! exner_layer_centres
 ,p_layer_centres(np_field,0:nlev)        & ! Pressure (Pa)
 ,p_layer_boundaries(np_field,0:nlev)       ! Pressure at half level above
                                            ! p_layer_centres (Pa)

REAL, INTENT(IN) :: pstar(npnts) ! Surface pressure (Pa)

REAL, INTENT(IN) ::      &
  z_theta(np_field,nlev) & ! height of theta levels above surface (m)
 ,z_rho(np_field,nlev)     ! height of rho levels above surface (m)


REAL, INTENT(IN) ::   &
  t1_sd(npnts)        & ! Standard deviation of turbulent fluctuations of
                        ! layer 1 temp. (K)
 ,q1_sd(npnts)          ! Standard deviation of turbulent fluctuations of
                        ! layer 1 q (kg/kg)

REAL, INTENT(IN) ::  &
  th(np_field,nlev)  & ! Model potential temperature (K)
 ,q(np_field,nlev)     ! Model water vapour (kg/kg)

REAL, INTENT(IN) :: timestep    ! Model timestep (s)

REAL, INTENT(IN) ::  &
  uw0(npnts)         & ! U-comp of surface stress(N/m2)
 ,vw0(npnts)         & ! V-comp of surface stress(N/m2)
 ,w_max(npnts)       & ! max w in column
 ,zlcl(npnts)        & ! Lifting condensation level accurate
                       ! height (m) not a model level.
 ,zlcl_uv(npnts)     & ! Lifting condensation level
                       ! defined for the uv grid (m)
 ,ztop_uv(npnts)       ! Top of cloud layer defined for the uv grid (m)

REAL, INTENT(IN) ::  &
  u(np_field,nlev)   & ! Model U field (m/s)
 ,v(np_field,nlev)   & ! Model V field (m/s)
 ,w(np_field,nlev)     ! Model W field (m/s)


REAL, INTENT(IN) ::    &
  entrain_coef(npnts)     ! entrainment coefficients

LOGICAL, INTENT(IN) :: l_tracer ! Switch for inclusion of tracers

LOGICAL, INTENT(IN) ::     &
  l_shallow_bl(npnts)      & ! Shallow cumulus indicator
 ,l_pc2_diag_sh_pts(npnts) & ! Carry diagnostic shallow convective
                             ! information for PC2
 ,l_congestus(npnts)       & ! congestus cumulus
 ,l_mid(npnts)             & ! possible mid-level cnv
 ,cumulus_bl(npnts)          ! Cumulus indicator

LOGICAL, INTENT(IN) :: &
  l_mcr_qrain          & ! true - prognostic rain
 ,l_mcr_qgraup         & ! true - prognostic graupel
 ,l_mcr_qcf2           & ! true - prognostic qcf2
 ,l_calc_dxek          & ! Switch for calculation of condensate increment
 ,l_q_interact           ! Switch allows overwriting parcel variables when
                         ! calculating condensate incr.

LOGICAL, INTENT(IN) :: bland(npnts) ! Land/sea mask


! Arguments with intent INOUT:

! History prognostics only in use if l_conv_hist = .true.

REAL, INTENT(INOUT) :: &
 deep_flag(npnts)      & ! 0.0-1.0, 1. indicates deep last time step
,past_precip(npnts)    & ! convective precip rate or a decayed value
                         ! (kg/m2/s)
,past_conv_ht(npnts)     ! convective height (m)

! NOTE - All moist variables passed down to this routine are :
! specific humidities if l_mr_physics = .false.
! mixing ratios       if l_mr_physics = .true.

REAL, INTENT(INOUT) ::      &
  qcl(np_field,nlev)        & ! Liq condensate  (kg/kg)
 ,qcf(np_field,nlev)        & ! Ice condensate  (kg/kg)
 ,qrain(np_field,nlev)      & ! rain (kg/kg)
 ,qgraup(np_field,nlev)     & ! graupel (kg/kg)
 ,qcf2(np_field,nlev)       & ! 2nd ice type (kg/kg)
 ,cf_liquid(np_field,nlev)  & ! Liq water cloud volume (fraction)
 ,cf_frozen(np_field,nlev)  & ! Frozen water cloud volume (fraction?)
 ,bulk_cf(np_field,nlev)      ! Bulk total cloud volume ( )

REAL, INTENT(INOUT) :: tracer(np_field,trlev,ntra) !Model tracer
                                                   ! fields (kg/kg)

REAL, INTENT(INOUT) ::     &
  w2p(np_field,nlev)         ! (Parcel vertical velocity)^2 [(m/s)^2]

! Arguments with intent OUT:

REAL, INTENT(OUT) ::       &
  dthbydt(np_field,nlev)   & ! Increments to potential temp. due to
                             ! convection (K/s)
 ,dqbydt(np_field,nlev)    & ! Increments to q due to convection (kg/kg/s)

 ,dqclbydt(np_field,nlev)  & ! Increments to liq condensate due to convection
                             ! (kg/kg/s)
 ,dqcfbydt(np_field,nlev)  & ! Increments to ice condensate due to convection
                             ! (kg/kg/s)
 ,dcflbydt(np_field,nlev)  & ! Increments to liq cloud volume due to convection
                             ! (/s)
 ,dcffbydt(np_field,nlev)  & ! Increments to ice cloud volume due to convection
                             ! (/s)
 ,dbcfbydt(np_field,nlev)  & ! Increments to total cld volume due to convection
                             ! (/s)
 ,dubydt(np_field,nlev+1)  & ! Increments to U due to CMT (m/s2)

 ,dvbydt(np_field,nlev+1)    ! Increments to V due to CMT (m/s2)


REAL, INTENT(OUT) ::       &
  rain(npnts)              & ! Surface convective rainfall (kg/m2/s)
 ,snow(npnts)              & ! Surface convective snowfall (kg/m2/s)
 ,rain_3d(np_field,nlev)   & ! convective rainfall flux (kg/m2/s)
 ,snow_3d(np_field,nlev)     ! convective snowfall flux (kg/m2/s)

! Section 5 Convective Cloud properties
REAL, INTENT(OUT) ::       &
  cca(np_field,n_cca_lev)  &! Cnv. cld amount (0-1)
, ccw(np_field,nlev)       &! Cnv. in-cld water (kg/kg)
, cclwp(npnts)             &! Cond. water path (kg/m^2)
, lcca(npnts)               ! Cnv. cld amount at base of lowest
                            ! cloud layer (no anvil). (0-1)

INTEGER, INTENT(OUT) ::    &
  iccb(npnts)              &! Cnv. cld base level (Highest Cld layer)
, icct(npnts)              &! Cnv. cld top  level (Highest Cld layer)
, lcbase(npnts)            &! Cnv. cld base level (Lowest  Cld layer)
, lctop(npnts)              ! Cnv. cld top  level (Lowest  Cld layer)


! Section 0 Convective Cloud properties for Radiative impacts
REAL ::       &
  cca0(np_field,n_cca_lev) &! Cnv. cld amount (0-1)
, cca0_sh(np_field,n_cca_lev) &! Shallow cnv. cld amount (0-1)
, cca0_md(np_field,n_cca_lev) &! Mid-level cnv. cld amount (0-1)
, cca0_dp(np_field,n_cca_lev) &! Deep cnv. cld amount (0-1)
, ccw0(np_field,nlev)      &! Cnv. in-cld water (kg/kg)
, cclwp0(npnts)             ! Cond. water path (kg/m^2)

INTEGER, INTENT(OUT) ::    &
  iccb0(npnts)             &! Cnv. cld base level (Highest Cld layer)
, icct0(npnts)             &! Cnv. cld top  level (Highest Cld layer)
, lcbase0(npnts)            ! Cnv. cld base level (Lowest  Cld layer)

INTEGER, INTENT(OUT) :: freeze_lev(npnts) ! index for freezing lev

INTEGER, INTENT(OUT) :: kterm_deep(npnts) ! index deep conv
INTEGER, INTENT(OUT) :: kterm_shall(npnts) ! Level shallow terminates

LOGICAL, INTENT(OUT) :: l_mid_all(npnts)  ! on exit true if mid level
                                          ! convection has triggered
REAL, INTENT(OUT) ::      &
  deep_cfl_limited(npnts) & ! indicator for cfl limited deep convection
 ,mid_cfl_limited(npnts)    ! indicator for cfl limited mid-level convection

REAL, INTENT(OUT) ::   &
  precip_deep(npnts)   & ! deep precip (kg/m2/s)
 ,precip_shall(npnts)  & ! shallow precip(kg/m2/s)
 ,precip_mid(npnts)    & ! mid precip (kg/m2/s)
 ,precip_cong(npnts)     ! congest precip (kg/m2/s)

REAL, INTENT(OUT) ::  &
  wstar_dn(npnts)     & ! subcloud layer convective velocity scale(m/s)
, wstar_up(npnts)     & ! cumulus layer convective velocity scale (m/s)
, mb1(npnts)          & ! cloud base mass flux from wstar_dn (m/s)
, mb2(npnts)            ! cloud base mass flux for cloud layer (m/s)

INTEGER, INTENT(OUT) :: kterm_congest(npnts)
                        ! termination level for congestus

! Meaning of this diagnostic depends on scheme
!  Plume model - updraught mass flux  (Pa/s)
!  turbulence model - mass flux (not exactly updraught) (m/s)


REAL, INTENT(OUT) ::          &
  up_flux_half(np_field,nlev) & ! mass flux on rho
                                ! dummy variable not used in turbulence
 ,up_flux(np_field,nlev)        ! mass flux

! Diagnostics with no meaning for turbulence based schemes

REAL, INTENT(OUT) ::          &
  dwn_flux(np_field,nlev)     & ! Downdraught mass flux (Pa/s)
 ,entrain_up(np_field,nlev)   & ! Fractional entrainment rate into updraughts
                                ! (Pa/s)
 ,detrain_up(np_field,nlev)   & ! Fractional detrainment rate into updraughts
                                ! (Pa/s)
 ,entrain_dwn(np_field,nlev)  & ! Fractional entrainment rate into
                                ! downdraughts (Pa/s)
 ,detrain_dwn(np_field,nlev)    ! Fractional detrainment rate into
                                ! downdraughts (Pa/s)

! Diagnostics relating to momentum fluxes

REAL, INTENT(OUT) ::     &
 uw_deep(np_field,nlev)  & ! X-comp. of stress from deep convection (kg/m/s2)
,vw_deep(np_field,nlev)  & ! Y-comp. of stress from deep convection (kg/m/s2)
,uw_shall(np_field,nlev) & ! X-comp. of stress from shallow convection (kg/m/s2)
,vw_shall(np_field,nlev) & ! Y-comp. of stress from shallow convection (kg/m/s2)
,uw_mid(np_field,nlev)   & ! U comp of stress from mid convection (kg/m/s2)
,vw_mid(np_field,nlev)     ! V comp of stress from mid convection (kg/m/s2)

REAL, INTENT(OUT) :: cape_out(npnts) ! Saved convective available
                                     ! potential energy for diagnostic
                                     ! output (Jkg-1)

! Fluxes from turbulence based convection schemes

REAL, INTENT(OUT) ::             &
  wqt_flux_sh(np_field,nlev)     & ! w'qt' flux (m/s kg/kg)
, wthetal_flux_sh(np_field,nlev) & ! w'thetal' flux  (m/s K)
, wthetav_flux_sh(np_field,nlev) & ! w'thetav' flux  (m/s K)
, wql_flux_sh(np_field,nlev)     & ! w'ql' flux  (m/s kg/kg)
, mf_deep(np_field,nlev)         & ! mass flux deep
, mf_congest(np_field,nlev)      & ! mass flux congestus
, mf_shall(np_field,nlev)        & ! mass flux shallow
, mf_midlev(np_field,nlev)       & ! mass flux mid-lev
, dt_deep(np_field,nlev)         & ! dt increment deep   (K/s)
, dt_congest(np_field,nlev)      & ! dt increment congestus (K/s)
, dt_shall(np_field,nlev)        & ! dt increment shallow (K/s)
, dt_midlev(np_field,nlev)       & ! dt increment mid-level (K/s)
, dq_deep(np_field,nlev)         & ! dq increment deep (kg/kg/s)
, dq_congest(np_field,nlev)      & ! dq increment congestus (kg/kg/s)
, dq_shall(np_field,nlev)        & ! dq increment shallow (kg/kg/s)
, dq_midlev(np_field,nlev)       & ! dq increment mid-level (kg/kg/s)
, du_deep(np_field,nlev+1)       & ! du increment deep (m/s)
, du_congest(np_field,nlev+1)    & ! du increment congestus (m/s)
, du_shall(np_field,nlev+1)      & ! du increment shallow (m/s)
, du_midlev(np_field,nlev+1)     & ! du increment mid-level (m/s)
, dv_deep(np_field,nlev+1)       & ! dv increment deep (m/s)
, dv_congest(np_field,nlev+1)    & ! dv increment congestus (m/s)
, dv_shall(np_field,nlev+1)      & ! dv increment shallow (m/s)
, dv_midlev(np_field,nlev+1)       ! dv increment mid-level (m/s)

REAL, INTENT(OUT) ::         &
 ind_cape_reduced(npnts)     & ! indicates cape timescale reduced
,cape_ts_used(npnts)         & ! cape timescale used for deep (s)
,ind_deep(npnts)             & ! indicator of real deep
,ind_shall(npnts)              ! indicator of real shallow

REAL, INTENT(OUT) ::         &
 dt_dd(np_field,nlev)        & ! dT/dt from DD & evap below cloud base (K/s)
,dq_dd(np_field,nlev)        & ! dq/dt from DD & evap below cloud base (kg/kg/s)
,du_dd(np_field,nlev)        & ! du/dt from DD  (m/s/s)
,dv_dd(np_field,nlev)        & ! dv/dt from DD  (m/s/s)
,area_ud(np_field,nlev)      & ! updraught fractional area
,area_dd(np_field,nlev)        ! downdraught fractional area

!-----------------------------------------------------------------------
! Redundant arguments
!-----------------------------------------------------------------------

INTEGER, INTENT(IN) :: n_cumulus          ! (NOT USED)
INTEGER, INTENT(IN) :: conv_type(npnts)   ! (NOT USED)


! local required in move mixing to glue

REAL :: dtrabydt(npnts,nlev,ntra) ! Increment to tracer due to
                                  ! convection (kg/kg/s)

CHARACTER (LEN=*), PARAMETER ::  RoutineName = 'GLUE_CONV_5A'
!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to DEEP convection scheme.
! Arrays are identified by underscore DP (_dp) and are of length
! n_dp where n_dp is the number of points diagnosed as deep in the
! boundary layer diagnosis routine. For full desriptions of variables
! see above.
!-----------------------------------------------------------------------
INTEGER :: error_point      ! location of problem deep point

INTEGER :: dpi(n_dp)        ! index for deep points in full grid

INTEGER ::                                                        &
  ntml_dp(n_dp)                                                   &
, ntpar_dp(n_dp)

REAL ::                 &
  pstar_dp(n_dp)        &
, recip_pstar_dp(n_dp)  &
, t1_sd_dp(n_dp)        &
, q1_sd_dp(n_dp)        &
, uw0_dp(n_dp)          &
, vw0_dp(n_dp)          &
, fqt_dp(n_dp)          &!  (Deep_turb_conv)
, ftl_dp(n_dp)          &!  (Deep_turb_conv)
, zlcl_dp(n_dp)         &!  (Deep_turb_conv)
, zlcl_uv_dp(n_dp)      &
, wth0_dp(n_dp)         &!  (Deep_turb_conv)
, wq0_dp(n_dp)          &!  (Deep_turb_conv)
, wstar_dp(n_dp)        &
, wthvs_dp(n_dp)        &!  (Deep_turb_conv)
, delthvu_dp(n_dp)      &
, entrain_coef_dp(n_dp) &
, qsat_lcl_dp(n_dp)     &
, mb_dp(n_dp)

REAL ::                                                            &
  p_layer_centres_dp(n_dp,0:nlev)                                  &
, p_layer_boundaries_dp(n_dp,0:nlev)                               &
, exner_layer_centres_dp(n_dp,0:nlev)                              &
, exner_layer_boundaries_dp(n_dp,0:nlev)                           &
, r2rho_dp(n_dp,nlev)                                              &
, r2rho_th_dp(n_dp,nlev)                                           &
, rho_dp(n_dp,nlev)                                                &
, rho_theta_dp(n_dp,nlev)                                          &
, dr_across_th_dp(n_dp,nlev)                                       &
, dr_across_rh_dp(n_dp,nlev)                                       &
, z_theta_dp(n_dp,nlev)                                            &
, z_rho_dp(n_dp,nlev)                                              &
, r_rho_dp(n_dp,nlev)                                              &
, r_theta_dp(n_dp,0:nlev)

LOGICAL :: bland_dp(n_dp)

REAL ::                &
  u_dp(n_dp,nlev)      &
, v_dp(n_dp,nlev)      &
, w_dp(n_dp,nlev)      &
, th_dp(n_dp,nlev)     &
, q_dp(n_dp,nlev)      &
, qse_dp(n_dp,nlev)

REAL :: tracer_dp(n_dp,trlev,ntra)
! increments

REAL ::                                                           &
  dthbydt_dp(n_dp,nlev)                                           &
, dqbydt_dp(n_dp,nlev)                                            &
, dubydt_dp(n_dp,nlev)                                            &
, dvbydt_dp(n_dp,nlev)

! output variables
REAL ::                                                           &
  rain_dp(n_dp)                                                   &
, snow_dp(n_dp)                                                   &
, rain_3d_dp(n_dp,nlev)                                           &
, snow_3d_dp(n_dp,nlev)                                           &
, tcw_dp(n_dp)                                                    &
, cclwp_dp(n_dp)                                                  &
, lcca_dp(n_dp)                                                   &
, cape_out_dp(n_dp)

INTEGER ::                                                        &
  iccb_dp(n_dp)                                                   &
, icct_dp(n_dp)                                                   &
, lcbase_dp(n_dp)                                                 &
, lctop_dp(n_dp)                                                  &
, freeze_lev_dp(n_dp)

REAL ::                                                           &
  ccw_dp(n_dp,nlev)                                               &
, up_flux_dp(n_dp,nlev)                                           &
, up_flux_half_dp(n_dp,nlev)                                      &
, dwn_flux_dp(n_dp,nlev)                                          &
, entrain_up_dp(n_dp,nlev)                                        &
, detrain_up_dp(n_dp,nlev)                                        &
, entrain_dwn_dp(n_dp,nlev)                                       &
, detrain_dwn_dp(n_dp,nlev)                                       &
, uw_deep_dp(n_dp,nlev)                                           &
, vw_deep_dp(n_dp,nlev)                                           &
, w_max_dp(n_dp)

!PC2
REAL ::                                                           &
  qcl_dp(n_dp,nlev)                                               &
, qcf_dp(n_dp,nlev)                                               &
, cf_liquid_dp(n_dp,nlev)                                         &
, cf_frozen_dp(n_dp,nlev)                                         &
, bulk_cf_dp(n_dp,nlev)                                           &
, dqclbydt_dp(n_dp,nlev)                                          &
, dqcfbydt_dp(n_dp,nlev)                                          &
, dcflbydt_dp(n_dp,nlev)                                          &
, dcffbydt_dp(n_dp,nlev)                                          &
, dbcfbydt_dp(n_dp,nlev)

REAL :: dtrabydt_dp(n_dp,nlev,ntra) ! Increment to tracer due to
                                    ! convection (kg/kg/s)
INTEGER :: kterm_dp(n_dp)    ! required by mid level scheme

REAL ::                     &
  cca_2d_dp(n_dp)           & ! required by mid level scheme
, cca_dp(n_dp,n_cca_lev)    & ! 3d CCA
, ind_cape_reduced_dp(n_dp) & ! indicates reduced cape timescale
, cape_ts_used_dp(n_dp)     & ! cape timescale deep
, cfl_limited_dp(n_dp)      & ! CFL limited conv
, ind_deep_dp(n_dp)           ! indicator of real deep

! temporary arrays contianing copies of parts of full arrays for deep,
! shallow, congestus or mid-level schemes

REAL, ALLOCATABLE ::        &
  dt_dd_temp(:,:)           & ! dT/dt from DD & evap below cloud from scheme
 ,dq_dd_temp(:,:)           & ! dq/dt from DD & evap below cloud from scheme
 ,du_dd_temp(:,:)           & ! du/dt from DD 
 ,dv_dd_temp(:,:)           & ! dv/dt from DD
 ,area_ud_temp(:,:)         & ! updraught fractional area
 ,area_dd_temp(:,:)           ! downdraught fractional area

REAL, ALLOCATABLE ::        &
  rho_dry_temp(:,:)         & ! dry density on rho levels (kg/m3)
 ,rho_dry_theta_temp(:,:)     ! dry density on theta levels (kg/m3)


INTEGER :: ishall_precip  ! flag for precip in shallow turbulence
                          ! scheme
                          !  0 - no precip
                          !  1 - precip

! Compressed SCM diagnostics for adaptive modset
REAL :: rbuoy_p_out_dp (ScmRowLen*ScmRow,nlev)
REAL :: the_out_dp     (ScmRowLen*ScmRow,nlev)
REAL :: thp_out_dp     (ScmRowLen*ScmRow,nlev)
REAL :: qe_out_dp      (ScmRowLen*ScmRow,nlev)
REAL :: qp_out_dp      (ScmRowLen*ScmRow,nlev)

!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to congestus convection scheme.
! Arrays are identified by underscore cg (_cg) and are of length
! n_cg where n_cg is the number of points diagnosed as congestus in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

INTEGER :: cgi(n_cg)     ! index for congestus points in full grid

INTEGER ::                                                        &
  ntml_cg(n_cg)                                                   &
, ntpar_cg(n_cg)

REAL ::                                                           &
  pstar_cg(n_cg)                                                  &
, recip_pstar_cg(n_cg)                                            &
, delthvu_cg(n_cg)                                                &
, t1_sd_cg(n_cg)                                                  &
, q1_sd_cg(n_cg)                                                  &
, uw0_cg(n_cg)                                                    &
, vw0_cg(n_cg)                                                    &
, wstar_cg(n_cg)                                                  &
, wthvs_cg(n_cg)                                                  &
, zlcl_uv_cg(n_cg)                                                &
, ztop_uv_cg(n_cg)                                                &
, entrain_coef_cg(n_cg)

REAL ::                                                           &
  p_layer_centres_cg(n_cg,0:nlev)                                 &
, p_layer_boundaries_cg(n_cg,0:nlev)                              &
, exner_layer_centres_cg(n_cg,0:nlev)                             &
, exner_layer_boundaries_cg(n_cg,0:nlev)                          &
, z_theta_cg(n_cg,nlev)                                           &
, z_rho_cg(n_cg,nlev)                                             &
, u_cg(n_cg,nlev)                                                 &
, v_cg(n_cg,nlev)                                                 &
, th_cg(n_cg,nlev)                                                &
, q_cg(n_cg,nlev)                                                 &
, ccw_cg(n_cg,nlev)                                               &
, qse_cg(n_cg,nlev)                                               &
, rho_cg(n_cg,nlev)                                               &
, r_rho_cg(n_cg,nlev)                                             &
, r_theta_cg(n_cg,0:nlev)

REAL ::                      &
  r2rho_th_cg(n_cg,nlev)     & ! radius**2 density theta lev (kg/m)
, r2rho_cg(n_cg,nlev)        & ! radius**2 density rho lev (kg/m)
, dr_across_th_cg(n_cg,nlev) & ! thickness of theta levels (m)
, dr_across_rh_cg(n_cg,nlev) & ! thickness of rho levels (m)
, rho_th_cg(n_cg,nlev)         ! rho on theta levels (kg/m3)

LOGICAL :: bland_cg(n_cg)

! tracers
REAL :: tracer_cg(n_cg,trlev,ntra)                                &
, dtrabydt_cg(n_cg,nlev,ntra)

REAL ::                                                           &
  dthbydt_cg(n_cg,nlev)                                           &
, dqbydt_cg(n_cg,nlev)                                            &
, dubydt_cg(n_cg,nlev)                                            &
, dvbydt_cg(n_cg,nlev)

REAL ::                                                           &
  rain_cg(n_cg)                                                   &
, snow_cg(n_cg)                                                   &
, rain_3d_cg(n_cg,nlev)                                           &
, snow_3d_cg(n_cg,nlev)                                           &
, tcw_cg(n_cg)                                                    &
, cclwp_cg(n_cg)                                                  &
, lcca_cg(n_cg)                                                   &
, cape_out_cg(n_cg)

REAL ::                  &
 cca_2d_cg(n_cg)         & ! required by mid level scheme
,cca_cg(n_cg,n_cca_lev)    ! 3d CCA

INTEGER ::                                                        &
  iccb_cg(n_cg)                                                   &
, icct_cg(n_cg)                                                   &
, lcbase_cg(n_cg)                                                 &
, lctop_cg(n_cg)                                                  &
, freeze_lev_cg(n_cg)                                             &
, kterm_cg(n_cg)    ! required by mid level scheme

! diagnostics
REAL ::                                                           &
  up_flux_cg(n_cg,nlev)                                           &
, up_flux_half_cg(n_cg,nlev)                                      &
, dwn_flux_cg(n_cg,nlev)                                          &
, entrain_up_cg(n_cg,nlev)                                        &
, detrain_up_cg(n_cg,nlev)                                        &
, entrain_dwn_cg(n_cg,nlev)                                       &
, detrain_dwn_cg(n_cg,nlev)                                       &
, uw_shall_cg(n_cg,nlev)                                          &
, vw_shall_cg(n_cg,nlev)

!PC2
REAL ::                                                           &
  qcl_cg(n_cg,nlev)                                               &
, qcf_cg(n_cg,nlev)                                               &
, cf_liquid_cg(n_cg,nlev)                                         &
, cf_frozen_cg(n_cg,nlev)                                         &
, bulk_cf_cg(n_cg,nlev)                                           &
, dqclbydt_cg(n_cg,nlev)                                          &
, dqcfbydt_cg(n_cg,nlev)                                          &
, dcflbydt_cg(n_cg,nlev)                                          &
, dcffbydt_cg(n_cg,nlev)                                          &
, dbcfbydt_cg(n_cg,nlev)


! Compressed SCM diagnostics for adaptive modset
REAL :: rbuoy_p_out_cg (ScmRowLen*ScmRow,nlev)
REAL :: the_out_cg     (ScmRowLen*ScmRow,nlev)
REAL :: thp_out_cg     (ScmRowLen*ScmRow,nlev)
REAL :: qe_out_cg      (ScmRowLen*ScmRow,nlev)
REAL :: qp_out_cg      (ScmRowLen*ScmRow,nlev)


!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to SHALLOW convection scheme.
! Arrays are identified by underscore SH (_sh) and are of length
! n_sh where n_sh is the number of points diagnosed as shallow in the
! boundary layer diagnosis routine.For full desriptions of variables
! see above.
!-----------------------------------------------------------------------

INTEGER :: shi(n_sh)      ! index for shallow points in full grid

! compressed inputs

INTEGER ::          &
  ntml_sh(n_sh)     &
, ntpar_sh(n_sh)    &
, kterm_sh(n_sh)      ! May be needed in adaptive option 4 or 6

REAL ::                                                           &
  pstar_sh(n_sh)                                                  &
, recip_pstar_sh(n_sh)                                            &
, delthvu_sh(n_sh)                                                &
, ql_ad_sh(n_sh)                                                  &
, t1_sd_sh(n_sh)                                                  &
, q1_sd_sh(n_sh)                                                  &
, uw0_sh(n_sh)                                                    &
, vw0_sh(n_sh)                                                    &
, wstar_sh(n_sh)                                                  &
, wthvs_sh(n_sh)                                                  &
, zlcl_sh(n_sh)                                                   &
, zlcl_uv_sh(n_sh)                                                &
, ztop_uv_sh(n_sh)                                                &
, entrain_coef_sh(n_sh)                                           &
, delta_smag_sh(n_sh)                                             &
, ind_shall_sh(n_sh)  ! indicator of real shallow

REAL ::                                  &
  p_layer_centres_sh(n_sh,0:nlev)        &
, p_layer_boundaries_sh(n_sh,0:nlev)     &
, exner_layer_centres_sh(n_sh,0:nlev)    &
, exner_layer_boundaries_sh(n_sh,0:nlev) &
, z_theta_sh(n_sh,nlev)                  &
, z_rho_sh(n_sh,nlev)                    &
, rho_sh(n_sh,nlev)                      &
, r_rho_sh(n_sh,nlev)                    &
, r_theta_sh(n_sh,0:nlev)                &
, r2rho_th_sh(n_sh,nlev)                 & ! radius**2 density theta lev (kg/m)
, r2rho_sh(n_sh,nlev)                    & ! radius**2 density rho lev (kg/m)
, dr_across_th_sh(n_sh,nlev)             & ! thickness of theta levels (m)
, dr_across_rh_sh(n_sh,nlev)             & ! thickness of rho levels (m)
, rho_th_sh(n_sh,nlev)                     ! rho on theta levels (kg/m3)

LOGICAL :: bland_sh(n_sh)

REAL ::                                                           &
  u_sh(n_sh,nlev)                                                 &
, v_sh(n_sh,nlev)                                                 &
, w_sh(n_sh,nlev)                                                 &
, th_sh(n_sh,nlev)                                                &
, q_sh(n_sh,nlev)                                                 &
, qse_sh(n_sh,nlev)

! output

REAL ::                                                           &
  dthbydt_sh(n_sh,nlev)                                           &
, dqbydt_sh(n_sh,nlev)                                            &
, dubydt_sh(n_sh,nlev)                                            &
, dvbydt_sh(n_sh,nlev)

REAL ::                 &
  rain_sh(n_sh)         &
, snow_sh(n_sh)         &
, rain_3d_sh(n_sh,nlev) &
, snow_3d_sh(n_sh,nlev) &
, tcw_sh(n_sh)          &
, cclwp_sh(n_sh)        &
, lcca_sh(n_sh)         &
, cape_out_sh(n_sh)     &
, cca_2d_sh(n_sh)       &   ! required by mid level scheme
, cca_sh(n_sh,n_cca_lev)    ! 3d CCA

REAL ::                                                           &
  wqt_sh(n_sh,nlev)                                               &
, wthetal_sh(n_sh,nlev)                                           &
, wthetav_sh(n_sh,nlev)                                           &
, wql_sh(n_sh,nlev)                                               &
, wstar_up_sh(n_sh)                                               &
, mb1_sh(n_sh)                                                    &
, mb2_sh(n_sh)

INTEGER :: iccb_sh(n_sh)                                          &
, icct_sh(n_sh)                                                   &
, lcbase_sh(n_sh)                                                 &
, lctop_sh(n_sh)                                                  &
, freeze_lev_sh(n_sh)

REAL ::                                                           &
  ccw_sh(n_sh,nlev)                                               &
, up_flux_sh(n_sh,nlev)                                           &
, up_flux_half_sh(n_sh,nlev)                                      &
, dwn_flux_sh(n_sh,nlev)                                          &
, entrain_up_sh(n_sh,nlev)                                        &
, detrain_up_sh(n_sh,nlev)                                        &
, entrain_dwn_sh(n_sh,nlev)                                       &
, detrain_dwn_sh(n_sh,nlev)                                       &
, uw_shall_sh(n_sh,nlev)                                          &
, vw_shall_sh(n_sh,nlev)

! PC2  input
REAL ::                                                           &
  qcl_sh(n_sh,nlev)                                               &
, qcf_sh(n_sh,nlev)                                               &
, cf_liquid_sh(n_sh,nlev)                                         &
, cf_frozen_sh(n_sh,nlev)                                         &
, bulk_cf_sh(n_sh,nlev)
! PC2 output
REAL ::                                                           &
  dqclbydt_sh(n_sh,nlev)                                          &
, dqcfbydt_sh(n_sh,nlev)                                          &
, dcflbydt_sh(n_sh,nlev)                                          &
, dcffbydt_sh(n_sh,nlev)                                          &
, dbcfbydt_sh(n_sh,nlev)

! Tracers in and out
REAL :: tracer_sh(n_sh,trlev,ntra)
REAL :: dtrabydt_sh(n_sh,nlev,ntra) ! Increment to tracer due to
                                    ! convection (kg/kg/s)

! Compressed SCM diagnostics for adaptive modset
REAL :: rbuoy_p_out_sh (ScmRowLen*ScmRow,nlev)
REAL :: the_out_sh     (ScmRowLen*ScmRow,nlev)
REAL :: thp_out_sh     (ScmRowLen*ScmRow,nlev)
REAL :: qe_out_sh      (ScmRowLen*ScmRow,nlev)
REAL :: qp_out_sh      (ScmRowLen*ScmRow,nlev)


! Local arrays for (Parcel vertical velocity)^2 [(m/s)^2] on ...
REAL :: w2p_sh(n_sh,nlev) ! ...shallow   convection points
REAL :: w2p_cg(n_cg,nlev) ! ...congestus convection points
REAL :: w2p_md(n_md,nlev) ! ...mid-level convection points
REAL :: w2p_dp(n_dp,nlev) ! ...deep      convection points


!-----------------------------------------------------------------------
! LOCAL compressed arrays to be passed to MID-LEVEL convection scheme.
! Arrays are identified by underscore MD (_md) and are of length
! npnts where npnts is the total number of points in the grid (since
! mid-level convection may occur on any point previously diagnosed as
! shallow or deep).
!-----------------------------------------------------------------------


INTEGER :: midtrig(n_md)   ! Level at which mid level convection
                           ! may start
INTEGER :: mdi(n_md)       ! index for mid points in full grid

INTEGER ::                                                        &
  ntml_md(n_md)                                                   &
 ,ntpar_md(n_md)

REAL ::                                                           &
  pstar_md(n_md)                                                  &
 ,recip_pstar_md(n_md)                                            &
 ,t1_sd_md(n_md)                                                  &
 ,q1_sd_md(n_md)

LOGICAL :: bland_md(n_md)
LOGICAL :: l_mid_md(n_md)     ! true if mid level convection
                              ! compressed version of l_mid_all

REAL ::                 &
 cca_2d_md(n_md)        & ! required by mid level scheme
,cca_md(n_md,n_cca_lev)   ! 3d CCA

REAL ::                                                            &
  p_layer_centres_md(n_md,0:nlev)                                  &
, p_layer_boundaries_md(n_md,0:nlev)                               &
, exner_layer_centres_md(n_md,0:nlev)                              &
, exner_layer_boundaries_md(n_md,0:nlev)                           &
, z_theta_md(n_md,nlev), z_rho_md(n_md,nlev)                       &
, r_theta_md(n_md,0:nlev), r_rho_md(n_md,nlev)                     &
, r2rho_th_md(n_md,nlev)                                           &
, r2rho_md(n_md,nlev)                                              &
, rho_md(n_md,nlev)                                                &
, rho_theta_md(n_md,nlev)                                          &
, dr_across_th_md(n_md,nlev)                                       &
, dr_across_rh_md(n_md,nlev)                                       &
, u_md(n_md,nlev)                                                  &
, v_md(n_md,nlev)                                                  &
, w_md(n_md,nlev)                                                  &
, th_md(n_md,nlev)                                                 &
, q_md(n_md,nlev)

! output variables

REAL ::                                                            &
  dthbydt_md(n_md,nlev)                                            &
, dqbydt_md(n_md,nlev)                                             &
, dubydt_md(n_md,nlev)                                             &
, dvbydt_md(n_md,nlev)

REAL ::                                                            &
  rain_md(n_md)                                                    &
, snow_md(n_md)                                                    &
, rain_3d_md(n_md,nlev)                                            &
, snow_3d_md(n_md,nlev)                                            &
, uw_mid_md(n_md,nlev)                                             &
, vw_mid_md(n_md,nlev)                                             &
, tcw_md(n_md)                                                     &
, cclwp_md(n_md)                                                   &
, ccw_md(n_md,nlev)                                                &
, lcca_md(n_md)                                                    &
, cape_out_md(n_md)                                                &
, cfl_limited_md(n_md)

INTEGER ::                                                         &
  iccb_md(n_md)                                                    &
, icct_md(n_md)                                                    &
, lcbase_md(n_md)                                                  &
, lctop_md(n_md)                                                   &
, freeze_lev_md(n_md)

REAL ::                                                            &
  up_flux_md(n_md,nlev)                                            &
, up_flux_half_md(n_md,nlev)                                       &
, dwn_flux_md(n_md,nlev)                                           &
, entrain_up_md(n_md,nlev)                                         &
, detrain_up_md(n_md,nlev)                                         &
, entrain_dwn_md(n_md,nlev)                                        &
, detrain_dwn_md(n_md,nlev)                                        &
, w_max_md(n_md)

! PC2
REAL ::                                                           &
  qcl_md(n_md,nlev)                                               &
, qcf_md(n_md,nlev)                                               &
, cf_liquid_md(n_md,nlev)                                         &
, cf_frozen_md(n_md,nlev)                                         &
, bulk_cf_md(n_md,nlev)                                           &
, qse_md(n_md,nlev)

REAL ::                                                           &
  dqclbydt_md(n_md,nlev)                                          &
, dqcfbydt_md(n_md,nlev)                                          &
, dcflbydt_md(n_md,nlev)                                          &
, dcffbydt_md(n_md,nlev)                                          &
, dbcfbydt_md(n_md,nlev)

! tracers
REAL :: tracer_md(n_md,trlev,ntra)
REAL :: dtrabydt_md(n_md,nlev,ntra) ! Increment to tracer due to
                                     ! convection (kg/kg/s)

REAL :: rbuoy_p_out_md (ScmRowLen*ScmRow,nlev)
REAL :: the_out_md     (ScmRowLen*ScmRow,nlev)
REAL :: thp_out_md     (ScmRowLen*ScmRow,nlev)
REAL :: qe_out_md      (ScmRowLen*ScmRow,nlev)
REAL :: qp_out_md      (ScmRowLen*ScmRow,nlev)

! Compressed SCM diagnostics for adaptive modset
REAL :: rbuoy_p_out (ScmRowLen*ScmRow,nlev)
REAL :: the_out     (ScmRowLen*ScmRow,nlev)
REAL :: thp_out     (ScmRowLen*ScmRow,nlev)
REAL :: qe_out      (ScmRowLen*ScmRow,nlev)
REAL :: qp_out      (ScmRowLen*ScmRow,nlev)

! Copies so we can see if mid convection modifies the profile
REAL :: rbuoy_p_out_copy (ScmRowLen*ScmRow,nlev)
REAL :: the_out_copy     (ScmRowLen*ScmRow,nlev)
REAL :: thp_out_copy     (ScmRowLen*ScmRow,nlev)
REAL :: qe_out_copy      (ScmRowLen*ScmRow,nlev)
REAL :: qp_out_copy      (ScmRowLen*ScmRow,nlev)


INTEGER :: errorstatus       ! error status

CHARACTER (LEN=errormessagelength) ::                                       &
  cmessage                   ! error message

! added for deep flag
LOGICAL ::      &
  l_deep           ! indicates deep convection

REAL ::         &
  decay_amount  &  ! fraction to reduce deep flag by
 ,cld_depth        ! cloud depth


! should be in a include file
REAL, PARAMETER :: deep_depth=2500.0     ! Minimum depth to be classed
                                        ! as a deep cloud

!-----------------------------------------------------------------------
! Mixing ratio / specific humidity conversions
!-----------------------------------------------------------------------
!  The full model holds q, qcl and qcf, all of which should be used in the
!  conversions and evaluation of qt. (High resolution versions can also
!  hold other moist variables.) The current mass flux scheme without
!  PC2 assume qt=q i.e. it ignores the presence of qcl and qcf.
!  This assumption is also true of the new shallow turbulence based scheme.
!
!  Let m indicate mixing ratio and q specific then
!
!  Now taking account of all moist variables
!      qt = qv + qcl + qcf + (qrain + qgraup + qcf2)
!      mt = mv + mcl + mcf + (mrain + mgraup + mcf2)
!
! Note there could be problems if qcf2 is present as in this case qcf holds
! the snow like ice and qcf2 the ice crystals so what is "seen" and incremented
! by convection (i.e. qcf) could be inappropriate.
!
!  Not PC2
!      dqcl = 0.0,      dqcf = 0.0
!      dmcl = 0.0,      dmcf = 0.0
!
!  Before scheme
!  mv(t) = qv(t)/(1-qt(t))          or qv(t) = mv(t)/(1+mt(t))
!
!  after scheme
!  mv(t+dt) = qv(t+dt)/(1-qt(t+dt)) or qv(t+dt) = mv(t+dt)/(1+mt(t+dt))
!
!  where
!  dmv= mv(t+dt) -  mv(t)           or dqv  =  qv(t+dt)- qv(t)
!
! How to convert increments
!
!  dmv = (dqv + qv*dqt - qt*dqv)/[(1-qt)(1-qt-dqt)]
!
!   or
!
!  dqv = (dmv + mt*dmv - mv*dmt)/[(1+mt)(1+mt+dmt)]
!
! Conversion of qsat
! ------------------
! rsat - mixing ratio saturation value
!
! qsat=rsat/(1+rsat)      and rsat=qsat(1-qsat)
!
!-----------------------------------------------------------------------
!  Arrays required for mixing ratio to specific conversions

REAL :: qt(npnts,nlev)    ! total water variable either mixing
                          ! ratio or specific humidity (kg/kg)
                          ! qt = q + qcl + qcf + (qrain+qgraup+qcf2)
                          ! ( ) extra moist variables not in most runs.
REAL ::                &
  dqt                  &  ! total moisture increment per timestep
, dqtt                 &  ! total moisture increment per s
, denom                   ! 1./denominator

!-----------------------------------------------------------------------
!  Arrays required by 3d CCA calculation

REAL :: cca_2d(npnts)        ! 2d convective cloud (section 5)
REAL :: cca0_2d(npnts)       ! 2d convective cloud (section 0)

REAL :: tcw(npnts)           ! total column water on all points

! Saturation mixing ratio calulation and 1/pstar

REAL ::                 &
  recip_pstar(npnts)    & ! Reciprocal of pstar array
, qse(npnts,nlev)       & ! Saturation specific q of cloud
                          ! cloud environment (kg/kg)
, pt(npnts)             & ! Temporary store for P in calc. of saturation
                          ! value. (Pa)
, tt(npnts)             & ! Temporary store for T in calc.
                          ! of saturation value. (K)
, ttkm1(npnts)            ! Temporary store for T in layer
                          ! k-1 for use in freezing level
                          ! calc. for anvil. (K)
REAL, ALLOCATABLE  ::   &
 qse_mix(:,:)             ! Saturation mixing ratio of cloud
                          ! cloud environment (kg/kg) ONLY used if turbulence
                          ! schemes in use

!  arrays required by grid calculations

REAL ::                    &
 r2rho_th(npnts,nlev)      & ! radius**2 density theta lev (kg/m)
,r2rho(npnts,nlev)         & ! radius**2 density rho lev (kg/m)
,dr_across_th(npnts,nlev)  & ! thickness of theta levels (m)
,dr_across_rh(npnts,nlev)    ! thickness of rho levels (m)

!  arrays required by for identifying convective points

INTEGER :: index1(npnts)
INTEGER :: nconv_all         ! Total number of points convecting

! array required by tracers at end

REAL :: limited_step(npnts)     ! Reduced step size for tracer
                                ! mixing

REAL :: step_test(npnts)        ! Work array used in reducing step


! Parameters


REAL, PARAMETER :: qmin = 1.0e-8 ! Global minimum allowed Q

REAL, PARAMETER :: safety_margin = 1.0e-100 ! Small number used in
                                ! tracer step reduction

!---------------------------------------------------------------------

! Loop counters


INTEGER :: i,j,k,ktra


REAL :: TmpScm1d(ScmRowLen*ScmRow)
REAL :: TmpScm2d(ScmRowLen*ScmRow,nlev+1)
CHARACTER (LEN=15)  :: sname
CHARACTER (LEN=100) :: lname
CHARACTER (LEN=11)  :: units


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Initialise output variables not already initialised
!-----------------------------------------------------------------------

! Initialise 3d rain variables

DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO i=1, n_sh
    rain_3d_sh(i,k) = 0.0
    snow_3d_sh(i,k) = 0.0
  END DO
END DO

DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO i=1, n_dp
    rain_3d_dp(i,k) = 0.0
    snow_3d_dp(i,k) = 0.0
  END DO
END DO

DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO i=1, n_cg
    rain_3d_cg(i,k) = 0.0
    snow_3d_cg(i,k) = 0.0
  END DO
END DO

DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO i=1, n_md
    rain_3d_md(i,k) = 0.0
    snow_3d_md(i,k) = 0.0
  END DO
END DO

DO i = 1,npnts
  wstar_dn(i)                = wstar(i)   ! set to BL value
  ind_cape_reduced(i)        = 0.0
  cape_ts_used(i)            = 0.0
  ind_deep(i)                = 0.0
  ind_shall(i)               = 0.0
END DO

IF (model_type == mt_single_column) THEN
  ! initialize adaptive en/detrainment variables
  ! full arrays initialized to environmental profile
  DO i=1, npnts
    DO k=1, nlev
      the_out (i,k) = th (i,k)
      thp_out (i,k) = th (i,k)
      qe_out  (i,k) = q  (i,k)
      qp_out  (i,k) = q  (i,k)

      the_out_copy (i,k) = th (i,k)
      thp_out_copy (i,k) = th (i,k)
      qe_out_copy  (i,k) = q  (i,k)
      qp_out_copy  (i,k) = q  (i,k)
    END DO
  END DO

  rbuoy_p_out(:,:)      = 0.0
  rbuoy_p_out_copy(:,:) = 0.0

  ! compressed arrays initialised to zero (will be
  ! initialised to compressed versions of environmental
  ! profile within DPCONV, SHCONV and MDCONV
  rbuoy_p_out_dp (:,:) = 0.0
  the_out_dp     (:,:) = 0.0
  thp_out_dp     (:,:) = 0.0
  qe_out_dp      (:,:) = 0.0
  qp_out_dp      (:,:) = 0.0

  rbuoy_p_out_sh (:,:) = 0.0
  the_out_sh     (:,:) = 0.0
  thp_out_sh     (:,:) = 0.0
  qe_out_sh      (:,:) = 0.0
  qp_out_sh      (:,:) = 0.0

  rbuoy_p_out_cg (:,:) = 0.0
  the_out_cg     (:,:) = 0.0
  thp_out_cg     (:,:) = 0.0
  qe_out_cg      (:,:) = 0.0
  qp_out_cg      (:,:) = 0.0

  rbuoy_p_out_md (:,:) = 0.0
  the_out_md     (:,:) = 0.0
  thp_out_md     (:,:) = 0.0
  qe_out_md      (:,:) = 0.0
  qp_out_md      (:,:) = 0.0

END IF ! model_type


!==========================================
! Initialise local arrays for w-eqn
!
IF (flg_w_eqn .OR. l_new_dd) THEN
  DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_sh
      w2p_sh(i,k) = 0.0
    END DO
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_md
      w2p_md(i,k) = 0.0
    END DO
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_dp
      w2p_dp(i,k) = 0.0
    END DO
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_cg
      w2p_cg(i,k) = 0.0
    END DO
  END DO
END IF
!
!==========================================

 ! compressed arrays initialised to zero (will be
 ! initialised to compressed versions of environmental
 ! profile within DPCONV, SHCONV and MDCONV

DO k = 1,nlev
  DO i = 1,npnts
    dqbydt(i,k)              = 0.0
    dthbydt(i,k)             = 0.0
    dqclbydt(i,k)            = 0.0
    dqcfbydt(i,k)            = 0.0
    dcflbydt(i,k)            = 0.0
    dcffbydt(i,k)            = 0.0
    dbcfbydt(i,k)            = 0.0
    ccw(i,k)                 = 0.0
    cca(i,k)                 = 0.0
    dt_dd(i,k)               = 0.0
    dq_dd(i,k)               = 0.0
    du_dd(i,k)               = 0.0
    dv_dd(i,k)               = 0.0
    area_ud(i,k)             = 0.0
    area_dd(i,k)             = 0.0
  END DO
END DO

! Initialisation of qt used for mixing/specific conversions
! Taking account of all possible moist variables even if not use in the
! convection scheme

DO k=1,nlev
  DO i = 1,npnts
    qt(i,k) = q(i,k) + qcl(i,k) + qcf(i,k)
  END DO
END DO

IF (l_mcr_qrain) THEN
  DO k=1,nlev
    DO i = 1,npnts
      qt(i,k) = qt(i,k) + qrain(i,k)
    END DO
  END DO
END IF
IF (l_mcr_qgraup) THEN
  DO k=1,nlev
    DO i = 1,npnts
      qt(i,k) = qt(i,k) + qgraup(i,k)
    END DO
  END DO
END IF
IF (l_mcr_qcf2) THEN
  DO k=1,nlev
    DO i = 1,npnts
      qt(i,k) = qt(i,k) + qcf2(i,k)
    END DO
  END DO
END IF


IF (l_mom) THEN
  DO k=1, nlev
    DO i=1, npnts
      dubydt(i,k)              = 0.0
      dvbydt(i,k)              = 0.0
    END DO
  END DO
END IF
IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k = 1,nlev
      DO i = 1,npnts
        dtrabydt(i,k,ktra)  = 0.0   ! tracer increments
      END DO
    END DO
  END DO
END IF

! Initialise diagnostics if requested

IF (flg_up_flx) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux(i,k)             = 0.0
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      up_flux_half(i,k)        = 0.0
    END DO
  END DO
END IF
IF (flg_dwn_flx) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dwn_flux(i,k)            = 0.0
    END DO
  END DO
END IF
IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_up(i,k)          = 0.0
    END DO
  END DO
END IF
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_up(i,k)          = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      entrain_dwn(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      detrain_dwn(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_mf_deep) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      mf_deep(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_mf_congest) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      mf_congest(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_mf_shall) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      mf_shall(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_mf_midlev) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      mf_midlev(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dt_deep) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dt_deep(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dt_congest) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dt_congest(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dt_shall) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dt_shall(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dt_midlev) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dt_midlev(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dq_deep) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dq_deep(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dq_congest) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dq_congest(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dq_shall) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dq_shall(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dq_midlev) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      dq_midlev(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_du_deep) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      du_deep(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_du_congest) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      du_congest(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_du_shall) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      du_shall(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_du_midlev) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      du_midlev(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dv_deep) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      dv_deep(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dv_congest) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      dv_congest(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dv_shall) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      dv_shall(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_dv_midlev) THEN
  DO k = 1,nlev+1
    DO i = 1,npnts
      dv_midlev(i,k)         = 0.0
    END DO
  END DO
END IF

IF (flg_wqt_flux) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      wqt_flux_sh(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_wql_flux) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      wql_flux_sh(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_wthetal_flux) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      wthetal_flux_sh(i,k)         = 0.0
    END DO
  END DO
END IF
IF (flg_wthetav_flux) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      wthetav_flux_sh(i,k)         = 0.0
    END DO
  END DO
END IF

! Extra variable initialisation for safety

DO i = 1,npnts
  tcw(i)      = 0.0
END DO

IF (flg_uw_dp) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      uw_deep(i,k)  = 0.0
    END DO
  END DO
END IF
IF (flg_vw_dp) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      vw_deep(i,k)  = 0.0
    END DO
  END DO
END IF
IF (flg_uw_shall) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      uw_shall(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_vw_shall) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      vw_shall(i,k) = 0.0
    END DO
  END DO
END IF

IF (flg_uw_mid) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      uw_mid(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_vw_mid) THEN
  DO k = 1,nlev
    DO i = 1,npnts
      vw_mid(i,k) = 0.0
    END DO
  END DO
END IF

DO i = 1,npnts
  kterm_shall(i) = 0
END DO

! Required to get same convective cloud as old scheme
! Need to pass values from deep and shallow to mid level.
! Currently conversion of 2d to 3d makes use of total 2d from
! shallow/deep and mid in a column. If in future the conversion
! does not work on a column total this could be removed.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
DO i=1,n_md
  tcw_md(i) = 0.0
  iccb_md(i) = 0
  icct_md(i) = 0
  cca_2d_md(i) = 0.0
END DO


!-----------------------------------------------------------------------
!  grid information
!-----------------------------------------------------------------------
! Calculate quantities involving density etc for future use

DO k=1,nlev
  DO i=1,npnts
    dr_across_rh(i,k) = r_theta(i,k) - r_theta(i,k-1)
    r2rho(i,k)        = r_rho(i,k)*r_rho(i,k)*rho(i,k)
  END DO
END DO

! rho_theta only set for nlev-1

k=1     ! bottom theta level thicker
DO i=1,npnts
  dr_across_th(i,k) = r_rho(i,k+1) - r_theta(i,0)
  r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
END DO

DO k=2,nlev-1
  DO i=1,npnts
    dr_across_th(i,k) = r_rho(i,k+1) - r_rho(i,k)
    r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k)
  END DO
END DO

k=nlev     ! top layer  (hope not used ?
           !             assume density as layer below)
DO i=1,npnts
  dr_across_th(i,k) = r_theta(i,nlev) - r_rho(i,k)
  r2rho_th(i,k)     = r_theta(i,k)*r_theta(i,k)*rho_theta(i,k-1)
END DO

!-----------------------------------------------------------------------
! 1.0 Section to calculate fields used by all convection types.
! Note cheaper to do here than in individual routines.
! Create saturation mixing ratio arrays  & calculate freeze level
!-----------------------------------------------------------------------

! Calculate 1/pstar and initialize freeze_lev array.

DO i = 1,npnts
  recip_pstar(i)=1.0 / pstar(i)
  freeze_lev(i) = 1
END DO

! Loop over levels

DO k = 1,nlev

  ! Find freezing level

  IF (k  ==  1) THEN
    DO i = 1,npnts
      tt(i) = th(i,k) * exner_layer_centres(i,k)
      pt(i) = p_layer_centres(i,k)

      ! Commented out as initialisation sets freeze_lev to 1.
      ! Code left incase altered in future.

      !            If (tt(i)  <   TM) then
      !              freeze_lev(i) = k
      !            End If

    END DO
  ELSE
    DO i = 1,npnts
      ttkm1(i) = tt(i)
      tt(i) = th(i,k) * exner_layer_centres(i,k)
      pt(i) = p_layer_centres(i,k)
      IF (tt(i)  <   tm .AND. ttkm1(i)  >=  tm) THEN
        IF (freeze_lev(i) == 1) THEN
          freeze_lev(i) = k
        END IF
      END IF
    END DO
  END IF


  ! Calculate saturation specific humidity/mixing ratio  lq_mix=.false.
  ! Note currently mass flux scheme work in specific humidity therefore
  ! call qsat_mix asking for specific values

  IF ( l_new_qsat_conv ) THEN
    CALL qsat_new(qse(:,k),tt,pt,npnts)
  ELSE
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qse(1,k),tt,pt,npnts,.FALSE.)
  END IF

  IF (iconv_deep == 2 .OR. (iconv_shallow >= 2) ) THEN
    ! Using a turbulence scheme also require a mixing ratio version
    IF (k == 1) THEN
      ALLOCATE( qse_mix(npnts,nlev) )
    END IF

    IF ( l_new_qsat_conv ) THEN
      CALL qsat_mix_new(qse_mix(:,k),tt,pt,npnts)
    ELSE
      ! DEPENDS ON: qsat_mix
      CALL qsat_mix(qse_mix(1,k),tt,pt,npnts,.TRUE.)
    END IF

  END IF
END DO  ! nlev

!-----------------------------------------------------------------------
! 1.0 DEEP Convection
! 1.1 Compress input variable arrays for deep convection scheme to
!     length n_dp (deep points only)
!-----------------------------------------------------------------------

IF (n_dp  >   0 .AND. iconv_deep >  0 ) THEN

  ALLOCATE(rho_dry_temp(n_dp,nlev))
  ALLOCATE(rho_dry_theta_temp(n_dp,nlev))

  j = 0
  DO i = 1,npnts
    IF (cumulus_bl(i) .AND. .NOT. l_shallow_bl(i) .AND.               &
                                       .NOT. l_congestus(i)) THEN
      j                        = j+1
      dpi(j)                   = i
    END IF
  END DO


  ! In only variables

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO j=1,n_dp
    bland_dp(j)           = bland(dpi(j))
    ntml_dp(j)            = ntml(dpi(j))
    ntpar_dp(j)           = ntpar(dpi(j))
    pstar_dp(j)           = pstar(dpi(j))
    recip_pstar_dp(j)     = recip_pstar(dpi(j))
    q1_sd_dp(j)           = q1_sd(dpi(j))
    t1_sd_dp(j)           = t1_sd(dpi(j))
    uw0_dp(j)             = uw0(dpi(j))
    vw0_dp(j)             = vw0(dpi(j))
    w_max_dp(j)           = w_max(dpi(j))
    wstar_dp(j)           = wstar(dpi(j))
    qsat_lcl_dp(j)        = qsat_lcl(dpi(j))
    zlcl_uv_dp(j)         = zlcl_uv(dpi(j))
    freeze_lev_dp(j)      = freeze_lev(dpi(j))
    delthvu_dp(j)         = delthvu_bl(dpi(j))
    entrain_coef_dp(j)    = entrain_coef(dpi(j))
  END DO

  DO k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j=1,n_dp
      p_layer_centres_dp(j,k)    = p_layer_centres(dpi(j),k)
      p_layer_boundaries_dp(j,k) = p_layer_boundaries(dpi(j),k)

      exner_layer_centres_dp(j,k)    = exner_layer_centres(dpi(j),k)
      exner_layer_boundaries_dp(j,k) = exner_layer_boundaries(dpi(j),k)
      r_theta_dp(j,k)  = r_theta(dpi(j),k)
    END DO
  END DO

  DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j=1,n_dp
      u_dp(j,k)           = u(dpi(j),k)
      v_dp(j,k)           = v(dpi(j),k)
      w_dp(j,k)           = w(dpi(j),k)
      th_dp(j,k)          = th(dpi(j),k)
      z_theta_dp(j,k)     = z_theta(dpi(j),k)
      z_rho_dp(j,k)       = z_rho(dpi(j),k)
      r_rho_dp(j,k)       = r_rho(dpi(j),k)
      r2rho_th_dp(j,k)     = r2rho_th(dpi(j),k)
      r2rho_dp(j,k)        = r2rho(dpi(j),k)
      rho_theta_dp(j,k)    = rho_theta(dpi(j),k)
      rho_dry_theta_temp(j,k) = rho_dry_theta(dpi(j),k)
      rho_dry_temp(j,k)       = rho_dry(dpi(j),k)
      rho_dp(j,k)          = rho(dpi(j),k)
      dr_across_th_dp(j,k) = dr_across_th(dpi(j),k)
      dr_across_rh_dp(j,k) = dr_across_rh(dpi(j),k)
    END DO
  END DO


  IF (iconv_deep == 1) THEN  ! G-R mass flux scheme

    IF (l_mr_physics) THEN  !  Conversion required

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_dp
          q_dp(j,k)   = q(dpi(j),k)  /(1.0+qt(dpi(j),k))
          qse_dp(j,k) = qse(dpi(j),k)   ! No conversion copy specific value
          qcl_dp(j,k) = qcl(dpi(j),k)/(1.0+qt(dpi(j),k))
          qcf_dp(j,k) = qcf(dpi(j),k)/(1.0+qt(dpi(j),k))
        END DO
      END DO

    ELSE         ! input is specific humidity therefore no problems

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_dp
          q_dp(j,k)   = q(dpi(j),k)
          qse_dp(j,k) = qse(dpi(j),k)
          qcl_dp(j,k) = qcl(dpi(j),k)
          qcf_dp(j,k) = qcf(dpi(j),k)
        END DO
      END DO

    END IF      ! Test on l_mr_physics

  ELSE        ! Future schemes to be code in mixing ratio

    IF (l_mr_physics) THEN   ! Input as required

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_dp
          q_dp(j,k)   = q(dpi(j),k)
          qse_dp(j,k) = qse_mix(dpi(j),k) ! copy mixing ratio value
          qcl_dp(j,k) = qcl(dpi(j),k)
          qcf_dp(j,k) = qcf(dpi(j),k)
        END DO
      END DO

    ELSE      !  Conversion required
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_dp
          q_dp(j,k)   = q(dpi(j),k)  /(1.0-qt(dpi(j),k))
          ! Copy mixing ratio value
          qse_dp(j,k) = qse_mix(dpi(j),k)
          qcl_dp(j,k) = qcl(dpi(j),k)/(1.0-qt(dpi(j),k))
          qcf_dp(j,k) = qcf(dpi(j),k)/(1.0-qt(dpi(j),k))
        END DO
      END DO

    END IF    ! Test on l_mr_physics


  END IF   ! test on deep scheme

  ! In/out variables

  DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j=1,n_dp
      cf_liquid_dp(j,k)   = cf_liquid(dpi(j),k)
      cf_frozen_dp(j,k)   = cf_frozen(dpi(j),k)
      bulk_cf_dp(j,k)     = bulk_cf(dpi(j),k)
    END DO
  END DO

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j=1,n_dp
          tracer_dp(j,k,ktra)  = tracer(dpi(j),k,ktra)
        END DO
      END DO
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j=1,n_dp
          dtrabydt_dp(j,k,ktra)  = 0.0
        END DO
      END DO
    END DO
  END IF

  ALLOCATE(dt_dd_temp(n_dp,nlev))
  ALLOCATE(dq_dd_temp(n_dp,nlev))
  ALLOCATE(du_dd_temp(n_dp,nlev))
  ALLOCATE(dv_dd_temp(n_dp,nlev))
  ALLOCATE(area_ud_temp(n_dp,nlev))
  ALLOCATE(area_dd_temp(n_dp,nlev))

  !-----------------------------------------------------------------------
  ! 1.2 Call deep convection code
  !-----------------------------------------------------------------------

  IF (iconv_deep == 1) THEN  ! 4A like Gregory Rowntree deep conv

    CALL deep_conv_5a(                                            &
                 !IN
                 nbl,nlev,ntra,n_cca_lev,n_dp,trlev,              &
                 bland_dp, delthvu_dp, exner_layer_centres_dp,    &
                 exner_layer_boundaries_dp,                       &
                 l_calc_dxek, l_q_interact,                       &
                 l_tracer, ntml_dp, ntpar_dp,                     &
                 pstar_dp,p_layer_centres_dp,                     &
                 p_layer_boundaries_dp,z_theta_dp,z_rho_dp,       &
                 r_theta_dp, r_rho_dp, rho_theta_dp, rho_dp,      &
                 rho_dry_theta_temp, rho_dry_temp,                &
                 r2rho_th_dp, r2rho_dp,                           &
                 dr_across_th_dp, dr_across_rh_dp,                &
                 q_dp,q1_sd_dp,                                   &
                 t1_sd_dp,th_dp,timestep,                         &
                 u_dp,v_dp,w_dp,                                  &
                 uw0_dp,vw0_dp,w_max_dp,wstar_dp,qsat_lcl_dp,     &
                 entrain_coef_dp,                                 &
                 zlcl_uv_dp,freeze_lev_dp,                        &
                 recip_pstar_dp,qse_dp,                           &
                 !INOUT
                 bulk_cf_dp,cf_frozen_dp,cf_liquid_dp,            &
                 qcf_dp,qcl_dp,tracer_dp,w2p_dp,dtrabydt_dp,      &
                 !OUT
                 cape_out_dp,cclwp_dp,ccw_dp,cca_dp,              &
                 dbcfbydt_dp,dcffbydt_dp,dcflbydt_dp,             &
                 dqbydt_dp,dqcfbydt_dp,dqclbydt_dp,               &
                 dthbydt_dp,dubydt_dp,dvbydt_dp,                  &
                 detrain_up_dp,detrain_dwn_dp,entrain_up_dp,      &
                 entrain_dwn_dp,                                  &
                 iccb_dp,icct_dp,lcca_dp,lcbase_dp,lctop_dp,      &
                 rain_dp,snow_dp, rain_3d_dp, snow_3d_dp,         &
                 up_flux_dp, up_flux_half_dp,                     &
                 dwn_flux_dp,uw_deep_dp,vw_deep_dp,kterm_dp,      &
                 tcw_dp,cca_2d_dp,                                &
                 rbuoy_p_out_dp,the_out_dp,thp_out_dp,            &
                 qe_out_dp,qp_out_dp,ind_cape_reduced_dp,         &
                 cape_ts_used_dp,cfl_limited_dp,ind_deep_dp,      &
                 dt_dd_temp, dq_dd_temp, du_dd_temp, dv_dd_temp,  &
                 area_ud_temp, area_dd_temp,                      &
                 error_point                                      &
                 )

    IF (error_point /= 0) THEN
      errorstatus = 2   ! will cause model to fail
      WRITE(cmessage,'(a37,i12,a8,i3,a9,i2)')                      &
      'Deep conv went to model top at point ',                     &
       dpi(error_point),' in seg ',seg_num,' on call ',call_number

      CALL ereport(routinename, errorstatus, cmessage )

    END IF

  ELSE IF (iconv_deep == 2) THEN ! Turbulence based deep convection

    ! only required by deep turbulence scheme
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j = 1,n_dp
      wthvs_dp(j) = wthvs(dpi(j))
      wth0_dp(j)  = ftl(dpi(j))/rho(dpi(j),1)/cp
      wq0_dp(j)   = fqt(dpi(j))/rho(dpi(j),1)
      zlcl_dp(j)  = zlcl(dpi(j))
      ftl_dp(j)   = ftl(dpi(j))
      fqt_dp(j)   = fqt(dpi(j))
    END DO

    CALL deep_turb_conv                                                 &
     ! Intent IN
     ( call_number, nbl, nlev, n_cca_lev, ntra, n_dp, trlev             &
     , ntml_dp, ntpar_dp, freeze_lev_dp                                 &
     , l_calc_dxek, l_q_interact, l_tracer                              &
     , bland_dp                                                         &
     , timestep, pstar_dp, p_layer_centres_dp, p_layer_boundaries_dp    &
     , exner_layer_centres_dp, exner_layer_boundaries_dp, z_theta_dp    &
     , z_rho_dp, dr_across_th_dp, dr_across_rh_dp, rho_dp, rho_theta_dp &
     , r_rho_dp, r_theta_dp, r2rho_dp, r2rho_th_dp                      &
     , qse_dp, th_dp, q_dp, qcl_dp, qcf_dp                              &
     , bulk_cf_dp, cf_frozen_dp, cf_liquid_dp                           &
     , u_dp, v_dp, w_dp                                                 &
     , uw0_dp, vw0_dp, fqt_dp, ftl_dp                                   &
     , zlcl_dp, zlcl_uv_dp, wth0_dp, wq0_dp, wstar_dp                   &
     , wthvs_dp, q1_sd_dp,t1_sd_dp                                      &
     , w_max_dp                                                         &

     ! Intent INOUT
     , tracer_dp                                                        &

     ! Intent OUT
     , kterm_dp, iccb_dp, icct_dp                                       &
     , dthbydt_dp, dqbydt_dp, dqclbydt_dp, dqcfbydt_dp, dbcfbydt_dp     &
     , dcffbydt_dp, dcflbydt_dp, dubydt_dp, dvbydt_dp, dtrabydt_dp      &
     , rain_dp, snow_dp, rain_3d_dp, snow_3d_dp, up_flux_dp, cca_2d_dp  &
     , cca_dp, ccw_dp, cclwp_dp                                         &
     , uw_deep_dp,vw_deep_dp, mb_dp )

    ! Not set by new scheme therefore zero to ensure safety and reproducible
    ! answers from mid.
    lcbase_dp(:) = iccb_dp(:)
    lctop_dp(:)  = icct_dp(:)
    tcw_dp(:) = 0.0
    ind_cape_reduced_dp(:) = 0.0
    cape_ts_used_dp(:) = 0.0

  END IF

  !-----------------------------------------------------------------------
  ! 1.3 Write data from deep convection points to full arrays
  !-----------------------------------------------------------------------
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO i = 1,n_dp
    cape_out(dpi(i))           = cape_out_dp(i)
    cclwp(dpi(i))              = cclwp_dp(i)
    iccb(dpi(i))               = iccb_dp(i)
    icct(dpi(i))               = icct_dp(i)

    lcca(dpi(i))               = lcca_dp(i)
    lcbase(dpi(i))             = lcbase_dp(i)
    lctop(dpi(i))              = lctop_dp(i)

    rain(dpi(i))               = rain_dp(i)
    snow(dpi(i))               = snow_dp(i)
    precip_deep(dpi(i))        = rain_dp(i) + snow_dp(i)
    kterm_deep(dpi(i))         = kterm_dp(i)
    ind_cape_reduced(dpi(i))   = ind_cape_reduced_dp(i)
    cape_ts_used(dpi(i))       = cape_ts_used_dp(i)
    deep_cfl_limited(dpi(i))   = cfl_limited_dp(i)
    ind_deep(dpi(i))           = ind_deep_dp(i)

    tcw   (dpi(i)) = tcw_dp   (i)
    cca_2d(dpi(i)) = cca_2d_dp(i)
  END DO

  IF (l_mom) THEN
    DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        dubydt(dpi(i),k)         = dubydt_dp(i,k)
        dvbydt(dpi(i),k)         = dvbydt_dp(i,k)
        du_dd(dpi(i),k)          = du_dd_temp(i,k)
        dv_dd(dpi(i),k)          = dv_dd_temp(i,k)
      END DO
    END DO
  END IF

  IF (flg_w_eqn) THEN
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_dp
        w2p(dpi(i),k) = w2p_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_area_ud) THEN
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_dp
        area_ud(dpi(i),k)        = area_ud_temp(i,k)
      END DO
    END DO
  END IF
  IF (flg_area_dd) THEN
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_dp
        area_dd(dpi(i),k)        = area_dd_temp(i,k)
      END DO
    END DO
  END IF

  DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i = 1,n_dp
      dqbydt(dpi(i),k)         = dqbydt_dp(i,k)
      dthbydt(dpi(i),k)        = dthbydt_dp(i,k)
      dqclbydt(dpi(i),k)       = dqclbydt_dp(i,k)
      dqcfbydt(dpi(i),k)       = dqcfbydt_dp(i,k)
      dcflbydt(dpi(i),k)       = dcflbydt_dp(i,k)
      dcffbydt(dpi(i),k)       = dcffbydt_dp(i,k)
      dbcfbydt(dpi(i),k)       = dbcfbydt_dp(i,k)
      ccw(dpi(i),k)            = ccw_dp(i,k)
      dt_dd(dpi(i),k)          = dt_dd_temp(i,k)
    END DO
  END DO

  DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_dp
      cca(dpi(i),k) = cca_dp(i,k)
    END DO
  END DO


  IF (l_ccrad) THEN

    ! Then scaling is done by CCRad, even if PC2 = .true.
    ! To zero CCA (section 0) for PC2 use the CCRAD knobs
    ! cca_knobs should be set to 0.0 in gui/namelist for PC2

    ! Use of Anvils will affect both section 0/5 diags if ccrad knobs
    ! are not set to 0.0, this may required further attention in future
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_dp
        ccw0(dpi(i),k) = ccw_dp(i,k) * ccw_dp_knob
      END DO
    END DO

    DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_dp
        cca0(dpi(i),k) = cca_dp(i,k) * cca_dp_knob
      END DO
    END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_dp
      cclwp0 (dpi(i)) = cclwp  (dpi(i))
      iccb0  (dpi(i)) = iccb   (dpi(i))
      icct0  (dpi(i)) = icct   (dpi(i))
      lcbase0(dpi(i)) = lcbase (dpi(i))
      !         cca0_2d(dpi(i)) = cca0_2d(dpi(i))
    END DO

  ELSE          ! Non-CCRad Code

    IF (l_q_interact) THEN
      ! PC2 Code (PC2 has no deep contribution)
    ELSE
      ! Non-PC2 Code
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_dp
        cclwp0 (dpi(i)) = cclwp (dpi(i))
        iccb0  (dpi(i)) = iccb  (dpi(i))
        icct0  (dpi(i)) = icct  (dpi(i))
        lcbase0(dpi(i)) = lcbase(dpi(i))
      END DO

      DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_dp
          cca0(dpi(i),k) = cca_dp(i,k)
        END DO
      END DO
    END IF

  END IF        ! l_ccrad

  IF (l_cca_dp_prog) THEN
    DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_dp
        cca0_dp(dpi(i),k) = cca_dp(i,k)
      END DO
    END DO
  END IF


  IF (model_type == mt_single_column) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        rbuoy_p_out(dpi(i),k) = rbuoy_p_out_dp(i,k)
        the_out(dpi(i),k)     = the_out_dp(i,k)
        thp_out(dpi(i),k)     = thp_out_dp(i,k)
        qe_out(dpi(i),k)      = qe_out_dp(i,k)
        qp_out(dpi(i),k)      = qp_out_dp(i,k)
      END DO
    END DO
  END IF ! model_type

  IF (iconv_deep == 1) THEN  !G-R mass flux scheme
    IF (l_mr_physics) THEN ! requires conversion

      IF (l_q_interact) THEN   ! PC2

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_dp
            dqtt  = dqbydt_dp(i,k)+dqclbydt_dp(i,k)+dqcfbydt_dp(i,k)
            dqt   = dqtt*timestep
            denom = 1.0/((1.0-qt(dpi(i),k))*(1.0-qt(dpi(i),k)-dqt))
            dqbydt(dpi(i),k)   = denom *                               &
                ( dqbydt_dp(i,k)*(1.0-qt(dpi(i),k))+q_dp(i,k)*dqtt )
            dqclbydt(dpi(i),k) = denom *                               &
              ( dqclbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcl_dp(i,k)*dqtt )
            dqcfbydt(dpi(i),k) = denom *                               &
              ( dqcfbydt_dp(i,k)*(1.0-qt(dpi(i),k))+qcf_dp(i,k)*dqtt )
            dq_dd(dpi(i),k)     = denom *                               &
                ( dq_dd_temp(i,k)*(1.0-qt(dpi(i),k))+q_dp(i,k)*dqtt )
          END DO
        END DO

      ELSE                     ! Not PC2
                               ! No qcl and qcf increments anyway
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_dp
            dqt   = dqbydt_dp(i,k)*timestep
            denom = 1.0/((1.0-qt(dpi(i),k))*(1.0-qt(dpi(i),k)-dqt))
            dqbydt(dpi(i),k)    = dqbydt_dp(i,k)*denom
            dqclbydt(dpi(i),k)  = 0.0
            dqcfbydt(dpi(i),k)  = 0.0
            dq_dd(dpi(i),k)     = dq_dd_temp(i,k)*denom
          END DO
        END DO

      END IF                   ! end test on PC2

    ELSE   ! output is specific humidity therefore no problems

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_dp
          dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
          dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
          dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
          dq_dd(dpi(i),k)    = dq_dd_temp(i,k)
        END DO
      END DO

    END IF      ! Test on l_mr_physics

  ELSE         ! Future schemes output as  mixing ratio

    IF (l_mr_physics) THEN ! ok
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_dp
          dqbydt(dpi(i),k)   = dqbydt_dp(i,k)
          dqclbydt(dpi(i),k) = dqclbydt_dp(i,k)
          dqcfbydt(dpi(i),k) = dqcfbydt_dp(i,k)
          dq_dd(dpi(i),k)    = dq_dd_temp(i,k)
        END DO
      END DO

    ELSE   ! output is specific humidity

      IF (l_q_interact) THEN      ! PC2

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_dp
            dqtt= dqbydt_dp(i,k)+dqclbydt_dp(i,k)+dqcfbydt_dp(i,k)
            dqt = dqtt*timestep
            denom = 1.0/((1.0+qt(dpi(i),k))*(1.0+qt(dpi(i),k)+dqt))
            dqbydt(dpi(i),k)   = denom *                             &
               (dqbydt_dp(i,k)*(1.0+qt(dpi(i),k))-q_dp(i,k)*dqtt)
            dqclbydt(dpi(i),k) = denom *                             &
             (dqclbydt_dp(i,k)*(1.0+qt(dpi(i),k))-qcl_dp(i,k)*dqtt)
            dqcfbydt(dpi(i),k) = denom *                             &
             (dqcfbydt_dp(i,k)*(1.0+qt(dpi(i),k))-qcf_dp(i,k)*dqtt)

            dq_dd(dpi(i),k)    = denom *                            &
               (dq_dd_temp(i,k)*(1.0+qt(dpi(i),k))-q_dp(i,k)*dqtt)
          END DO
        END DO

      ELSE                         ! Not PC2
                                   ! No qcl and qcf increments
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_dp
            dqt   = dqbydt_dp(i,k)*timestep
            denom = 1.0/((1.0+qt(dpi(i),k))*(1.0+qt(dpi(i),k)+dqt))
            dqbydt(dpi(i),k)   =  dqbydt_dp(i,k)*denom
            dqclbydt(dpi(i),k) = 0.0
            dqcfbydt(dpi(i),k) = 0.0
            dq_dd(dpi(i),k)    = denom *dq_dd_temp(i,k)
          END DO
        END DO

      END IF                       ! end test on PC2

    END IF      ! Test on l_mr_physics
  END IF        ! test on deep scheme

  DEALLOCATE(area_dd_temp)
  DEALLOCATE(area_ud_temp)
  DEALLOCATE(dv_dd_temp)
  DEALLOCATE(du_dd_temp)
  DEALLOCATE(dq_dd_temp)
  DEALLOCATE(dt_dd_temp)
  DEALLOCATE(rho_dry_theta_temp)
  DEALLOCATE(rho_dry_temp)

  IF (flg_up_flx) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        up_flux(dpi(i),k)        = up_flux_dp(i,k)
      END DO
    END DO
  END IF
  ! No values set for turbulence option
  IF (flg_up_flx_half .AND. iconv_deep == 1) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        up_flux_half(dpi(i),k)        = up_flux_half_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_dwn_flx) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        dwn_flux(dpi(i),k)       = dwn_flux_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_up) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        entrain_up(dpi(i),k)     = entrain_up_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_detr_up) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        detrain_up(dpi(i),k)     = detrain_up_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_entr_dwn) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        entrain_dwn(dpi(i),k)    = entrain_dwn_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_detr_dwn) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        detrain_dwn(dpi(i),k)    = detrain_dwn_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_uw_dp) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        uw_deep(dpi(i),k)        = uw_deep_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_vw_dp) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        vw_deep(dpi(i),k)        = vw_deep_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_mf_deep) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        mf_deep(dpi(i),k)        = up_flux_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_dt_deep) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        dt_deep(dpi(i),k)        = dthbydt_dp(i,k)
      END DO
    END DO
  END IF
  IF (flg_dq_deep) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_dp
        dq_deep(dpi(i),k)        = dqbydt_dp(i,k)
      END DO
    END DO
  END IF
  IF (l_mom) THEN
    IF (flg_du_deep) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_dp
          du_deep(dpi(i),k)        = dubydt_dp(i,k)
        END DO
      END DO
    END IF
    IF (flg_dv_deep) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_dp
          dv_deep(dpi(i),k)        = dvbydt_dp(i,k)
        END DO
      END DO
    END IF
  END IF

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1, n_dp
          dtrabydt(dpi(i),k,ktra)  = dtrabydt_dp(i,k,ktra)
        END DO
      END DO
    END DO
  END IF

  ! Merge dp 3d rain & snow profiles

  DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_dp
      rain_3d(dpi(i),k) = rain_3d_dp(i,k)
      snow_3d(dpi(i),k) = snow_3d_dp(i,k)
    END DO
  END DO

  !----------------------------------------------------------------------
  ! Setting of deep convective history flag
  ! Only class as deep convection if depth of convection greater than
  ! a set value.
  !----------------------------------------------------------------------

  IF (l_conv_hist) THEN

    decay_amount= timestep/decay_period

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i = 1,n_dp
      cld_depth = z_rho_dp(i,kterm_dp(i)+1) - z_rho_dp(i,ntml_dp(i)+1)

      ! Store depth of convection
      past_conv_ht(dpi(i)) = cld_depth

      IF (cld_depth > deep_depth) THEN
        deep_flag(dpi(i)) = 1.0      ! deep convection occurred
      ELSE   ! failed deep convection
        deep_flag(dpi(i)) = deep_flag(dpi(i)) - decay_amount
        deep_flag(dpi(i)) = MAX(deep_flag(dpi(i)),0.0)
      END IF
    END DO       ! loop over deep poins

  END IF

END IF   ! test on n_dp >0


!-----------------------------------------------------------------------
! If not deep decay deep flag
!-----------------------------------------------------------------------

IF (l_conv_hist) THEN

  decay_amount= timestep/decay_period

  DO i = 1,npnts
    l_deep=cumulus_bl(i) .AND. .NOT. l_shallow_bl(i) .AND.               &
                                      .NOT. l_congestus(i)
    ! Decay flag ensuring 0.0 is the minimum allowed value.
    IF (.NOT. l_deep) THEN
      deep_flag(i) = deep_flag(i) - decay_amount
      deep_flag(i) = MAX(deep_flag(i),0.0)

      ! Do I want to decay old convective depth here? At present NO
      !            past_conv_ht(i) = (1.0-decay_amount) * past_conv_ht(i)

    END IF
  END DO

END IF    ! test on l_conv_hist

!-----------------------------------------------------------------------
! 2.0 SHALLOW convection
! 2.1 Compress input variable arrays for shallow convection scheme to
!     length n_sh (shallow points only)
!-----------------------------------------------------------------------

IF (n_sh  >   0 .AND. iconv_shallow >  0) THEN

  ALLOCATE(rho_dry_temp(n_sh,nlev))
  ALLOCATE(rho_dry_theta_temp(n_sh,nlev))

  j = 0
  DO i = 1,npnts
    IF (cumulus_bl(i) .AND. l_shallow_bl(i)) THEN
      j                        = j+1
      shi(j)                   = i
    END IF
  END DO

  ! In only variables

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO j = 1,n_sh
    bland_sh(j)           = bland(shi(j))
    delthvu_sh(j)         = delthvu_bl(shi(j))
    ql_ad_sh(j)           = ql_ad(shi(j))
    ntml_sh(j)            = ntml(shi(j))
    ntpar_sh(j)           = ntpar(shi(j))
    pstar_sh(j)           = pstar(shi(j))
    recip_pstar_sh(j)     = recip_pstar(shi(j))
    q1_sd_sh(j)           = q1_sd(shi(j))
    t1_sd_sh(j)           = t1_sd(shi(j))
    uw0_sh(j)             = uw0(shi(j))
    vw0_sh(j)             = vw0(shi(j))
    wstar_sh(j)           = wstar(shi(j))
    wthvs_sh(j)           = wthvs(shi(j))
    zlcl_uv_sh(j)         = zlcl_uv(shi(j))
    ztop_uv_sh(j)         = ztop_uv(shi(j))
    freeze_lev_sh(j)      = freeze_lev(shi(j))
    entrain_coef_sh(j)    = entrain_coef(shi(j))
    delta_smag_sh(j)      = delta_smag(shi(j))

    ! initialise to zero as not used in every option

    wstar_up_sh(j) = 0.0
    mb1_sh(j) = 0.0
    mb2_sh(j) = 0.0
  END DO

  DO k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j = 1,n_sh
      p_layer_centres_sh(j,k)     = p_layer_centres(shi(j),k)
      p_layer_boundaries_sh(j,k)  = p_layer_boundaries(shi(j),k)
      exner_layer_centres_sh(j,k) = exner_layer_centres(shi(j),k)
      exner_layer_boundaries_sh(j,k) = exner_layer_boundaries(shi(j),k)
      r_theta_sh(j,k)             = r_theta(shi(j),k)
    END DO
  END DO

  DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j = 1,n_sh
      u_sh(j,k)           = u(shi(j),k)
      v_sh(j,k)           = v(shi(j),k)
      w_sh(j,k)           = w(shi(j),k)
      th_sh(j,k)          = th(shi(j),k)
      z_theta_sh(j,k)     = z_theta(shi(j),k)
      z_rho_sh(j,k)       = z_rho(shi(j),k)
      rho_sh(j,k)         = rho(shi(j),k)
      rho_th_sh(j,k)      = rho_theta(shi(j),k)
      rho_dry_theta_temp(j,k) = rho_dry_theta(shi(j),k)
      rho_dry_temp(j,k)       = rho_dry(shi(j),k)
      r2rho_sh(j,k)       = r2rho(shi(j),k)
      r2rho_th_sh(j,k)    = r2rho_th(shi(j),k)
      r_rho_sh(j,k)       = r_rho(shi(j),k)
      dr_across_rh_sh(j,k)       = dr_across_rh(shi(j),k)
      dr_across_th_sh(j,k)       = dr_across_th(shi(j),k)
    END DO
  END DO

  ! moisture  - input depends on scheme

  IF (iconv_shallow == 1) THEN     ! G-R scheme

    ! G-R requires input of specific humidity

    IF (l_mr_physics) THEN

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_sh
          q_sh(j,k)   = q(shi(j),k)  /(1.0+qt(shi(j),k))
          qse_sh(j,k) = qse(shi(j),k)  ! copy specific value
          qcl_sh(j,k) = qcl(shi(j),k)/(1.0+qt(shi(j),k))
          qcf_sh(j,k) = qcf(shi(j),k)/(1.0+qt(shi(j),k))
        END DO
      END DO

    ELSE       ! Input is specific humidity therefore no problems

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_sh
          q_sh(j,k)   = q(shi(j),k)
          qse_sh(j,k) = qse(shi(j),k) ! copy specific value
          qcl_sh(j,k) = qcl(shi(j),k)
          qcf_sh(j,k) = qcf(shi(j),k)
        END DO
      END DO

    END IF     ! Test on l_mr_physics

  ELSE         ! turbulence scheme requires mixing ratio

    ! turbulence scheme requires mixing ratios as input

    IF (l_mr_physics) THEN

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_sh
          q_sh(j,k)   = q(shi(j),k)
          qse_sh(j,k) = qse_mix(shi(j),k) ! copy mixing ratio value
          qcl_sh(j,k) = qcl(shi(j),k)
          qcf_sh(j,k) = qcf(shi(j),k)
        END DO
      END DO

    ELSE       ! Require conversion

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_sh
          q_sh(j,k)     = q(shi(j),k)/(1.0 -qt(shi(j),k))

          ! check for negative values
          IF (q_sh(j,k) <  0.0) THEN
            IF (printstatus >= prstatus_normal) THEN
              WRITE(umMessage,'(a18,g26.18,a5,2i6)') 'problem q_mix -ve ',    &
                             q_sh(j,k),' j,k ',j,k
              CALL umPrint(umMessage,src='glue_conv-gconv5a')

            END IF
            q_sh(j,k) = 0.0
          END IF
          qse_sh(j,k) = qse_mix(shi(j),k) ! copy mixing ratio value

          ! Currently turbulence shallow scheme does not alter or use qcl or qcf
          ! so no checks are being applied to the values.
          qcl_sh(j,k) = qcl(shi(j),k)/(1.0-qt(shi(j),k))
          qcf_sh(j,k) = qcf(shi(j),k)/(1.0-qt(shi(j),k))
        END DO  ! j
      END DO    ! k

    END IF      ! Test on l_mr_physics

  END IF        ! test on scheme


  ! In/out variables


  DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j = 1,n_sh
      cf_liquid_sh(j,k)   = cf_liquid(shi(j),k)
      cf_frozen_sh(j,k)   = cf_frozen(shi(j),k)
      bulk_cf_sh(j,k)     = bulk_cf(shi(j),k)
    END DO
  END DO

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_sh
          tracer_sh(j,k,ktra)  = tracer(shi(j),k,ktra)
        END DO
      END DO
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO j = 1,n_sh
          dtrabydt_sh(j,k,ktra)  = 0.0
        END DO
      END DO
    END DO
  END IF

  ALLOCATE(dt_dd_temp(n_sh,nlev))
  ALLOCATE(dq_dd_temp(n_sh,nlev))
  ALLOCATE(du_dd_temp(n_sh,nlev))
  ALLOCATE(dv_dd_temp(n_sh,nlev))
  ALLOCATE(area_ud_temp(n_sh,nlev))
  ALLOCATE(area_dd_temp(n_sh,nlev))

  !-----------------------------------------------------------------------
  ! 2.2 Call shallow convection code
  !-----------------------------------------------------------------------

  IF (iconv_shallow == 1) THEN !  G-R mass flux type scheme

    CALL shallow_conv_5a(                                         &
                     !IN
                     nbl,nlev,ntra,n_cca_lev,n_sh,trlev,          &
                     bland_sh,                                    &
                     delthvu_sh,exner_layer_centres_sh,           &
                     exner_layer_boundaries_sh, l_calc_dxek,      &
                     l_q_interact,l_tracer,ntml_sh,ntpar_sh,      &
                     pstar_sh,p_layer_centres_sh,                 &
                     p_layer_boundaries_sh,z_theta_sh,z_rho_sh,   &
                     r_theta_sh,r_rho_sh,rho_th_sh,               &
                     rho_dry_theta_temp, rho_sh, rho_dry_temp,    &
                     r2rho_th_sh,r2rho_sh,                        &
                     dr_across_th_sh,dr_across_rh_sh,             &
                     q_sh,q1_sd_sh,t1_sd_sh,th_sh,timestep,       &
                     u_sh,v_sh,w_sh,uw0_sh,vw0_sh,wstar_sh,       &
                     wthvs_sh,entrain_coef_sh,delta_smag_sh,      &
                     zlcl_uv_sh,ztop_uv_sh,freeze_lev_sh,         &
                     recip_pstar_sh,qse_sh,                       &
                     !INOUT
                     bulk_cf_sh,cf_frozen_sh,cf_liquid_sh,        &
                     qcf_sh,qcl_sh,tracer_sh,w2p_sh,              &
                     !OUT
                     cape_out_sh,cclwp_sh,ccw_sh,cca_sh,          &
                     dbcfbydt_sh,dcffbydt_sh,dcflbydt_sh,         &
                     dqbydt_sh,dqcfbydt_sh,dqclbydt_sh,           &
                     dthbydt_sh,dubydt_sh,dvbydt_sh,dtrabydt_sh,  &
                     detrain_up_sh,detrain_dwn_sh,entrain_up_sh,  &
                     entrain_dwn_sh,                              &
                     iccb_sh,icct_sh,lcca_sh,lcbase_sh,lctop_sh,  &
                     rain_sh,snow_sh,rain_3d_sh,snow_3d_sh,       &
                     up_flux_sh,up_flux_half_sh,                  &
                     dwn_flux_sh,uw_shall_sh,vw_shall_sh,         &
                     tcw_sh,cca_2d_sh,kterm_sh,ind_shall_sh,      &
                     rbuoy_p_out_sh,the_out_sh,thp_out_sh,        &
                     qe_out_sh,qp_out_sh,dt_dd_temp,dq_dd_temp,   &
                     du_dd_temp,dv_dd_temp,                       &
                     area_ud_temp, area_dd_temp                   &
                     )

    ! Initialisations required by some SCM configurations
    DO k = 1,nlev
      DO j = 1,n_sh
        wqt_sh(j,k)     = 0.0
        wql_sh(j,k)     = 0.0
        wthetal_sh(j,k) = 0.0
        wthetav_sh(j,k) = 0.0
      END DO
    END DO
 
  ELSE IF (iconv_shallow == 2 .OR. iconv_shallow == 3) THEN

    ! turbulence based shallow convection scheme
    !              2 - version with no shallow precipitation
    !              3 - version with shallow precipitation
    IF (iconv_shallow == 2) THEN
      ishall_precip = 0
    ELSE
      ishall_precip = 1
    END IF

    ! Information only required by turbulence version

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j = 1,n_sh
      zlcl_sh(j)            = zlcl(shi(j))
    END DO


    CALL shallow_turb_conv(                                       &
             !IN
                      nlev, n_sh, ntra,trlev                      &
,                     ntml_sh, ntpar_sh                           &
,                     ishall_precip, l_calc_dxek, l_q_interact    &
,                     l_tracer                                    &
,                     bland_sh, timestep                          &
,                     delthvu_sh,ql_ad_sh,uw0_sh,vw0_sh           &
,                     wstar_sh,wthvs_sh,zlcl_sh,zlcl_uv_sh        &
,                     ztop_uv_sh                                  &
,                     pstar_sh,p_layer_centres_sh                 &
,                     p_layer_boundaries_sh                       &
,                     exner_layer_centres_sh                      &
,                     exner_layer_boundaries_sh                   &
,                     z_theta_sh,z_rho_sh,r_rho_sh,r_theta_sh     &
,                     rho_sh,rho_th_sh,r2rho_sh,r2rho_th_sh       &
,                     dr_across_rh_sh,dr_across_th_sh             &
,                     q_sh,th_sh,u_sh,v_sh                        &
,                     qse_sh                                      &
             !INOUT
,                     bulk_cf_sh,cf_frozen_sh,cf_liquid_sh        &
,                     qcf_sh,qcl_sh,tracer_sh                     &
             !OUT
,                     cape_out_sh,cclwp_sh,ccw_sh,cca_sh          &
,                     dbcfbydt_sh,dcffbydt_sh,dcflbydt_sh         &
,                     dqbydt_sh,dqcfbydt_sh,dqclbydt_sh           &
,                     dthbydt_sh,dubydt_sh,dvbydt_sh,dtrabydt_sh  &
,                     iccb_sh,icct_sh,lcca_sh,lcbase_sh,lctop_sh  &
,                     rain_sh,snow_sh,up_flux_sh                  &
,                     uw_shall_sh,vw_shall_sh                     &
,                     wqt_sh,wthetal_sh,wthetav_sh,wql_sh         &
,                     wstar_up_sh,mb1_sh,mb2_sh                   &
,                     tcw_sh,cca_2d_sh                            &
                     )

  END IF        ! test on shallow convection scheme


  !-----------------------------------------------------------------------
  ! 2.3 Write data from shallow convection points to full arrays
  !-----------------------------------------------------------------------
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO i = 1,n_sh
    cape_out(shi(i))           = cape_out_sh(i)
    cclwp(shi(i))              = cclwp_sh(i)
    iccb(shi(i))               = iccb_sh(i)
    icct(shi(i))               = icct_sh(i)

    lcca(shi(i))               = lcca_sh(i)
    lcbase(shi(i))             = lcbase_sh(i)
    lctop(shi(i))              = lctop_sh(i)

    rain(shi(i))               = rain_sh(i)
    snow(shi(i))               = snow_sh(i)
    precip_shall(shi(i))       = rain_sh(i) + snow_sh(i)
    wstar_dn(shi(i))           = wstar_sh(i)  ! wstar_dn
    wstar_up(shi(i))           = wstar_up_sh(i)
    mb1(shi(i))                = mb1_sh(i)
    mb2(shi(i))                = mb2_sh(i)
    ind_shall(shi(i))          = ind_shall_sh(i)

    tcw   (shi(i)) = tcw_sh   (i)
    cca_2d(shi(i)) = cca_2d_sh(i)

    kterm_shall(shi(i)) = kterm_sh(i)
  END DO

  IF (l_mom) THEN
    DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        dubydt(shi(i),k)         = dubydt_sh(i,k)
        dvbydt(shi(i),k)         = dvbydt_sh(i,k)
      END DO
    END DO
  END IF

  IF (flg_w_eqn) THEN
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        w2p(shi(i),k) = w2p_sh(i,k)
      END DO
    END DO
  END IF

  IF (flg_area_ud) THEN
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        area_ud(shi(i),k) = area_ud_temp(i,k)
      END DO
    END DO
  END IF

  IF (flg_area_dd) THEN
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        area_dd(shi(i),k) = area_dd_temp(i,k)
      END DO
    END DO
  END IF

  DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i = 1,n_sh
      dthbydt(shi(i),k)        = dthbydt_sh(i,k)
      dcflbydt(shi(i),k)       = dcflbydt_sh(i,k)
      dcffbydt(shi(i),k)       = dcffbydt_sh(i,k)
      dbcfbydt(shi(i),k)       = dbcfbydt_sh(i,k)
      ccw(shi(i),k)            = ccw_sh(i,k)
      dt_dd(shi(i),k)          = dt_dd_temp(i,k)
      du_dd(shi(i),k)          = du_dd_temp(i,k)
      dv_dd(shi(i),k)          = dv_dd_temp(i,k)
    END DO
  END DO

  DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_sh
      cca(shi(i),k) = cca_sh(i,k)
    END DO
  END DO

  IF (l_ccrad) THEN

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_sh
      iccb0  (shi(i)) = iccb  (shi(i))
      icct0  (shi(i)) = icct  (shi(i))
      cclwp0 (shi(i)) = cclwp (shi(i))
      lcbase0(shi(i)) = lcbase(shi(i))
    END DO

    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        ccw0(shi(i),k) = ccw_sh(i,k) * ccw_sh_knob
      END DO
    END DO

    DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        cca0(shi(i),k) = cca_sh(i,k) * cca_sh_knob
      END DO
    END DO

  ELSE ! Non-CCRad code

    IF ( l_q_interact .AND. (.NOT. l_pc2_diag_sh) ) THEN
      ! Do nothing as PC2 only keeps shallow component
      ! if l_pc2_diag_sh=.TRUE.
    ELSE
      ! Either PC2=.FALSE. OR PC2=.TRUE. and it wants
      ! Shallow cloud
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        iccb0   (shi(i)) = iccb   (shi(i))
        icct0   (shi(i)) = icct   (shi(i))
        cclwp0  (shi(i)) = cclwp  (shi(i))
        lcbase0 (shi(i)) = lcbase (shi(i))
      END DO

      DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_sh
          cca0(shi(i),k) = cca_sh(i,k)
        END DO
      END DO

      !         ! Although ccw0 is not passed up to atm_step
      !         ! because (l_ccrad=.false.), shallow points
      !         ! are kept because they may needed in ni_imp_ctl
      !         ! by l_pc2_diag_sh. This needs looking into/possibly
      !         ! fixing.

      !         DO k=1, nlev
      !           DO i=1, n_sh
      !             ccw0(shi(i),k) = ccw_sh(i,k)
      !           END DO
      !         END DO
      !
      !         DO i=1, n_sh
      !           cca0_2d(shi(i)) = cca_2d_sh(i)
      !         END DO

    END IF

  END IF ! l_ccrad

  IF (l_cca_sh_prog) THEN
    DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        cca0_sh(shi(i),k) = cca_sh(i,k)
      END DO
    END DO
  END IF


  IF (model_type == mt_single_column) THEN
    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_sh
        rbuoy_p_out(shi(i),k) = rbuoy_p_out_sh(i,k)
        the_out(shi(i),k)     = the_out_sh(i,k)
        thp_out(shi(i),k)     = thp_out_sh(i,k)
        qe_out(shi(i),k)      = qe_out_sh(i,k)
        qp_out(shi(i),k)      = qp_out_sh(i,k)
      END DO
    END DO
  END IF


  ! Set past convective depth BUT only if dump holds space for history
  ! prognostics

  IF (l_conv_hist) THEN
    IF (iconv_shallow == 1) THEN

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        cld_depth = z_rho_sh(i,kterm_sh(i)+1) - z_rho_sh(i,ntml_sh(i)+1)

        ! In some cases shallow fails and kterm remains set to zero so cld_depth
        ! will be negative.
        IF ( cld_depth < 0.0 ) THEN
          cld_depth = 0.0
        END IF
        past_conv_ht(shi(i)) = cld_depth
        !          write(6,*) ' shallow cld depth 2 ',cld_depth
      END DO
    ELSE      ! shallow turbulence scheme so cloud ntpar
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        cld_depth = z_rho_sh(i,ntpar_sh(i)+1) - z_rho_sh(i,ntml_sh(i)+1)
        past_conv_ht(shi(i)) = cld_depth
      END DO
    END IF
  END IF

  ! moisture  - output depends on scheme

  IF (iconv_shallow == 1) THEN     ! G-R scheme

    ! G-R requires input of specific humidity

    IF (l_mr_physics) THEN ! requires conversion

      IF (l_q_interact) THEN    ! PC2

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_sh
            dqtt  = dqbydt_sh(i,k)+dqclbydt_sh(i,k)+dqcfbydt_sh(i,k)
            dqt   = dqtt*timestep
            denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
            dqbydt(shi(i),k)    = denom *                            &
                ( dqbydt_sh(i,k)*(1.0-qt(shi(i),k))+q_sh(i,k)*dqtt )
            dqclbydt(shi(i),k)  = denom *                            &
              ( dqclbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcl_sh(i,k)*dqtt )
            dqcfbydt(shi(i),k)  = denom *                            &
              ( dqcfbydt_sh(i,k)*(1.0-qt(shi(i),k))+qcf_sh(i,k)*dqtt )
            dq_dd(shi(i),k)     = denom *                            &
              ( dq_dd_temp(i,k)* (1.0-qt(shi(i),k))+q_sh(i,k)*dqtt )
          END DO
        END DO

      ELSE                      ! Not PC2
                                ! No qcl and qcf increments anyway
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_sh
            dqt = dqbydt_sh(i,k)*timestep
            denom = 1.0/((1.0-qt(shi(i),k))*(1.0-qt(shi(i),k)-dqt))
            dqbydt(shi(i),k)    = dqbydt_sh(i,k)* denom
            dqclbydt(shi(i),k)  = 0.0
            dqcfbydt(shi(i),k)  = 0.0
            dq_dd(shi(i),k)          = dq_dd_temp(i,k)* denom
          END DO
        END DO
      END IF                    ! End test on PC2

    ELSE   ! output is specific humidity therefore no problems

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          dqbydt(shi(i),k)   = dqbydt_sh(i,k)
          dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
          dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
          dq_dd(shi(i),k)          = dq_dd_temp(i,k)
        END DO
      END DO

    END IF     ! Test on l_mr_physics

  ELSE         ! turbulence scheme requires mixing ratio

    ! turbulence scheme outputs mixing ratios

    IF (l_mr_physics) THEN

      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          dqbydt(shi(i),k)   = dqbydt_sh(i,k)
          dqclbydt(shi(i),k) = dqclbydt_sh(i,k)
          dqcfbydt(shi(i),k) = dqcfbydt_sh(i,k)
          dq_dd(shi(i),k)    = dq_dd_temp(i,k)
        END DO
      END DO

    ELSE    ! Requires conversion increment to mixing ratio not q
      IF (l_q_interact) THEN     ! PC2
        ! At present turbulence scheme has no dqcl and dqcf so dqt=dq
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_sh
            dqtt = dqbydt_sh(i,k)
            dqt = dqtt*timestep
            denom = 1.0/((1.0+qt(shi(i),k))*(1.0+qt(shi(i),k)+dqt))
            dqbydt(shi(i),k)   = denom *                           &
                (dqbydt_sh(i,k)*(1.0+qt(shi(i),k))-q_sh(i,k)*dqtt)
            dqclbydt(shi(i),k) = denom *                           &
              (dqclbydt_sh(i,k)*(1.0+qt(shi(i),k))-qcl_sh(i,k)*dqtt)
            dqcfbydt(shi(i),k) = denom *                           &
              (dqcfbydt_sh(i,k)*(1.0+qt(shi(i),k))-qcf_sh(i,k)*dqtt)
            dq_dd(shi(i),k)    = denom *                           &
                (dq_dd_temp(i,k)*(1.0+qt(shi(i),k))-q_sh(i,k)*dqtt)
          END DO
        END DO

      ELSE                       ! Not PC2
                                 !  No qcl and qcf increments
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_sh
            dqt   = dqbydt_sh(i,k)*timestep
            denom = 1.0/((1.0+qt(shi(i),k))*(1.0+qt(shi(i),k)+dqt))
            dqbydt(shi(i),k)   =  dqbydt_sh(i,k)*denom
            dqclbydt(shi(i),k) = 0.0
            dqcfbydt(shi(i),k) = 0.0
            dq_dd(shi(i),k)    = dq_dd_temp(i,k)*denom
          END DO
        END DO
      END IF                     ! Test on PC2

    END IF       ! Test on l_mr_physics

  END IF        ! test on scheme

  DEALLOCATE(area_ud_temp)
  DEALLOCATE(area_dd_temp)
  DEALLOCATE(dv_dd_temp)
  DEALLOCATE(du_dd_temp)
  DEALLOCATE(dq_dd_temp)
  DEALLOCATE(dt_dd_temp)
  DEALLOCATE(rho_dry_theta_temp)
  DEALLOCATE(rho_dry_temp)

  ! diagnostics output

  IF (flg_up_flx) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        up_flux(shi(i),k)        = up_flux_sh(i,k)
      END DO
    END DO
  END IF
  ! Only values for mass flux option
  IF (flg_up_flx_half .AND. iconv_shallow == 1) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        up_flux_half(shi(i),k)        = up_flux_half_sh(i,k)
      END DO
    END DO
  END IF
  ! Shallow scheme only allowed to have downdraughts in new scheme
  IF (flg_dwn_flx .AND. l_new_dd) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        dwn_flux(shi(i),k)       = dwn_flux_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_uw_shall) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        uw_shall(shi(i),k)       = uw_shall_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_vw_shall) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        vw_shall(shi(i),k)       = vw_shall_sh(i,k)
      END DO
    END DO
  END IF

  IF (flg_wqt_flux) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        wqt_flux_sh(shi(i),k)       = wqt_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_wql_flux) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        wql_flux_sh(shi(i),k)       = wql_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_wthetal_flux) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        wthetal_flux_sh(shi(i),k)    = wthetal_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_wthetav_flux) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        wthetav_flux_sh(shi(i),k)    = wthetav_sh(i,k)
      END DO
    END DO
  END IF

  IF (flg_mf_shall) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        mf_shall(shi(i),k)        = up_flux_sh(i,k)
      END DO
    END DO
  END IF
  IF (iconv_shallow == 1) THEN
    IF (flg_entr_up) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          entrain_up(shi(i),k)     = entrain_up_sh(i,k)
        END DO
      END DO
    END IF
    IF (flg_detr_up) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          detrain_up(shi(i),k)     = detrain_up_sh(i,k)
        END DO
      END DO
    END IF
    IF (flg_entr_dwn) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          entrain_dwn(shi(i),k)    = entrain_dwn_sh(i,k)
        END DO
      END DO
    END IF
    IF (flg_detr_dwn) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          detrain_dwn(shi(i),k)    = detrain_dwn_sh(i,k)
        END DO
      END DO
    END IF
  END IF  ! iconv_shallow equal to 1
  IF (flg_dt_shall) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        dt_shall(shi(i),k)        = dthbydt_sh(i,k)
      END DO
    END DO
  END IF
  IF (flg_dq_shall) THEN
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_sh
        dq_shall(shi(i),k)        = dqbydt_sh(i,k)
      END DO
    END DO
  END IF
  IF (l_mom) THEN
    IF (flg_du_shall) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          du_shall(shi(i),k)        = dubydt_sh(i,k)
        END DO
      END DO
    END IF
    IF (flg_dv_shall) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          dv_shall(shi(i),k)        = dvbydt_sh(i,k)
        END DO
      END DO
    END IF
  END IF

  IF (l_tracer) THEN
    DO ktra = 1,ntra
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_sh
          dtrabydt(shi(i),k,ktra)  = dtrabydt_sh(i,k,ktra)
        END DO
      END DO
    END DO
  END IF

  ! Merge sh 3d rain & snow profiles

  DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1,n_sh
      rain_3d(shi(i),k) = rain_3d_sh(i,k)
      snow_3d(shi(i),k) = snow_3d_sh(i,k)
    END DO
  END DO

END IF ! n_sh > 0

IF ( (model_type == mt_single_column) .AND. &
     (n_sh == 0) ) THEN

  ! SCM diagnostics if no shallow

  IF (iconv_shallow == 2 .OR. iconv_shallow == 3) THEN

    ! Zero output arrays
    TmpScm1d(:)   = 0.0
    TmpScm2d(:,:) = 0.0

    !  cloud base fluxes etc
    sname='mb_cb'
    lname='cloud base mass flux'
    units='m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='mb_new_cb'
    lname='revised cloud base mass flux'
    units='m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wstar_up'
    lname='wstar_up'
    units='m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wthetal_cb'
    lname='wthetal at cloud base'
    units='Km/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wqt_cb'
    lname='wqt at cloud base'
    units='Kg/kg m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wql_cb'
    lname='wql at cloud base'
    units='Kg/kg m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wql_inv'
    lname='wql at inversion'
    units='Kg/kg m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wqr_inv'
    lname='wqr at inversion'
    units='Kg/kg m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='dthetav_cb'
    lname='dthetav across cloud base'
    units='K'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='dtheta_cb'
    lname='dtheta across cloud base'
    units='K'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)
    sname='drv_cb'
    lname='drv across cloud base'
    units='kg/kg'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wthetav_m_cb'
    lname='wthetav on lower side of cloud base'
    units='Km/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wqt_inv'
    lname='wqt at inversion'
    units='kg/kg m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='wthetal_inv'
    lname='wthetal at inversion'
    units='kg/kg m/s'
    CALL scmoutput(TmpScm1d,sname,lname,units,       &
         t_inst,d_point,default_streams,'',routinename)

    sname='dwqldz'
    lname='dql/dz  on model levels'
    units='Kg/kg s'
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)

    sname='qlup'
    lname='qlup  on model levels'
    units='Kg/kg '
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)
    sname='wup'
    lname='wup  on model levels'
    units='m/s'
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)

    sname='wthetav'
    lname='wthetav  on model levels'
    units='K m/s '
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)

    sname='wql'
    lname='wql flux in convective cloud'
    units=' '
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)

    sname='dqsat_dt'
    lname='dqsat_dt'
    units='K m/s'
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)
    sname='wthetal'
    lname='wthetal'
    units='K m/s'
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)
    sname='wqt'
    lname='wqt'
    units='Kg/kg m/s'
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)
    sname='wqr'
    lname='wqr on uv levels'
    units='Kg/kgm/s'
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)
    sname='precip_prod'
    lname='precip production on theta levels'
    units='/s'
    CALL scmoutput(TmpScm2d,sname,lname,units,       &
         t_inst,d_wet,default_streams,'',routinename)

  END IF ! type of shallow convection
END IF ! model_type





!-----------------------------------------------------------------------
! 3.0 CONGESTUS Convection   - call depending on switches
!-----------------------------------------------------------------------

IF (iconv_congestus >  0) THEN

  IF (n_cg  >   0) THEN   ! there must be congestus points

    ALLOCATE(rho_dry_temp(n_cg,nlev))
    ALLOCATE(rho_dry_theta_temp(n_cg,nlev))

    !-----------------------------------------------------------------------
    ! 3.1 Compress input variable arrays for congestus convection scheme to
    !     length n_cg (congestus points only)
    !-----------------------------------------------------------------------

    j = 0
    DO i = 1,npnts
      IF (cumulus_bl(i) .AND. l_congestus(i)) THEN
        j                        = j+1
        cgi(j)                   = i
      END IF
    END DO

    ! In only variables

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j = 1,n_cg
      bland_cg(j)           = bland(cgi(j))
      delthvu_cg(j)         = delthvu_bl(cgi(j))
      ntml_cg(j)            = ntml(cgi(j))
      ntpar_cg(j)           = ntpar(cgi(j))
      pstar_cg(j)           = pstar(cgi(j))
      recip_pstar_cg(j)     = recip_pstar(cgi(j))
      q1_sd_cg(j)           = q1_sd(cgi(j))
      t1_sd_cg(j)           = t1_sd(cgi(j))
      uw0_cg(j)             = uw0(cgi(j))
      vw0_cg(j)             = vw0(cgi(j))
      wstar_cg(j)           = wstar(cgi(j))
      wthvs_cg(j)           = wthvs(cgi(j))
      zlcl_uv_cg(j)         = zlcl_uv(cgi(j))
      ztop_uv_cg(j)         = ztop_uv(cgi(j))
      freeze_lev_cg(j)      = freeze_lev(cgi(j))
      entrain_coef_cg(j)    = entrain_coef(cgi(j))
      kterm_cg(j) = 0     ! initialise array Will only get set
                          ! if call deep scheme

    END DO

    DO k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO j = 1,n_cg
        p_layer_centres_cg(j,k)    = p_layer_centres(cgi(j),k)
        p_layer_boundaries_cg(j,k) = p_layer_boundaries(cgi(j),k)
        exner_layer_centres_cg(j,k)= exner_layer_centres(cgi(j),k)
        exner_layer_boundaries_cg(j,k)                            &
                                = exner_layer_boundaries(cgi(j),k)
        r_theta_cg(j,k)            = r_theta(cgi(j),k)
      END DO
    END DO

    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO j = 1,n_cg
        u_cg(j,k)           = u(cgi(j),k)
        v_cg(j,k)           = v(cgi(j),k)
        th_cg(j,k)          = th(cgi(j),k)
        z_theta_cg(j,k)     = z_theta(cgi(j),k)
        z_rho_cg(j,k)       = z_rho(cgi(j),k)
        r_rho_cg(j,k)       = r_rho(cgi(j),k)
        rho_cg(j,k)         = rho(cgi(j),k)
        rho_th_cg(j,k)      = rho_theta(cgi(j),k)
        rho_dry_theta_temp(j,k) = rho_dry_theta(cgi(j),k)
        rho_dry_temp(j,k)       = rho_dry(cgi(j),k)
        r2rho_cg(j,k)       = r2rho(cgi(j),k)
        r2rho_th_cg(j,k)    = r2rho_th(cgi(j),k)
        dr_across_rh_cg(j,k)       = dr_across_rh(cgi(j),k)
        dr_across_th_cg(j,k)       = dr_across_th(cgi(j),k)
      END DO
    END DO

    ! moisture  - input depends on scheme

    IF (iconv_congestus == 1) THEN     ! G-R scheme

      ! G-R requires input of specific humidity

      IF (l_mr_physics) THEN

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_cg
            q_cg(j,k)   = q(cgi(j),k)  /(1.0+qt(cgi(j),k))
            qse_cg(j,k) = qse(cgi(j),k)  ! copy specific value
            qcl_cg(j,k) = qcl(cgi(j),k)/(1.0+qt(cgi(j),k))
            qcf_cg(j,k) = qcf(cgi(j),k)/(1.0+qt(cgi(j),k))
          END DO
        END DO

      ELSE   ! input is specific humidity therefore no problems

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_cg
            q_cg(j,k)   = q(cgi(j),k)
            qse_cg(j,k) = qse(cgi(j),k)  ! copy specific value
            qcl_cg(j,k) = qcl(cgi(j),k)
            qcf_cg(j,k) = qcf(cgi(j),k)
          END DO
        END DO

      END IF     ! l_mr_physics

    ELSE         ! turbulence scheme requires mixing ratio

      ! turbulence scheme requires mixing ratios as input

      IF (l_mr_physics) THEN

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_cg
            q_cg(j,k)   = q(cgi(j),k)
            qse_cg(j,k) = qse_mix(cgi(j),k)  ! copy mixing ratio value
            qcl_cg(j,k) = qcl(cgi(j),k)
            qcf_cg(j,k) = qcf(cgi(j),k)
          END DO
        END DO

      ELSE       ! Require conversion

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_cg

            q_cg(j,k)     = q(cgi(j),k)/(1.0 -qt(cgi(j),k))

            ! check for negative values
            IF (q_cg(j,k) <  0.0) THEN
              IF (printstatus >= prstatus_normal) THEN
                WRITE(umMessage,'(a17,g26.18,a5,2i6)') 'problem q_mix -ve ',       &
                          q_cg(j,k),' j,k ',j,k
                CALL umPrint(umMessage,src='glue_conv-gconv5a')
              END IF
              q_cg(j,k) = 0.0
            END IF

            qse_cg(j,k) = qse_mix(cgi(j),k) ! copy mixing ratio value

            ! No checks on qcl and qcf as not expected to be used or altered.

            qcl_cg(j,k) = qcl(cgi(j),k)/(1.0-qt(cgi(j),k))
            qcf_cg(j,k) = qcf(cgi(j),k)/(1.0-qt(cgi(j),k))
          END DO
        END DO

      END IF      ! Test on l_mr_physics

    END IF        ! Test on scheme (iconv_congestus)


    ! In/out variables

    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO j = 1,n_cg
        cf_liquid_cg(j,k)   = cf_liquid(cgi(j),k)
        cf_frozen_cg(j,k)   = cf_frozen(cgi(j),k)
        bulk_cf_cg(j,k)     = bulk_cf(cgi(j),k)
      END DO
    END DO

    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_cg
            tracer_cg(j,k,ktra)  = tracer(cgi(j),k,ktra)
          END DO
        END DO
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_cg
            dtrabydt_cg(j,k,ktra)  = 0.0
          END DO
        END DO
      END DO
    END IF

    ALLOCATE(dt_dd_temp(n_cg,nlev))
    ALLOCATE(dq_dd_temp(n_cg,nlev))
    ALLOCATE(du_dd_temp(n_cg,nlev))
    ALLOCATE(dv_dd_temp(n_cg,nlev))
    ALLOCATE(area_ud_temp(n_cg,nlev))
    ALLOCATE(area_dd_temp(n_cg,nlev))

    !-----------------------------------------------------------------------
    ! 3.2 Call congestus convection code
    !-----------------------------------------------------------------------

    IF (iconv_congestus == 1) THEN   ! G-R type scheme

      ! new mass flux routine (at present a copy of shallow scheme)

      CALL congest_conv(nbl,nlev,ntra,n_cca_lev,n_cg,trlev        &
,                       bland_cg                                  &
,                       delthvu_cg,exner_layer_centres_cg         &
,                       exner_layer_boundaries_cg                 &
,                       l_calc_dxek, l_q_interact                 &
,                       l_tracer, ntml_cg, ntpar_cg               &
,                       pstar_cg,p_layer_centres_cg               &
,                       p_layer_boundaries_cg,z_theta_cg,z_rho_cg &
,                       r_theta_cg,r_rho_cg, rho_th_cg, rho_cg    &
,                       rho_dry_theta_temp, rho_dry_temp          &
,                       r2rho_th_cg,r2rho_cg, dr_across_th_cg     &
,                       dr_across_rh_cg                           &
,                       q_cg,q1_sd_cg,t1_sd_cg,th_cg,timestep     &
,                       u_cg,v_cg,uw0_cg,vw0_cg,wstar_cg,wthvs_cg &
,                       entrain_coef_cg                           &
,                       zlcl_uv_cg,ztop_uv_cg,freeze_lev_cg       &
,                       recip_pstar_cg,qse_cg                     &

                      ! InOut
,                     bulk_cf_cg,cf_frozen_cg,cf_liquid_cg,qcf_cg &
,                     qcl_cg,tracer_cg,w2p_cg                     &

                      ! Out
,                     cape_out_cg,cclwp_cg,ccw_cg,cca_cg          &
,                     dbcfbydt_cg,dcffbydt_cg,dcflbydt_cg         &
,                     dqbydt_cg,dqcfbydt_cg,dqclbydt_cg,dthbydt_cg&
,                     dubydt_cg,dvbydt_cg,dtrabydt_cg             &
,                     detrain_up_cg,detrain_dwn_cg                &
,                     entrain_up_cg,entrain_dwn_cg                &
,                     iccb_cg,icct_cg,lcca_cg,lcbase_cg,lctop_cg  &
,                     rain_cg,snow_cg,rain_3d_cg,snow_3d_cg       &
,                     up_flux_cg, up_flux_half_cg                 &
,                     dwn_flux_cg,uw_shall_cg,vw_shall_cg,kterm_cg&
,                     tcw_cg,cca_2d_cg                            &
,                     rbuoy_p_out_cg,the_out_cg,thp_out_cg        &
,                     qe_out_cg,qp_out_cg, dt_dd_temp,dq_dd_temp  &
,                     du_dd_temp,dv_dd_temp                       &
,                     area_ud_temp, area_dd_temp)


      !         Else If (iconv_congestus == 2) then  ! turbulence scheme

      !            Call TURB_CONGESTUS_CONV(

      !        Not available yet - still to be written

    END IF   ! test of type of congestus scheme

    !-----------------------------------------------------------------------
    ! 3.3 Write data from congestus convection points to full arrays
    !-----------------------------------------------------------------------
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i = 1,n_cg
      cape_out(cgi(i))           = cape_out_cg(i)
      cclwp(cgi(i))              = cclwp_cg(i)
      iccb(cgi(i))               = iccb_cg(i)
      icct(cgi(i))               = icct_cg(i)

      lcca(cgi(i))               = lcca_cg(i)
      lcbase(cgi(i))             = lcbase_cg(i)
      lctop(cgi(i))              = lctop_cg(i)

      rain(cgi(i))               = rain_cg(i)
      snow(cgi(i))               = snow_cg(i)
      precip_cong(cgi(i))        = rain_cg(i) + snow_cg(i)
      wstar_dn(cgi(i))           = wstar_cg(i)  ! wstar_dn
      ! may be required if congestus not forced to stop at ntpar
      kterm_congest(cgi(i))      = kterm_cg(i)
    END DO

    ! required at present to replicate old code
    IF (.NOT. l_ccrad) THEN
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_cg
        tcw(cgi(i))              = tcw_cg(i)
        cca_2d(cgi(i))           = cca_2d_cg(i)
      END DO
    END IF

    IF (l_mom) THEN
      DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_cg
          dubydt(cgi(i),k)         = dubydt_cg(i,k)
          dvbydt(cgi(i),k)         = dvbydt_cg(i,k)
        END DO
      END DO
    END IF

    IF (flg_w_eqn) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_cg
          w2p(cgi(i),k) = w2p_cg(i,k)
        END DO
      END DO
    END IF

    IF (flg_area_ud) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_cg
          area_ud(cgi(i),k) = area_ud_temp(i,k)
        END DO
      END DO
    END IF

    IF (flg_area_dd) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_cg
          area_dd(cgi(i),k) = area_dd_temp(i,k)
        END DO
      END DO
    END IF

    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_cg
        dthbydt(cgi(i),k)        = dthbydt_cg(i,k)
        dcflbydt(cgi(i),k)       = dcflbydt_cg(i,k)
        dcffbydt(cgi(i),k)       = dcffbydt_cg(i,k)
        dbcfbydt(cgi(i),k)       = dbcfbydt_cg(i,k)
        ccw(cgi(i),k)            = ccw_cg(i,k)
        dt_dd(cgi(i),k)          = dt_dd_temp(i,k)
        du_dd(cgi(i),k)          = du_dd_temp(i,k)
        dv_dd(cgi(i),k)          = dv_dd_temp(i,k)
      END DO
    END DO

    DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_cg
        cca(cgi(i),k) = cca_cg(i,k)
      END DO
    END DO

    IF (l_ccrad) THEN

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_cg
        iccb0  (cgi(i)) = iccb  (cgi(i))
        icct0  (cgi(i)) = icct  (cgi(i))
        cclwp0 (cgi(i)) = cclwp (cgi(i))
        lcbase0(cgi(i)) = lcbase(cgi(i))
      END DO

      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_cg
          ccw0(cgi(i),k) = ccw_cg(i,k) * ccw_sh_knob
        END DO
      END DO

      DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_cg
          cca0(cgi(i),k) = cca_cg(i,k) * cca_sh_knob
        END DO
      END DO

    ELSE ! Non-ccrad code

      IF ( l_q_interact ) THEN
        ! Do nothing as congest does not contribute to PC2
      ELSE

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_cg
          iccb0  (cgi(i)) = iccb  (cgi(i))
          icct0  (cgi(i)) = icct  (cgi(i))
          cclwp0 (cgi(i)) = cclwp (cgi(i))
          lcbase0(cgi(i)) = lcbase(cgi(i))
        END DO

        DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i=1, n_cg
            cca0(cgi(i),k) = cca_cg(i,k)
          END DO
        END DO

      END IF

    END IF ! l_ccrad

    IF (model_type == mt_single_column) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_cg
          rbuoy_p_out(cgi(i),k) = rbuoy_p_out_cg(i,k)
          the_out(cgi(i),k)     = the_out_cg(i,k)
          thp_out(cgi(i),k)     = thp_out_cg(i,k)
          qe_out(cgi(i),k)      = qe_out_cg(i,k)
          qp_out(cgi(i),k)      = qp_out_cg(i,k)
        END DO
      END DO
    END IF ! model_type


    ! moisture  - output depends on scheme

    IF (iconv_congestus == 1) THEN     ! G-R scheme

      ! Outputs specific humidity increments

      IF (l_mr_physics) THEN ! requires conversion

        IF (l_q_interact) THEN    ! PC2

          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_cg
              dqtt  = dqbydt_cg(i,k)+dqclbydt_cg(i,k)+dqcfbydt_cg(i,k)
              dqt   = dqtt*timestep
              denom = 1.0/((1.0-qt(cgi(i),k))*(1.0-qt(cgi(i),k)-dqt))
              dqbydt(cgi(i),k)   = denom *                              &
                 ( dqbydt_cg(i,k)*(1.0-qt(cgi(i),k))+q_cg(i,k)*dqtt )
              dqclbydt(cgi(i),k) = denom *                              &
               ( dqclbydt_cg(i,k)*(1.0-qt(cgi(i),k))+qcl_cg(i,k)*dqtt )
              dqcfbydt(cgi(i),k) = denom *                              &
               ( dqcfbydt_cg(i,k)*(1.0-qt(cgi(i),k))+qcf_cg(i,k)*dqtt )
              dq_dd(cgi(i),k)    = denom *                              &
                 ( dq_dd_temp(i,k)*(1.0-qt(cgi(i),k))+q_cg(i,k)*dqtt )
            END DO
          END DO

        ELSE                      ! Not PC2
                                ! No qcl and qcf increments anyway
          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_cg
              dqt = dqbydt_cg(i,k)*timestep
              denom = 1.0/((1.0-qt(cgi(i),k))*(1.0-qt(cgi(i),k)-dqt))
              dqbydt(cgi(i),k)    =  dqbydt_cg(i,k)*denom
              dqclbydt(cgi(i),k)  = 0.0
              dqcfbydt(cgi(i),k)  = 0.0
              dq_dd(cgi(i),k)    = dq_dd_temp(i,k)*denom
            END DO
          END DO
        END IF                     ! End test on PC2

      ELSE   ! output is specific humidity therefore no problems

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_cg
            dqbydt(cgi(i),k)   = dqbydt_cg(i,k)
            dqclbydt(cgi(i),k) = dqclbydt_cg(i,k)
            dqcfbydt(cgi(i),k) = dqcfbydt_cg(i,k)
            dq_dd(cgi(i),k)    = dq_dd_temp(i,k)
          END DO
        END DO

      END IF     ! l_mr_physics

    ELSE         ! turbulence scheme

      ! turbulence scheme outputs mixing ratios

      IF (l_mr_physics) THEN

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_cg
            dqbydt(cgi(i),k)   = dqbydt_cg(i,k)
            dqclbydt(cgi(i),k) = dqclbydt_cg(i,k)
            dqcfbydt(cgi(i),k) = dqcfbydt_cg(i,k)
            dq_dd(cgi(i),k)    = dq_dd_temp(i,k)
          END DO
        END DO

      ELSE    ! Requires conversion increment to mixing ratio not q

        IF (l_q_interact) THEN       ! PC2

          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_cg
              dqtt  = dqbydt_cg(i,k)+dqclbydt_cg(i,k)+dqcfbydt_cg(i,k)
              dqt   = dqtt*timestep
              denom = 1.0/((1.0+qt(cgi(i),k))*(1.0+qt(cgi(i),k)+dqt))
              dqbydt(cgi(i),k)   = denom *                           &
                   (dqbydt_cg(i,k)*(1.0+qt(cgi(i),k))-q_cg(i,k)*dqtt)
              dqclbydt(cgi(i),k) = denom *                           &
                 (dqclbydt_cg(i,k)*(1.0+qt(cgi(i),k))-qcl_cg(i,k)*dqtt)
              dqcfbydt(cgi(i),k) = denom *                           &
                 (dqcfbydt_cg(i,k)*(1.0+qt(cgi(i),k))-qcf_cg(i,k)*dqtt)
              dq_dd(cgi(i),k)    = denom *                              &
                 ( dq_dd_temp(i,k)*(1.0-qt(cgi(i),k))+q_cg(i,k)*dqtt )
            END DO
          END DO

        ELSE                         ! Not PC2
                                    ! No qcl and qcf increments
          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_cg
              dqt   = dqbydt_cg(i,k)*timestep
              denom = 1.0/((1.0+qt(cgi(i),k))*(1.0+qt(cgi(i),k)+dqt))
              dqbydt(cgi(i),k)   = dqbydt_cg(i,k)*denom
              dqclbydt(cgi(i),k) = 0.0
              dqcfbydt(cgi(i),k) = 0.0
              dq_dd(cgi(i),k)    = dq_dd_temp(i,k)*denom
            END DO
          END DO
        END IF                       ! End test on PC2

      END IF     ! Test on l_mr_physics

    END IF        ! Test on scheme (iconv_congestus)

    DEALLOCATE(area_dd_temp)
    DEALLOCATE(area_ud_temp)
    DEALLOCATE(dv_dd_temp)
    DEALLOCATE(du_dd_temp)
    DEALLOCATE(dq_dd_temp)
    DEALLOCATE(dt_dd_temp)
    DEALLOCATE(rho_dry_theta_temp)
    DEALLOCATE(rho_dry_temp)

    ! Merge dp 3d rain & snow profiles

    DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_cg
        rain_3d(cgi(i),k) = rain_3d_cg(i,k)
        snow_3d(cgi(i),k) = snow_3d_cg(i,k)
      END DO
    END DO

    IF (flg_up_flx) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_cg
          up_flux(cgi(i),k)        = up_flux_cg(i,k)
        END DO
      END DO
    END IF
    IF (flg_up_flx_half) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_cg
          up_flux_half(cgi(i),k)        = up_flux_half_cg(i,k)
        END DO
      END DO
    END IF

    IF (flg_mf_congest) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_cg
          mf_congest(cgi(i),k)        = up_flux_cg(i,k)
        END DO
      END DO
    END IF
    IF (flg_dt_congest) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_cg
          dt_congest(cgi(i),k)        = dthbydt_cg(i,k)
        END DO
      END DO
    END IF
    IF (flg_dq_congest) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_cg
          dq_congest(cgi(i),k)        = dqbydt_cg(i,k)
        END DO
      END DO
    END IF
    IF (l_mom) THEN
      IF (flg_du_congest) THEN
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_cg
            du_congest(cgi(i),k)        = dubydt_cg(i,k)
          END DO
        END DO
      END IF
      IF (flg_dv_congest) THEN
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_cg
            dv_congest(cgi(i),k)        = dvbydt_cg(i,k)
          END DO
        END DO
      END IF
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_cg
            dtrabydt(cgi(i),k,ktra)  = dtrabydt_cg(i,k,ktra)
          END DO
        END DO
      END DO
    END IF

  END IF    ! number of congestus points > 0

  ! set past convective depth

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  DO i = 1,n_cg
    cld_depth = z_rho_cg(i,kterm_cg(i)+1) - z_rho_cg(i,ntml_cg(i)+1)
    past_conv_ht(cgi(i)) = cld_depth
  END DO

END IF      ! iconv_congestus > 0


IF (model_type == mt_single_column) THEN
  !-----------------------------------------------------------------------
    ! make copies of rbuoy etc. for checking purposes
    ! Used after mid convection to see if mid has updated values
  DO k=1, nlev
    DO i=1, npnts
      rbuoy_p_out_copy(i,k) = rbuoy_p_out(i,k)
      the_out_copy(i,k)     = the_out(i,k)
      thp_out_copy(i,k)     = thp_out(i,k)
      qe_out_copy(i,k)      = qe_out(i,k)
      qp_out_copy(i,k)      = qp_out(i,k)
    END DO
  END DO
END IF ! model_type


!-----------------------------------------------------------------------
! 4.0 MID-LEVEL Convection
! 4.1 Set lowest level that mid level convection can trigger
!-----------------------------------------------------------------------

IF (iconv_mid >  0) THEN
  IF (n_md > 0) THEN       ! check there are mid points

    ALLOCATE(rho_dry_temp(n_md,nlev))
    ALLOCATE(rho_dry_theta_temp(n_md,nlev))

    ! Create index of points where mid-level is possible

    j = 0
    DO i = 1,npnts
      IF (l_mid(i)) THEN
        j      = j+1
        mdi(j) = i
      END IF
    END DO

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i = 1,n_md
      IF (.NOT. cumulus_bl(mdi(i))) THEN
        midtrig(i) = ntml(mdi(i)) + 1
        IF (ntml(mdi(i))  ==  nbl) THEN
          midtrig(i) = ntml(mdi(i))
        END IF
      ELSE

        !  Cumulus points ie deep or shallow convection has occurred

        midtrig(i) = ntpar(mdi(i)) + 2


        ! NTPAR has a maximum value which can be less than the top level for
        ! deep. Deep convection may terminate above this. Need to prevent mid
        ! level convection occurring at the same levels as deep.

        IF ( .NOT. l_shallow_bl(mdi(i)) .AND.                           &
            .NOT. l_congestus(mdi(i)) ) THEN
          IF (kterm_deep(mdi(i)) >  ntpar(mdi(i))+2) THEN
            midtrig(i) = kterm_deep(mdi(i))
          END IF
        END IF  ! deep points

        ! if congestus calls deep scheme can terminate above ntpar
        IF (l_congestus(mdi(i))) THEN
          IF (kterm_congest(mdi(i)) >  ntpar(mdi(i))+2) THEN
            midtrig(i) = kterm_congest(mdi(i))
          END IF
        END IF

      END IF
    END DO

    !-----------------------------------------------------------------------
    ! 4.2 Copy all input arrays to arrays ending in _md for passing to
    !     mid-level scheme
    !-----------------------------------------------------------------------

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO j=1,n_md
      bland_md(j)           = bland(mdi(j))
      ntml_md(j)            = ntml(mdi(j))
      ntpar_md(j)           = ntpar(mdi(j))
      pstar_md(j)           = pstar(mdi(j))
      recip_pstar_md(j)     = recip_pstar(mdi(j))
      q1_sd_md(j)           = q1_sd(mdi(j))
      t1_sd_md(j)           = t1_sd(mdi(j))
      w_max_md(j)           = w_max(mdi(j))
      freeze_lev_md(j)      = freeze_lev(mdi(j))
    END DO


    DO k = 0,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_md
        p_layer_centres_md(i,k)        = p_layer_centres(mdi(i),k)
        p_layer_boundaries_md(i,k)     = p_layer_boundaries(mdi(i),k)
        exner_layer_centres_md(i,k)    = exner_layer_centres(mdi(i),k)
        exner_layer_boundaries_md(i,k) = exner_layer_boundaries(mdi(i),k)
        r_theta_md(i,k)                = r_theta(mdi(i),k)
      END DO
    END DO
    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_md
        u_md(i,k)           = u(mdi(i),k)
        v_md(i,k)           = v(mdi(i),k)
        w_md(i,k)           = w(mdi(i),k)
        th_md(i,k)          = th(mdi(i),k)
        cf_liquid_md(i,k)   = cf_liquid(mdi(i),k)
        cf_frozen_md(i,k)   = cf_frozen(mdi(i),k)
        bulk_cf_md(i,k)     = bulk_cf(mdi(i),k)
        z_theta_md(i,k)     = z_theta(mdi(i),k)
        z_rho_md(i,k)       = z_rho(mdi(i),k)
        r_rho_md(i,k)       = r_rho(mdi(i),k)
        r2rho_th_md(i,k)     = r2rho_th(mdi(i),k)
        r2rho_md(i,k)        = r2rho(mdi(i),k)
        rho_theta_md(i,k)    = rho_theta(mdi(i),k)
        rho_md(i,k)          = rho(mdi(i),k)
        rho_dry_theta_temp(i,k) = rho_dry_theta(mdi(i),k)
        rho_dry_temp(i,k)       = rho_dry(mdi(i),k)
        dr_across_th_md(i,k) = dr_across_th(mdi(i),k)
        dr_across_rh_md(i,k) = dr_across_rh(mdi(i),k)
      END DO
    END DO

    IF (.NOT. l_ccrad) THEN
      ! cca_2d is generated from all cloud on a given gridbox,
      ! (i.e. Shallow/deep/congestus + mid-level) instead of CCRad
      ! which generates cloud for each cloud type.
      ! So current cca_2d array is copied to cca_2d_md
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_md
        tcw_md(i)    = tcw(mdi(i))
        cca_2d_md(i) = cca_2d(mdi(i))
      END DO
    END IF

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i=1, n_md
      iccb_md(i) = iccb(mdi(i))
      icct_md(i) = icct(mdi(i))
    END DO


    ! moisture  - input depends on scheme


    IF (iconv_mid == 1) THEN     ! G-R scheme

      ! G-R requires input of specific humidity

      IF (l_mr_physics) THEN

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_md
            q_md(j,k)   = q(mdi(j),k)  /(1.0+qt(mdi(j),k))
            qse_md(j,k) = qse(mdi(j),k)  ! copy specific value
            qcl_md(j,k) = qcl(mdi(j),k)/(1.0+qt(mdi(j),k))
            qcf_md(j,k) = qcf(mdi(j),k)/(1.0+qt(mdi(j),k))
          END DO
        END DO

      ELSE   ! input is specific humidity therefore no problems

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_md
            q_md(j,k)   = q(mdi(j),k)
            qse_md(j,k) = qse(mdi(j),k) ! copy specific value
            qcl_md(j,k) = qcl(mdi(j),k)
            qcf_md(j,k) = qcf(mdi(j),k)
          END DO
        END DO

      END IF      ! Test on l_mr_physics

    ELSE         ! future turbulence scheme requires mixing ratio

      ! turbulence scheme requires mixing ratios as input

      IF (l_mr_physics) THEN

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_md
            q_md(j,k)   = q(mdi(j),k)
            qse_md(j,k) = qse_mix(mdi(j),k) ! copy mixing ratio value
            qcl_md(j,k) = qcl(mdi(j),k)
            qcf_md(j,k) = qcf(mdi(j),k)
          END DO
        END DO

      ELSE       ! Require conversion

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO j = 1,n_md

            q_md(j,k)     = q(mdi(j),k)/(1.0 -qt(mdi(j),k))

            ! check for negative values
            IF (q_md(j,k) <  0.0) THEN
              IF (printstatus >= prstatus_normal) THEN
                WRITE(umMessage,'(a18,g26.18,a5,2i6)') 'problem q_mix -ve ',   &
                     q_md(j,k),' j,k ',j,k
                CALL umPrint(umMessage,src='glue_conv-gconv5a')
              END IF
              q_md(j,k) = 0.0
            END IF
            qse_md(j,k) = qse_mix(mdi(j),k) ! copy mixing ratio value
            ! Not sure whether qcl and qcf will be used or altered and don't
            ! know whether additional tests will be required if converting.

            qcl_md(j,k) = qcl(mdi(j),k)/(1.0-qt(mdi(j),k))
            qcf_md(j,k) = qcf(mdi(j),k)/(1.0-qt(mdi(j),k))
          END DO
        END DO

      END IF     ! Test on l_mr_physics

    END IF        ! test on scheme (iconv_mid)

    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO k = 1,trlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_md
            tracer_md(i,k,ktra)  = tracer(mdi(i),k,ktra)
          END DO
        END DO
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_md
            dtrabydt_md(i,k,ktra)  = 0.0
          END DO
        END DO
      END DO
    END IF

    ALLOCATE(dt_dd_temp(n_md,nlev))
    ALLOCATE(dq_dd_temp(n_md,nlev))
    ALLOCATE(du_dd_temp(n_md,nlev))
    ALLOCATE(dv_dd_temp(n_md,nlev))
    ALLOCATE(area_ud_temp(n_md,nlev))
    ALLOCATE(area_dd_temp(n_md,nlev))

    !-----------------------------------------------------------------------
    ! 4.3 Call mid-level convection code
    !-----------------------------------------------------------------------

    IF (iconv_mid == 1) THEN ! Gregory-Rowntree mid level convection

      CALL mid_conv_5a(                                           &
                !IN
                nbl,nlev,ntra,n_cca_lev,n_md,trlev,               &
                bland_md,w_max_md,exner_layer_centres_md,         &
                exner_layer_boundaries_md, l_calc_dxek,           &
                l_q_interact, l_tracer, midtrig, ntml_md,         &
                ntpar_md, freeze_lev_md,                          &
                pstar_md,p_layer_centres_md,                      &
                p_layer_boundaries_md,r_theta_md,r_rho_md,        &
                z_theta_md,z_rho_md,rho_md,rho_theta_md,          &
                rho_dry_theta_temp, rho_dry_temp,                 &
                r2rho_th_md,r2rho_md,dr_across_th_md,             &
                dr_across_rh_md,                                  &
                q_md,q1_sd_md,t1_sd_md,th_md,timestep,            &
                u_md,v_md,w_md,recip_pstar_md,qse_md,             &
                !INOUT
                bulk_cf_md,cf_frozen_md,cf_liquid_md,             &
                qcf_md,qcl_md,tracer_md,w2p_md,                   &
                !OUT
                cape_out_md,cclwp_md,ccw_md,cca_md,               &
                dbcfbydt_md,dcffbydt_md,dcflbydt_md,              &
                dqbydt_md,dqcfbydt_md,dqclbydt_md,                &
                dthbydt_md,dubydt_md,dvbydt_md,dtrabydt_md,       &
                detrain_up_md,detrain_dwn_md,entrain_up_md,       &
                entrain_dwn_md,                                   &
                iccb_md,icct_md,lcca_md,lcbase_md,lctop_md,       &
                rain_md,snow_md,rain_3d_md,snow_3d_md,            &
                up_flux_md,up_flux_half_md,                       &
                dwn_flux_md,tcw_md,l_mid_md,cca_2d_md,            &
                rbuoy_p_out_md,the_out_md,thp_out_md,             &
                qe_out_md,qp_out_md,                              &
                uw_mid_md,vw_mid_md, cfl_limited_md,              &
                dt_dd_temp, dq_dd_temp, du_dd_temp, dv_dd_temp,   &
                area_ud_temp, area_dd_temp,                       &
                error_point)

      IF (error_point /= 0) THEN
        errorstatus = 3   ! will cause model to fail
        WRITE(cmessage,'(a47,i12,a16,i2)')                                  &
        'Mid conv went to the top of the model at point ',                  &
        mdi(error_point),' in seg on call ',call_number

        CALL ereport(routinename, errorstatus, cmessage )

      END IF

    ELSE        ! turbulence based alternative

      errorstatus = 1   ! will cause model to fail
      cmessage='Alternative mid-level scheme not available yet'

      CALL ereport(routinename, errorstatus, cmessage )

      !          call mid_turb_conv( ) ?
    END IF

    !-----------------------------------------------------------------------
    ! 4.4 Write data from mid-level convection to full arrays
    !-----------------------------------------------------------------------

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i = 1,n_md

      ! Cloud variables - only overwrite deep or shallow values if
      ! iccb_md  and icct_md > 0

      IF (iccb_md(i)  >   0 .AND. icct_md(i)  >   0) THEN
        iccb(mdi(i))     = iccb_md(i)
        icct(mdi(i))     = icct_md(i)
      END IF

      ! Overwrite lowest cloud values only if cumulus = .F.

      IF (.NOT. cumulus_bl(mdi(i))) THEN
        lcca(mdi(i))     = lcca_md(i)
        lcbase(mdi(i))   = lcbase_md(i)
        lctop(mdi(i))    = lctop_md(i)
      END IF
    END DO

    ! Write remaining data to full arrays

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    DO i = 1,n_md
      cape_out(mdi(i))   = cape_out(mdi(i)) + cape_out_md(i)
      cclwp(mdi(i))      = cclwp(mdi(i))    + cclwp_md(i)
      rain(mdi(i))       = rain(mdi(i))     + rain_md(i)
      snow(mdi(i))       = snow(mdi(i))     + snow_md(i)
      precip_mid(mdi(i)) = rain_md(i)       + snow_md(i)
      mid_cfl_limited(mdi(i)) = cfl_limited_md(i)
      l_mid_all(mdi(i))  = l_mid_md(i)
      !===================================================================
      ! NOTE: At this point if l_ccrad = T, then cca_2d_md is ONLY equal
      !       to that from mid-level cloud on a given grid point.  The
      !       original code I.E. l_ccrad = F means that cca_2d_md will
      !       include that from sh/dp aswell.
      !===================================================================

    END DO

    ! Merge md 3d rain & snow profiles

    DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1,n_md
        rain_3d(mdi(i),k)     = rain_3d(mdi(i),k) + rain_3d_md(i,k)
        snow_3d(mdi(i),k)     = snow_3d(mdi(i),k) + snow_3d_md(i,k)
      END DO
    END DO
    IF (l_mom) THEN
      IF (flg_uw_mid) THEN
        DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i=1,n_md
            uw_mid(mdi(i),k)     = uw_mid(mdi(i),k) + uw_mid_md(i,k)
          END DO
        END DO
      END IF
      IF (flg_vw_mid) THEN
        DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i=1,n_md
            vw_mid(mdi(i),k)     = vw_mid(mdi(i),k) + vw_mid_md(i,k)
          END DO
        END DO
      END IF

      DO k=1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          dubydt(mdi(i),k)    = dubydt(mdi(i),k) + dubydt_md(i,k)
          dvbydt(mdi(i),k)    = dvbydt(mdi(i),k) + dvbydt_md(i,k)
        END DO
      END DO
    END IF

    DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i = 1,n_md
        dthbydt(mdi(i),k)     = dthbydt(mdi(i),k)  + dthbydt_md(i,k)
        dcflbydt(mdi(i),k)    = dcflbydt(mdi(i),k) + dcflbydt_md(i,k)
        dcffbydt(mdi(i),k)    = dcffbydt(mdi(i),k) + dcffbydt_md(i,k)
        dbcfbydt(mdi(i),k)    = dbcfbydt(mdi(i),k) + dbcfbydt_md(i,k)
        ccw(mdi(i),k)         = ccw(mdi(i),k)      + ccw_md(i,k)
        dt_dd(mdi(i),k)       = dt_dd(mdi(i),k)    + dt_dd_temp(i,k)
        du_dd(mdi(i),k)       = du_dd(mdi(i),k)    + du_dd_temp(i,k)
        dv_dd(mdi(i),k)       = dv_dd(mdi(i),k)    + dv_dd_temp(i,k)
      END DO
    END DO


    DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_md
        cca(mdi(i),k) = cca(mdi(i),k) + cca_md(i,k)
      END DO
    END DO

    IF (flg_w_eqn) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          w2p(mdi(i),k) = w2p(mdi(i),k) + w2p_md(i,k)
        END DO
      END DO
    END IF

    IF (flg_area_ud) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          area_ud(mdi(i),k) = area_ud(mdi(i),k) + area_ud_temp(i,k)
        END DO
      END DO
    END IF

    IF (flg_area_dd) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          area_dd(mdi(i),k) = area_dd(mdi(i),k) + area_dd_temp(i,k)
        END DO
      END DO
    END IF

    IF (l_ccrad) THEN

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_md
        cclwp0 (mdi(i)) = cclwp (mdi(i))
        iccb0  (mdi(i)) = iccb  (mdi(i))
        icct0  (mdi(i)) = icct  (mdi(i))
        lcbase0(mdi(i)) = lcbase(mdi(i))
      END DO


      ! Then scaling is done by CCRad, even if PC2 = .true.
      ! To zero CCA (section 0) for PC2 use the CCRAD knobs
      ! cca_knobs should be set to 0.0 in gui/namelist for PC2

      ! Anvils if selected will affect both section 0/5 diags
      ! if ccrads are not set to 0.0, this may required further
      ! attention in future
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          ccw0(mdi(i),k) = ccw0(mdi(i),k) + ccw_md(i,k) * ccw_md_knob
        END DO
      END DO

      DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          cca0(mdi(i),k) = cca0(mdi(i),k) + cca_md(i,k) * cca_md_knob
        END DO
      END DO

      ! If l_ccrad=T add on cca_2d_sh,cca_2d_dp to cca_2d.
      ! Contributions from cca_2d_sh/cca_2d_dp have already
      ! been applied to cca_sh/cca_dp

      ! However, contributions from cca_2d_sh/cca_2d_dp need
      ! to be included in the cca_2d diagnostic so as not to
      ! upset any Downstream products that use it.
      ! (at this cca_2d point only holds cca_2d_md)

      ! cca_2d_md entered mid_conv as an empty array,
      ! i.e. no contribution from shallow/congestus/deep.
      ! Add any contribution from mid level to preserve
      ! CCA_2d diagnostic.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_md
        cca_2d(mdi(i)) = cca_2d(mdi(i)) + cca_2d_md(i)
      END DO

    ELSE ! Non-CCRad

      IF (l_q_interact) THEN
        ! PC2 Code (PC2 has no mid-level contribution)
      ELSE
        ! Non-PC2 Code
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          cclwp0 (mdi(i)) = cclwp (mdi(i))
          iccb0  (mdi(i)) = iccb  (mdi(i))
          icct0  (mdi(i)) = icct  (mdi(i))
          lcbase0(mdi(i)) = lcbase(mdi(i))
        END DO

        DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i=1, n_md
            cca0(mdi(i),k) = cca0(mdi(i),k) + cca_md(i,k)
          END DO
        END DO
      END IF

      ! In non-ccrad configs, the cca_2d_md array is not empty when
      ! entering mid-conv; it hold contributions from deep_conv/
      ! congest_conv/shallow_conv, so we can just copy cca_2d_md
      !to cca_2d
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
      DO i=1, n_md
        cca_2d(mdi(i)) = cca_2d_md(i)
      END DO

    END IF ! l_ccrad

    IF (l_cca_md_prog) THEN
      DO k=1, n_cca_lev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          cca0_md(mdi(i),k) = cca_md(i,k)
        END DO
      END DO
    END IF


    IF (model_type == mt_single_column) THEN
      DO k=1, nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i=1, n_md
          ! with mid level arrays from parcel, check that values
          ! have changed from the values they were initialised to
          ! so that we only overwrite at points where mid level convection
          ! has taken place and don't inadvertently over-write shallow and
          ! deep values with the original model profiles.
          !
          ! If test should be on parcel values, 'cos environment values
          ! shouldn't change
          ! looking for points where shallow and deep didn't affect
          ! original profile, but mid did
          IF ((rbuoy_p_out_copy(mdi(i),k)  ==  0.0)                     &
               .AND. (rbuoy_p_out_md(i,k)  /=  0.0)) THEN
            rbuoy_p_out(mdi(i),k) = rbuoy_p_out_md(i,k)
          END IF

          IF ((thp_out_copy(mdi(i),k)  ==  th_md(i,k)) .AND.            &
               (thp_out_md(i,k)  /=  th_md(i,k))) THEN
            the_out(mdi(i),k) = the_out_md(i,k)
            thp_out(mdi(i),k) = thp_out_md(i,k)
          END IF
          IF ((qp_out_copy(i,k)  ==  q_md(i,k)) .AND.                   &
               (qp_out_md(i,k)  /=  q_md(i,k))) THEN
            qe_out(mdi(i),k) = qe_out_md(i,k)
            qp_out(mdi(i),k) = qp_out_md(i,k)
          END IF

        END DO
      END DO
    END IF ! model_type


    ! moisture  - output depends on scheme

    IF (iconv_mid == 1) THEN     ! G-R scheme

      ! Output in specific humidity

      IF (l_mr_physics) THEN ! requires conversion

        IF (l_q_interact) THEN    ! PC2

          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_md
              dqtt  = dqbydt_md(i,k)+dqclbydt_md(i,k)+dqcfbydt_md(i,k)
              dqt   = dqtt*timestep
              denom = 1.0/((1.0-qt(mdi(i),k))*(1.0-qt(mdi(i),k)-dqt))
              dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + denom *          &
                       ( dqbydt_md(i,k)*(1.0-qt(mdi(i),k))+               &
                       q_md(i,k)*dqtt )
              dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + denom *          &
                       ( dqclbydt_md(i,k)*(1.0-qt(mdi(i),k))+             &
                       qcl(mdi(i),k)*dqtt )
              dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + denom *          &
                       ( dqcfbydt_md(i,k)*(1.0-qt(mdi(i),k))+             &
                       qcf(mdi(i),k)*dqtt )
              dq_dd(mdi(i),k)      = dq_dd(mdi(i),k) + denom *           &
                      (dq_dd_temp(i,k)*(1.0+qt(mdi(i),k))- q_md(i,k)*dqtt)
            END DO
          END DO

        ELSE                       ! Not PC2
                                   ! No qcl and qcf increments anyway
          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_md
              dqt = dqbydt_md(i,k)*timestep
              denom = 1.0/((1.0-qt(mdi(i),k))*(1.0-qt(mdi(i),k)-dqt))
              dqbydt(mdi(i),k)    =  dqbydt(mdi(i),k) + dqbydt_md(i,k)*denom
              dqclbydt(mdi(i),k)  = 0.0
              dqcfbydt(mdi(i),k)  = 0.0
              dq_dd(mdi(i),k)      = dq_dd(mdi(i),k) + dq_dd_temp(i,k)*denom
            END DO
          END DO

        END IF                     ! End test on PC2

      ELSE   ! output is specific humidity therefore no problems

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_md
            dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + dqbydt_md(i,k)
            dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + dqclbydt_md(i,k)
            dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + dqcfbydt_md(i,k)
            dq_dd(mdi(i),k)    = dq_dd(mdi(i),k)    + dq_dd_temp(i,k)
          END DO
        END DO

      END IF     ! Test on scheme (iconv_mid)

    ELSE         ! turbulence scheme requires mixing ratio

      ! turbulence scheme outputs mixing ratios

      IF (l_mr_physics) THEN

        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_md
            dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + dqbydt_md(i,k)
            dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + dqclbydt_md(i,k)
            dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + dqcfbydt_md(i,k)
            dq_dd(mdi(i),k)    = dq_dd(mdi(i),k)    + dq_dd_temp(i,k)
          END DO
        END DO

      ELSE    ! Requires conversion increment to mixing ratio not q

        IF (l_q_interact) THEN    ! PC2

          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_md
              dqtt  = dqbydt_md(i,k)+dqclbydt_md(i,k)+dqcfbydt_md(i,k)
              dqt   = dqtt*timestep
              denom = 1.0/((1.0+qt(mdi(i),k))*(1.0+qt(mdi(i),k)+dqt))
              dqbydt(mdi(i),k)   = dqbydt(mdi(i),k)   + denom *      &
                     (dqbydt_md(i,k)*(1.0+qt(mdi(i),k))-             &
                     q_md(i,k)*dqtt)
              dqclbydt(mdi(i),k) = dqclbydt(mdi(i),k) + denom *      &
                     (dqclbydt_md(i,k)*(1.0+qt(mdi(i),k))-           &
                     qcl(mdi(i),k)*dqtt)
              dqcfbydt(mdi(i),k) = dqcfbydt(mdi(i),k) + denom *      &
                     (dqcfbydt_md(i,k)*(1.0+qt(mdi(i),k))-           &
                     qcf(i,k)*dqtt)
              dq_dd(mdi(i),k)      = dq_dd(mdi(i),k) + denom *      &
                     (dq_dd_temp(i,k)*(1.0+qt(mdi(i),k))- q_md(i,k)*dqtt)
            END DO
          END DO

        ELSE                      ! Not PC2
                                   ! No qcl and qcf increments
          DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
            DO i = 1,n_md
              dqt   = dqbydt_md(i,k)*timestep
              denom = 1.0/((1.0+qt(i,k))*(1.0+qt(mdi(i),k)+dqt))
              dqbydt(mdi(i),k)   = dqbydt(mdi(i),k) + dqbydt_md(i,k)*denom
              dqclbydt(mdi(i),k) = 0.0
              dqcfbydt(mdi(i),k) = 0.0
              dq_dd(mdi(i),k)      = dq_dd(mdi(i),k) + dq_dd_temp(i,k)*denom
            END DO
          END DO

        END IF                    ! End test on PC2

      END IF  ! Test on l_mr_physics

    END IF        ! test on scheme (iconv_mid)

    DEALLOCATE(area_dd_temp)
    DEALLOCATE(area_ud_temp)
    DEALLOCATE(dv_dd_temp)
    DEALLOCATE(du_dd_temp)
    DEALLOCATE(dq_dd_temp)
    DEALLOCATE(dt_dd_temp)
    DEALLOCATE(rho_dry_theta_temp)
    DEALLOCATE(rho_dry_temp)

    IF (flg_up_flx) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          up_flux(mdi(i),k)      = up_flux(mdi(i),k) + up_flux_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_up_flx_half) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          up_flux_half(mdi(i),k) = up_flux_half(mdi(i),k)+up_flux_half_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_dwn_flx) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          dwn_flux(mdi(i),k)      = dwn_flux(mdi(i),k) + dwn_flux_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_entr_up) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          entrain_up(mdi(i),k)   = entrain_up(mdi(i),k) + entrain_up_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_detr_up) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          detrain_up(mdi(i),k)   = detrain_up(mdi(i),k) + detrain_up_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_entr_dwn) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          entrain_dwn(mdi(i),k)  = entrain_dwn(mdi(i),k) + entrain_dwn_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_detr_dwn) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          detrain_dwn(mdi(i),k)  = detrain_dwn(mdi(i),k) + detrain_dwn_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_mf_midlev) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          mf_midlev(mdi(i),k)    = up_flux_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_dt_midlev) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          dt_midlev(mdi(i),k)    = dthbydt_md(i,k)
        END DO
      END DO
    END IF
    IF (flg_dq_midlev) THEN
      DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
        DO i = 1,n_md
          dq_midlev(mdi(i),k)    = dqbydt_md(i,k)
        END DO
      END DO
    END IF
    IF (l_mom) THEN
      IF (flg_du_midlev) THEN
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_md
            du_midlev(mdi(i),k)   = dubydt_md(i,k)
          END DO
        END DO
      END IF
      IF (flg_dv_midlev) THEN
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_md
            dv_midlev(mdi(i),k)   = dvbydt_md(i,k)
          END DO
        END DO
      END IF
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO k = 1,nlev
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          DO i = 1,n_md
            dtrabydt(mdi(i),k,ktra) = dtrabydt(mdi(i),k,ktra)        &
                                    + dtrabydt_md(i,k,ktra)
          END DO
        END DO
      END DO
    END IF

  END IF     ! test n_md >0
END IF     ! test on iconv_mid

!-----------------------------------------------------------------------
! 5.0 Work out which points convection has occurred at.
!-----------------------------------------------------------------------

nconv_all=0
DO i = 1,npnts
  IF (cumulus_bl(i) .OR. l_mid_all(i)) THEN
    nconv_all = nconv_all + 1
    index1(nconv_all) = i
  END IF
END DO

!-----------------------------------------------------------------------
! 6.0  update tracer field
!      More efficient to do here rather than in subroutines.
!      This seems to be an expensive bit of code on NEC.
!      Changed to operate on just convective points.
!-----------------------------------------------------------------------

IF (l_tracer .AND. (nconv_all >  0)) THEN


  ! Adjust timestep to prevent any negative values invading the tracer
  ! fields (adjusted timestep is a function of geographical  location and
  ! tracer type.

  DO ktra=1,ntra
    ! initialise step to timestep
    DO i = 1,nconv_all
      limited_step(i) =timestep
    END DO

    DO k=1,nlev
      DO j = 1,nconv_all
        i=index1(j)
        ! negative increments  may have a problem
        IF (dtrabydt(i,k,ktra)  <   0.0 ) THEN

          step_test(j) = (0.9999 * ABS(tracer(i,k,ktra))) /       &
                    (ABS(dtrabydt(i,k,ktra)) + safety_margin)

          IF (step_test(j)   <   limited_step(j) ) THEN
            ! then increment is bigger than tracer and timestep needs to be reduced
            limited_step (j) = step_test(j)
          END IF

        END IF
      END DO
    END DO

    ! Update tracer field using limited_step.

    DO k = 1,nlev
      ! NEC compiler directive
      !CDIR NODEP
      DO j = 1,nconv_all
        i=index1(j)
        tracer(i,k,ktra) = tracer(i,k,ktra) +                       &
                         dtrabydt(i,k,ktra) * limited_step(j)
      END DO
    END DO

  END DO  ! ktra loop
END IF   ! L_tracer

! Ensure no convective cloud in top layer

DO i=1,npnts
  cca(i,n_cca_lev) = 0.0
  cca0(i,n_cca_lev) = 0.0
END DO


! Cleanup allocation
IF (iconv_deep == 2 .OR. (iconv_shallow >= 2) ) THEN
  ! Using a turbulence scheme so a mixing ratio version
  DEALLOCATE( qse_mix )
END IF


!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE glue_conv_5a

END MODULE glue_conv_5a_mod
