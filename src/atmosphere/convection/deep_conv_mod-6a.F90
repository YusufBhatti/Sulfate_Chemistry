! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Deep convection scheme

MODULE deep_conv_6a_mod

IMPLICIT NONE

!
! Description:
!   Deep convection scheme
!   works on points diagnosed as deep in subroutine CONV_DIAG.
!   Called by GLUE_CONV.
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DEEP_CONV_6A_MOD'

CONTAINS

SUBROUTINE deep_conv_6a(nbl,nlev,ntra,n_cca_lev,n_dp,trlev,       &
                       bland, delthvu,                            &
                       exner_rho,                                 &
                       exner_layer_centres,                       &
                       exner_layer_boundaries,                    &
                       l_q_interact,                              &
                       l_tracer,ntml,ntpar,                       &
                       pstar,p_layer_centres,                     &
                       p_layer_boundaries,                        &
                       z_theta, z_rho,                            &
                       r_theta, r_rho,                            &
                       rho_theta, rho,                            &
                       rho_dry_theta, rho_dry,                    &
                       r2rho_th, r2rho,                           &
                       dr_across_th, dr_across_rh,                &
                       conv_prog_precip,                          &
                       q,th,                                      &
                       timestep,                                  &
                       u,v,w,uw0,vw0,w_max,wstar,qsat_lcl,        &
                       entrain_coef,zlcl_uv,                      &
                       freeze_lev,recip_pstar,qse,                &
                       l_scm_convss_dg,                           &

                       ! InOut
                       bulk_cf,cf_frozen,cf_liquid,qcf,           &
                       qcl,tracer,w2p,                            &
                       scm_convss_dg,                             &

                       ! Out
                       cape_out,cclwp,ccw,cca,                    &
                       dbcfbydt,dcffbydt,dcflbydt,dqbydt,dqcfbydt,&
                       dqclbydt,dthbydt,                          &
                       dubydt,dvbydt,dtrabydt,                    &
                       detrain_up,detrain_dwn,                    &
                       entrain_up,entrain_dwn,                    &
                       iccb,icct,lcca,                            &
                       lcbase,lctop,rain,snow,                    &
                       rain_3d, snow_3d, up_flux, up_flux_half,   &
                       dwn_flux,uw_deep,vw_deep,kterm,tcw,cca_2d, &
                       ind_cape_reduced,cape_ts_used,cfl_limited, &
                       ind_deep,dt_dd,dq_dd, du_dd, dv_dd,        &
                       area_ud, area_dd,                          &
                       error_point)

USE planet_constants_mod, ONLY:                                   &
    r, cp, kappa, pref, repsilon, c_virtual, g

USE cv_run_mod, ONLY:                                             &
    l_mom, l_eman_dd, l_new_dd, cldbase_opt_dp, cape_min,         &
    w_cape_limit, cape_timescale, deep_cmt_opt, bl_cnv_mix,       &
    cca2d_dp_opt, cca_dp_knob, ccw_dp_knob, limit_pert_opt,       &
    icvdiag, l_anvil, cnv_wat_load_opt, l_cv_conserve_check,      &
    l_3d_cca, l_cmt_heating, l_snow_rain, l_use_dd

USE cv_param_mod, ONLY:                                           &
    total_condensed_water, srf_precip,                            &
    a_land, a_sea, b_land, b_sea,                                 &
    thpixs_deep, qpixs_deep, c_mass, wcape_fac,                   &
    max_dp_thpert, min_dp_thpert, max_dp_qpert_fac,               &
    a_cape, b_cape, max_cape, a_cb, b_cb, c_cb

USE cv_dependent_switch_mod, ONLY:                                &
    dp_on, mdet_dp_on, dp_ent_on, dp_new_termc

USE cv_stash_flg_mod, ONLY:                                       &
    flg_up_flx,  flg_up_flx_half, flg_dwn_flx,                    &
    flg_entr_up, flg_detr_up, flg_detr_up, flg_detr_dwn,          &
    flg_entr_dwn, flg_uw_dp, flg_vw_dp, flg_mf_deep, flg_w_eqn,   &
    flg_area_ud

USE scm_convss_dg_mod, ONLY: scm_convss_dg_type

USE water_constants_mod, ONLY: lc, lf, tm

USE conversions_mod, ONLY: rsec_per_day

USE yomhook,     ONLY: lhook, dr_hook
USE parkind1,    ONLY: jprb, jpim
USE ereport_mod, ONLY: ereport
USE umPrintMgr,  ONLY: umPrint, umMessage, newline,               &
                       printstatus, prstatus_normal
USE lift_par_6a_mod
USE convec2_6a_mod
USE water_loading_mod
USE cor_engy_6a_mod
USE mix_ipert_6a_mod
USE cmt_heating_mod
USE layer_cn_6a_mod, ONLY: layer_cn_6a, deep
USE eman_dd_rev_mod, ONLY: eman_dd_rev
USE dd_evap_ctl_mod, ONLY: dd_evap_ctl

USE calc_3d_cca_mod, ONLY: calc_3d_cca
USE cmt_mass_mod, ONLY: cmt_mass
USE dd_all_call_6a_mod, ONLY: dd_all_call_6a
USE deep_cmt_incr_mod, ONLY: deep_cmt_incr
USE deep_grad_stress_mod, ONLY: deep_grad_stress
USE deep_ngrad_stress_mod, ONLY: deep_ngrad_stress
USE deep_turb_cmt_mod, ONLY: deep_turb_cmt
USE eman_dd_mod, ONLY: eman_dd
USE evap_bcb_nodd_all_mod, ONLY: evap_bcb_nodd_all
USE flag_wet_mod, ONLY: flag_wet

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    l_new_qsat_conv !Currently defaults to FALSE

IMPLICIT NONE

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

! Arguments with intent IN:


INTEGER, INTENT(IN) :: nbl      ! No. of boundary layer levels

INTEGER, INTENT(IN) :: nlev     ! No. of model layers

INTEGER, INTENT(IN) :: ntra     ! No. of tracer fields

INTEGER, INTENT(IN) :: n_cca_lev! No. of convective cloud
                                ! amount levels (1 for 2D,
                                               ! nlevs for 3D)

INTEGER, INTENT(IN) :: n_dp     ! No. of deep convection points

INTEGER, INTENT(IN) :: trlev    ! No. of model levels on which
                                ! tracers are included

LOGICAL, INTENT(IN) :: bland(n_dp) ! Land/sea mask

REAL, INTENT(IN)    :: delthvu(n_dp) ! a measure of CAPE used to cal wcld

REAL, INTENT(IN)    :: exner_rho(n_dp,nlev)  ! Exner on rho levels

REAL, INTENT(IN)    :: exner_layer_centres(n_dp,0:nlev) !Exner

REAL, INTENT(IN)    :: exner_layer_boundaries(n_dp,0:nlev)
                                ! Exner at half level above
                                ! exner_layer_centres

LOGICAL, INTENT(IN) :: l_q_interact ! Switch allows overwriting
                                    ! parcel variables when
                                    ! calculating condensate incr.

LOGICAL, INTENT(IN) :: l_tracer ! Switch for inclusion of tracers

INTEGER, INTENT(IN) :: ntml(n_dp) ! Top level of surface mixed
                                  ! layer defined relative to
                                  ! theta,q grid

INTEGER, INTENT(IN) :: ntpar(n_dp) ! Top level of initial parcel
                                   ! ascent in BL scheme defined
                                   ! relative to theta,q grid

REAL, INTENT(IN)    :: pstar(n_dp) ! Surface pressure (Pa)

REAL, INTENT(IN)    :: p_layer_centres(n_dp,0:nlev) ! Pressure (Pa)

REAL, INTENT(IN)    :: p_layer_boundaries(n_dp,0:nlev) ! Pressure
                                                       ! at half level above
                                                       ! p_layer_centres (Pa)

REAL, INTENT(IN) :: z_theta(n_dp,nlev)      ! height of theta levels (m)
REAL, INTENT(IN) :: z_rho(n_dp,nlev)        ! height of rho levels (m)
REAL, INTENT(IN) :: r_theta(n_dp,0:nlev)    ! radius of theta levels (m)
REAL, INTENT(IN) :: r_rho(n_dp,nlev)        ! radius of rho levels (m)
REAL, INTENT(IN) :: rho_theta(n_dp,nlev)    ! density for theta lev (kg/m3)
REAL, INTENT(IN) :: rho(n_dp,nlev)          ! density for rho lev (kg/m3)
REAL, INTENT(IN) ::         &
  rho_dry_theta(n_dp,nlev)  & ! dry density for theta lev (kg/m3)
 ,rho_dry(n_dp,nlev)          ! dry density for rho lev (kg/m3)

REAL, INTENT(IN) :: r2rho_th(n_dp,nlev)     ! radius**2 density for
                                            ! theta lev (kg/m)
REAL, INTENT(IN) :: r2rho(n_dp,nlev)        ! radius**2 density for
                                            ! rho lev (kg/m)
REAL, INTENT(IN) :: dr_across_th(n_dp,nlev) ! thickness of theta levels (m)
REAL, INTENT(IN) :: dr_across_rh(n_dp,nlev) ! thickness of rho levels (m)
REAL, INTENT(IN) :: conv_prog_precip(n_dp,nlev) ! Surface precipitation based
                                                ! 3d convective prognostic in
                                                ! kg/m2/s
REAL, INTENT(IN)    :: q(n_dp,nlev) ! Model mixing ratio (kg/kg)

REAL, INTENT(IN)    :: th(n_dp,nlev) !Model potential temperature (K)

REAL, INTENT(IN)    :: timestep ! Model timestep (s)

REAL, INTENT(IN)    :: &
  u(n_dp,nlev)         & ! Model u field (m/s)
 ,v(n_dp,nlev)         & ! Model v field (m/s)
 ,w(n_dp,nlev)           ! Model w field (m/s)

REAL, INTENT(IN)    :: uw0(n_dp) ! U-comp of surface stress (N/m2)

REAL, INTENT(IN)    :: vw0(n_dp) ! V-comp of surface stress (N/m2)

REAL, INTENT(IN)    :: wstar(n_dp) ! Convective velocity scale (m/s)

REAL, INTENT(IN)    :: w_max(n_dp) ! max w in column
                                   !for use in scale dependent cape timescale

REAL, INTENT(IN)    :: entrain_coef(n_dp) ! entrainment coefficient

REAL, INTENT(IN)    :: qsat_lcl(n_dp) ! qsat at cloud base (kg/kg)

REAL, INTENT(IN)    :: zlcl_uv(n_dp) !Lifting condensation level
                                     ! defined for the uv grid (m)

INTEGER, INTENT(IN) :: freeze_lev(n_dp) ! Level index for freezing level

REAL, INTENT(IN) :: recip_pstar(n_dp) ! Reciprocal of pstar array

REAL, INTENT(IN) :: qse(n_dp,nlev) ! Saturation mixing ratio of
                                   ! cloud environment (kg/kg)

! Arguments with intent INOUT:


REAL, INTENT(INOUT) :: bulk_cf(n_dp,nlev) ! Bulk total cloud volume ( )

REAL, INTENT(INOUT) :: cf_frozen(n_dp,nlev) ! Frozen water cloud volume ( )

REAL, INTENT(INOUT) :: cf_liquid(n_dp,nlev) ! Liq water cloud volume ( )

REAL, INTENT(INOUT) :: qcf(n_dp,nlev) ! Ice condensate mix ratio (kg/kg)

REAL, INTENT(INOUT) :: qcl(n_dp,nlev) ! Liq condensate mix ratio (kg/kg)

REAL, INTENT(INOUT) :: tracer(n_dp,trlev,ntra) !Model tracer fields (kg/kg)

REAL, INTENT(INOUT) :: w2p(n_dp,nlev) ! (Parcel vertical velocity)^2 [(m/s)^2]

! Structure containing SCM convection sub-step diagnostics
! (needs intent inout as contains allocatable arrays that need to
! retain their allocated status on input as well as output)
TYPE(scm_convss_dg_type), INTENT(INOUT) :: scm_convss_dg( n_dp )
! Flag for SCM convection sub-step diagnostics
LOGICAL, INTENT(IN) :: l_scm_convss_dg


! Arguments with intent OUT:

REAL, INTENT(OUT) :: cape_out(n_dp) ! Saved convective available
                                    ! potential energy for diagnostic
                                    ! output (J/kg)

REAL, INTENT(OUT) :: cclwp(n_dp)    ! Condensed water path (kg/m^2)

REAL, INTENT(OUT) :: ccw(n_dp,nlev)
                                ! Convective cloud liquid water
                                ! on model levels (kg/kg)

REAL, INTENT(OUT) :: cca(n_dp,n_cca_lev)
                                ! Convective cloud amount
                                ! on model levels (0-1)

REAL, INTENT(OUT) ::   &
  dbcfbydt(n_dp,nlev)  & ! Increments to total cld volume due to convection(/s)
 ,dcffbydt(n_dp,nlev)  & ! Increments to ice cloud volume due to convection (/s)
 ,dcflbydt(n_dp,nlev)  & ! Increments to liq cloud volume due to convection (/s)
 ,dqbydt(n_dp,nlev)    & ! Increments to q due to convection (kg/kg/s)
 ,dqcfbydt(n_dp,nlev)  & ! Increments to ice
                         ! condensate due to convection(kg/kg/s)
 ,dqclbydt(n_dp,nlev)  & ! Increments to liq condensate due to convection
                         ! (kg/kg/s)
 ,dthbydt(n_dp,nlev)   & ! Increments to potential temp. due to convection (K/s)
 ,dubydt(n_dp,nlev)    & ! Increments to U due to CMT (m/s2)
 ,dvbydt(n_dp,nlev)      ! Increments to V due to CMT (m/s2)

REAL, INTENT(OUT) ::       &
  dtrabydt(n_dp,nlev,ntra)   !Increment to tracer due to convection (kg/kg/s)

REAL, INTENT(OUT) ::     &
  detrain_up(n_dp,nlev)  & ! Fractional detrainment rate into updraughts
                           ! (Pa/s)
 ,detrain_dwn(n_dp,nlev) & ! Fractional detrainment rate into downdraughts
                           ! (Pa/s)
 ,entrain_up(n_dp,nlev)  & ! Fractional entrainment rate into updraughts
                           ! (Pa/s)
 ,entrain_dwn(n_dp,nlev)   ! Fractional entrainment rate into downdraughts
                           ! (Pa/s)

INTEGER, INTENT(OUT) :: iccb(n_dp) ! Convective cloud base level

INTEGER, INTENT(OUT) :: icct(n_dp) ! Convective cloud top level

REAL, INTENT(OUT) :: lcca(n_dp) ! Lowest conv. cloud amt. (%)

INTEGER, INTENT(OUT) :: lcbase(n_dp) ! Lowest conv. cloud base level

INTEGER, INTENT(OUT) :: lctop(n_dp) ! Lowest conv. cloud top level

REAL, INTENT(OUT) :: rain(n_dp) ! Surface convective rainfall (kg/m2/s)

REAL, INTENT(OUT) :: snow(n_dp) ! Surface convective snowfall (kg/m2/s)

REAL, INTENT(OUT) :: rain_3d(n_dp,nlev) ! Convective rainfall flux (kg/m2/s)

REAL, INTENT(OUT) :: snow_3d(n_dp,nlev) ! Convective snowfall flux (kg/m2/s)

REAL, INTENT(OUT) :: up_flux(n_dp,nlev) ! Updraught mass flux (Pa/s)

REAL, INTENT(OUT) :: up_flux_half(n_dp,nlev) ! Updraught mass flux on
                                             ! half levels (Pa/s)
REAL, INTENT(OUT) :: dwn_flux(n_dp,nlev) ! Downdraught mass flux (Pa/s)

REAL, INTENT(OUT) ::  &
  uw_deep(n_dp,nlev)  & ! X-comp. of stress from deep convection (kg/m/s2)
 ,vw_deep(n_dp,nlev)    ! Y-comp. of stress from deep convection (kg/m/s2)

INTEGER, INTENT(OUT) :: kterm(n_dp) ! Level at which deep
                                    ! convection terminates,
                                    ! required by mid level scheme
REAL, INTENT(OUT) :: tcw(n_dp)   ! Total condensed water(kg/m2/s)
                                 ! required by mid-level CCA cal.

REAL, INTENT(OUT) :: cca_2d(n_dp) ! 2D convective cloud amount (%)

REAL, INTENT(OUT) ::      &
  ind_cape_reduced(n_dp)  & ! 1.0 - if CAPE reduced applies to several
                            ! CAPE options
 ,cape_ts_used(n_dp)      & ! cape timescale used for deep convection (s)
 ,cfl_limited(n_dp)       & ! Indicator of CFL limited convection
 ,ind_deep(n_dp)          & ! 1.0 if real deep convection else 0.0
 ,dt_dd(n_dp,nlev)        & ! dT/dt from DD and evap below cloud base (K/s)
 ,dq_dd(n_dp,nlev)        & ! dq/dt from DD and evap below cloud base (kg/kg/s)
 ,du_dd(n_dp,nlev)        & ! du/dt from DD (m/s/s)
 ,dv_dd(n_dp,nlev)        & ! dv/dt from DD (m/s/s)
 ,area_ud(n_dp,nlev)      & ! fractional area of updraughts
 ,area_dd(n_dp,nlev)        ! fractional area of downdraughts

INTEGER, INTENT(OUT) :: error_point     ! 0 no error
                                        ! location of problem deep point

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

! Height above surface of model levels ...
REAL :: zkm1(n_dp)              ! ... k-1   [m]
REAL :: zk(n_dp)                ! ... k     [m]
REAL :: zkp12(n_dp)             ! ... k+1/2 [m]
REAL :: zkp1(n_dp)              ! ... k+1   [m]

INTEGER :: index1(n_dp),index2(n_dp)

INTEGER :: ncposs               ! No. of points which may convect

INTEGER :: nconv                ! No. of convecting points

REAL :: amdetk(n_dp)            ! Mixing detrainment coefficient
                                ! at level k multiplied by
                                ! appropriate layer thickness

REAL :: b_calc                  ! Coefficient in thpert calc.

REAL :: c_calc                  ! Coefficient in thpert calc.

REAL :: cape(n_dp)              ! Convective available potential
                                ! energy (J/kg)

REAL :: dq_sat_env              ! dqsat/dT  (kg/kg/K)

REAL ::                 &
  cca_2d_term(n_dp)     &  ! 2d CCA for termination level
 ,dcpbydt_term(n_dp)       ! Rate of change of cape (J/kg/s) at termination

REAL :: dcpbydt(n_dp)           ! Rate of change of cape (J/kg/s)

REAL :: depth(n_dp)             ! Depth of convective cloud (m)

REAL :: ekp14(n_dp)             ! Entrainment coefficients at
                                ! level k+1/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: ekp34(n_dp)             ! Entrainment coefficients at
                                ! level k+3/4 multiplied by
                                ! appropriate layer thickness
                                !(dimensionless)

REAL :: exk(n_dp)               ! Exner ratio at layer k

REAL :: exkp1(n_dp)             ! Exner ratio at layer k+1

REAL :: flxmax(n_dp)            ! Maximum initial convective
                                ! mass flux (Pa/s)

REAL :: flx_init(n_dp)          ! Initial mass flux at cloud base
                                ! (Pa/s)

REAL :: flx_init_new(n_dp)      ! flx_init scaled to destroy cape
                                ! over timescale cape_timescale (Pa/s)

REAL :: flxmax_init(n_dp)       ! Maximum possible initial mass
                                ! flux (limited to the mass in
                                ! the initial convecting layer
                                ! in Pa/s)

REAL :: max_cfl(n_dp)           ! Max cfl ratio over a convecting
                                ! layer

REAL :: p_lcl(n_dp)             ! Pressure at LCL (Pa)

REAL :: precip(n_dp,nlev)       ! Amount of precip from each layer
                                ! from each layer (kg/m/s)

REAL :: pk(n_dp)                ! Pressure at midpoint of layer
                                ! k (Pa)

REAL :: pkp1(n_dp)              ! Pressure at midpoint of layer
                                ! k+1 (Pa)

REAL :: delpk(n_dp)             ! Pressure difference over layer
                                ! k (Pa)

REAL :: delpkp1(n_dp)           ! Pressure difference over layer
                                ! k+1 (Pa)

REAL :: delpkp12(n_dp)          ! Pressure difference between
                                ! layers k and k+1 (Pa)

REAL :: delp_uv_k(n_dp)         ! Pressure difference across uv
                                ! layer k (Pa)

REAL :: delp_uv_kp1(n_dp)       ! Pressure difference across uv
                                ! layer k+1 (Pa)

REAL :: q_lcl(n_dp)             ! Mixing ratio at LCL (kg/kg)

REAL :: qse_lcl(n_dp)           ! Saturated q at LCL (kg/kg)

REAL :: rhum(n_dp)              ! Dummy relative humidity
                                ! (only used on shallow points)

REAL :: t_lcl(n_dp)             ! Temperature at LCL (K)

REAL :: th_lcl(n_dp)            ! Theta at LCL (K)

REAL :: thv_pert(n_dp)          ! Theta_v parcel perturbation (K)

REAL :: thpert(n_dp)            ! Theta parcel perturbation (K)

REAL :: qpert(n_dp)             ! q parcel perturbation (kg/kg)

REAL :: pstar_w_cape_limit(n_dp)! scaled critical vertical velocity

REAL :: wsc_o_mb(n_dp)          ! Dummy argument for layer_cn
                                ! Convective velocity scale divided
                                ! by cloud base mass flux mb

INTEGER :: start_lev(n_dp)      ! Convection initiation level
INTEGER :: start_lev_c2(n_dp)   ! Compressed convection initiation level

LOGICAL :: bgmk(n_dp)           ! Mask for points where parcel in
                                ! layer k is saturated
LOGICAL :: bgmk_term(n_dp)      ! Mask for points where parcel in
                                ! layer k is saturated at termination

LOGICAL :: blatent(n_dp)        ! Mask for points where latent heat has
                                ! been released

LOGICAL :: bwater(n_dp,2:nlev)  ! Mask for points at which
                                ! condensate is liquid

LOGICAL :: blowst(n_dp)         ! Dummy variable indicating low
                                ! enough stability for convection
                                ! to occur

LOGICAL :: bterm(n_dp)          ! Mask for points which have
                                ! stopped convecting

LOGICAL :: bconv(n_dp)          ! Mask for points at which
                                ! convection is occurring

LOGICAL :: bcposs(n_dp)         ! Mask for points passing
                                ! initial stability test


! Parcel variables


REAL :: qpi(n_dp)               ! Initial parcel mixing ratio
                                !(kg/kg)

REAL :: qp(n_dp,nlev)           ! Parcel mixing ratio (kg/kg)

REAL :: thpi(n_dp)              ! Initial parcel potential temp. (K)

REAL :: thp(n_dp,nlev)          ! Parcel potential temp (K)

REAL :: trap(n_dp,nlev,ntra)    ! Tracer content of parcel (kg/kg)

REAL :: expi(n_dp)              ! Initial parcel exner pressure

REAL :: flx(n_dp,nlev)          ! Parcel massflux (Pa/s)

REAL :: xsbmin_v(n_dp,nlev)     ! Minmum parcel buoyancy excess

REAL :: thpixs_v(n_dp,nlev)     ! Theta parcel excess (K)

REAL :: qpixs_v(n_dp,nlev)      ! Q parcel excess(kg/kg)

! PC2

REAL :: qclp(n_dp,nlev)         ! Parcel liquid condensated mixing
                                ! ratio in layer k (kg/kg)

REAL :: qcfp(n_dp,nlev)         ! Parcel frozen condensated mixing
                                ! ratio in layer k (kg/kg)


! Parameters
REAL, PARAMETER :: cfl_limit = 1.0 ! Max CFL ratio allowed
REAL, PARAMETER :: minflx = TINY(flx_init_new)  ! minimum allowable
                                                ! initial mass flux

! CMT variables  - those used depend on scheme

! Required by Gregory-Kershaw scheme operating in plume calculation

INTEGER ::         &
 nstart(n_dp)        ! Level for start of plume

REAL ::            &
 eflux_u_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                     !  top of layer (Pa m/s)
,eflux_v_ud(n_dp)  & ! Vertical eddy flux of momentum due to UD at
                     ! bottom of layer (Pa m/s)
,up(n_dp,nlev)     & ! Parcel U (m/s)
,vp(n_dp,nlev)     & ! Parcel V (m/s)
,zsurf(n_dp)         ! Height of start of plume = 0.1*zlcl

! Required by Turbulence base scheme called after plume calculation

INTEGER ::             &
 nlcl_uv(n_dp)         & ! Level index for LCL
,ntop_uv(n_dp)         & ! Level index for top of layer
,n_0degc(n_dp)         & ! Level index for zero degrees
,cu_term(n_dp)         & ! Indicies for CMT subs
,cu_tend(n_dp)           ! Indicies for CMT subs

REAL ::                &
 mass_dwn(nlev,n_dp)   & ! Downdraught mass flux (Pa/s)
,p_uv(nlev,n_dp)       & ! Pressure of model level (Pa)
,phalf_uv(nlev,n_dp)   & ! Pressure of half level (Pa)
,plcl_uv(n_dp)         & ! Pressure at LCL (Pa)
,ptop_uv(n_dp)         & ! Pressure at top of cloud layer (Pa)
,p_0degc_uv(n_dp)      & ! Pressure of zero degree level (Pa)
,rho_uv(nlev,n_dp)     & ! Density on uv level (kg/m3)
,visc(nlev,n_dp)       & ! CMT eddy viscosity (m2/s)
,uw(nlev,n_dp)         & ! U- comp stress profile (N/m2)
                         ! (units vary through calls)
,vw(nlev,n_dp)         & ! V-comp stress profile (N/m2)
,uw_base(nlev,n_dp)    & ! Cloud base U stress (N/m2)
,vw_base(nlev,n_dp)    & ! Cloud base V stress (N/m2)
,ue_p(nlev,n_dp)       & ! Environment U profile (m/s)
,ve_p(nlev,n_dp)         ! Environment V profile (m/s)

REAL :: exk_temp                ! Temporary exner

! Required by all version of CMT

REAL :: flxkp12(nlev,n_dp)      ! Mass flux on half level (Pa/s)

REAL :: mb(n_dp)                ! Cloud base mass flux (m/s or Pa/s depending
                                ! on location in code)

LOGICAL :: l_mom_gk             ! true if Gregory-Kershaw CMT required
LOGICAL :: l_mom_gk_stable      ! true if stabilized Gregory-Kershaw CMT
                                ! (different from the original) 

! Cape scaling/closure variables

INTEGER :: det_lev(n_dp)        ! Level at which split final
                                ! detrainment last occurred

INTEGER :: nterm                ! No. of points where conv.
                                ! has terminated

REAL :: tempnum                 ! Temporary variable for storage


! Arrays used by various CAPE closures
REAL ::                 &
  scale_f(n_dp)         & ! scale factor
 ,cape_ts_new(n_dp)     & ! Used as variable in RH-based closure
 ,relh(n_dp)            & ! RH integral (average when convection terminates)
 ,rh_mean(n_dp)         & ! mean RH over cloud on termination
 ,wls_mean(n_dp)        & ! wls - sum over convecting levels
 ,mass_mean(n_dp)       & ! mass - sum over convecting levels
 ,dptot(n_dp)             ! Delta P integral

! Original downdraught scheme variables

INTEGER :: npossdd              ! Max. no. of downdraughts
                                ! possible

INTEGER :: nnodd                ! No. of downdraughts not possible

INTEGER :: index_possdd(n_dp)   ! Index of downdraughts possible

INTEGER :: index_nodd(n_dp)     ! Index of downdraughts not
                                ! possible
INTEGER :: kmax_term            ! maximum termination level + 1

REAL :: deltap_cld              ! Pressure thickness of convective
                                ! cloud (Pa)
! New downdraught and evaporation scheme
! Uses npossdd to hold the number of real deep convecting points which are
! precipitating and so could have downdraughts
! index_possdd  the locations of these points  

REAL ::                  &
  precip_rain(n_dp,nlev) & ! Precipitation creation rate for rain (kg/m2/s)
 ,precip_snow(n_dp,nlev) & ! Precipitation creation rate for snow (kg/m2/s)
 ,tot_precip(n_dp)         ! total precipitation rate in a column (kg/m2/s)

! Local compressed arrays

LOGICAL :: bgmkp1_c(n_dp), bgmkp1_c2(n_dp) ! Mask for points
                                           ! where parcel in layer k+1
                                           ! is saturated

LOGICAL :: bwk_c(n_dp), bwk_c2(n_dp) ! bwater mask in layer k

LOGICAL :: bwkp1_c(n_dp), bwkp1_c2(n_dp) ! bwater mask in layer k+1

LOGICAL :: blatent_c2(n_dp)     ! Mask for points where latent heat has
                                ! been released

REAL :: thrk_c2(n_dp)           ! potential temperature of forced detrained air

REAL :: qrk_c2(n_dp)            ! specific humidity of forced detrained air

REAL :: deltak_c2(n_dp)         ! Parcel forced detrainment rate
                                ! in layer k multiplied by
                                ! appropriate layer thickness

REAL :: dqek_c2(n_dp)           ! Increment to q due to
                                ! convection in layer k (kg/kg)

REAL :: dqekp1_c2(n_dp)         ! Increment to q due to
                                ! convection in layer k+1 (kg/kg)

REAL :: dthek_c2(n_dp)          ! Increment to potential temp.
                                ! due to convection in layer k

REAL :: dthekp1_c2(n_dp)        ! Increment to potential temp.
                                ! due to convection in layer k+1

REAL :: dtraek_c2(n_dp,ntra)    ! Increment to model tracer due
                                ! to conv. at level k (kg/kg/s)

REAL :: dtraekp1_c2(n_dp,ntra)  ! Increment to model tracer due
                                ! to conv. at level k+1 (kg/kg/s)

REAL :: duek_c2(n_dp)           ! Increment to model U in layer k
                                ! due to CMT (m/s2)

REAL :: duekp1_c2(n_dp)         ! Increment to model U in layer
                                ! k+1 due to CMT (m/s2)

REAL :: dvek_c2(n_dp)           ! Increment to model V in layer k

REAL :: dvekp1_c2(n_dp)         ! Increment to model V in layer
                                ! k+1 due to CMT (m/s2)

REAL :: flxk_c(n_dp), flxk_c2(n_dp) !Parcel mass flux in layer k (Pa/s)

REAL :: flx_init_c2(n_dp)       ! Initial parcal mass flux at cloud base (Pa/s)

REAL :: flxkp12_c2(n_dp)        ! Half level mass flux (Pa/s)

REAL :: flxkp1_c2(n_dp)         ! Parcel mass flux in layer k+1
                                ! (Pa/s)

REAL :: prekp1_c2(n_dp)         ! Precip. from parcel as it rises
                                ! from layer k to k+1 (kg/m2/s)

REAL :: qpk_c(n_dp), qpk_c2(n_dp) ! Parcel mixing ratio in
                                  ! layer k(kg/kg)

REAL :: qpk(n_dp)

REAL :: qpkp1_c(n_dp), qpkp1_c2(n_dp) ! Parcel mixing ratio
                                      ! in layer k+1 (kg/kg)

REAL :: qek_c(n_dp), qek_c2(n_dp) ! Env. mixing ratio in
                                  ! layer k (kg/kg)

REAL :: qek(n_dp)
REAL :: qekp1_c(n_dp), qekp1_c2(n_dp) ! Env. mixing ratio in
                                      ! layer k+1 (kgkg-1)

REAL :: qekp1(n_dp)
REAL :: qsek_c2(n_dp)           ! Saturation mixing ratio of
                                ! cld. env. in layer k (kg/kg)

REAL :: qsekp1_c2(n_dp)         ! Saturation mixing ratio
                                ! of cld. env. in layer k+1 (kg/kg)

REAL :: thek_c(n_dp), thek_c2(n_dp) ! Env. potential temp in layer k (K)

REAL :: thek(n_dp)
REAL :: thekp1_c(n_dp), thekp1_c2(n_dp) ! Env. potential temp i in layer k (K)

REAL :: thekp1(n_dp)
REAL :: thpk_c(n_dp), thpk_c2(n_dp) ! Parcel potential temp in layer k (K)

REAL :: thpk(n_dp)
REAL :: thpkp1_c(n_dp), thpkp1_c2(n_dp)! Parcel potential temp in layer k (K)

REAL :: thpkp1(n_dp)

REAL :: traek_c2(n_dp,ntra)  ! Env. Tracer content in layer k (kg/kg)
REAL :: traekp1_c2(n_dp,ntra)! Env. Tracer content in layer k+1 (kg/kg)
REAL :: trapk_c2(n_dp,ntra)  ! Parcel Tracer content in layer k (kg/kg)
REAL :: trapkp1_c2(n_dp,ntra)! Parcel Tracer content in layer k+1 (kg/kg)

REAL :: rbuoyk_c(n_dp), rbuoyk_c2(n_dp)       ! Par. buoyancy at k (K)
REAL :: rbuoykp1_c(n_dp),rbuoykp1_c2(n_dp)    ! Par. buoyancy at k+1 (K)

REAL :: watldek_c(n_dp), watldek_c2(n_dp)     ! Env. water loading
                                              ! in layer k (kg/kg)
REAL :: watldpk_c(n_dp), watldpk_c2(n_dp)     ! Par. water loading
                                              ! in layer k (kg/kg)
REAL :: watldekp1_c(n_dp), watldekp1_c2(n_dp) ! Env. water loading
                                              ! in layer k+1 (kg/kg)
REAL :: watldpkp1_c(n_dp), watldpkp1_c2(n_dp) ! Par. water loading
                                              ! in layer k+1 (kg/kg)

REAL :: Qlkp1_c(n_dp),   &   ! Amount of condensation to liquid water
        Qlkp1_c2(n_dp)       ! in the parcel (kg/kg)
REAL :: Qfkp1_c(n_dp),   &   ! Amount of deposition to ice water
        Qfkp1_c2(n_dp)       ! in the parcel (kg/kg)
REAL :: Frezkp1_c(n_dp), &   ! Amount of freezing from liquid
        Frezkp1_c2(n_dp)     ! to frozen water in the parcel (kg/kg)

REAL :: uek_c2(n_dp)    ! Env. U in layer k (m/s)
REAL :: uekp1_c2(n_dp)  ! Env. U in layer k+1 (m/s)
REAL :: vek_c2(n_dp)    ! Env. V in layer k (m/s)
REAL :: vekp1_c2(n_dp)  ! Env. V in layer k+1 (m/s)
REAL :: upk_c2(n_dp)    ! Parcel U in layer k (m/s)
REAL :: upkp1_c2(n_dp)  ! Parcel U in layer k+1 (m/s)
REAL :: vpk_c2(n_dp)    ! Parcel V in layer k (m/s)
REAL :: vpkp1_c2(n_dp)  ! Parcel V in layer k+1 (m/s)

!====================================================
! Local compressed arrays for calculation of
! parcel vertical velocity

! Compressed arrays for (Parcel vertical velocity)^2 on ...
REAL :: w2p_km1_c2 (n_dp)   ! ...layer centre k-1 [(m/s)^2]
REAL :: w2p_k_c2   (n_dp)   ! ...layer centre k   [(m/s)^2]
REAL :: w2p_kp1_c2 (n_dp)   ! ...layer centre k+1 [(m/s)^2]

! Compressed arrays for height above surface of model level ...
REAL :: zkm1_c2    (n_dp)   ! ...k-1   [m]
REAL :: zk_c2      (n_dp)   ! ...k     [m]
REAL :: zkp12_c2   (n_dp)   ! ...k+1/2 [m]
REAL :: zkp1_c2    (n_dp)   ! ...k+1   [m]
!====================================================

! PC2 compression arrays

REAL :: qclek_c(n_dp), qclek_c2(n_dp) ! Environment liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qclekp1_c(n_dp), qclekp1_c2(n_dp) ! Environment liquid
                                          ! condensate mixing ratio in
                                          ! layer k+1 (kg/kg)

REAL :: qcfek_c(n_dp), qcfek_c2(n_dp) ! Environment frozen
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qcfekp1_c(n_dp), qcfekp1_c2(n_dp) ! Environment frozen
                                          ! condensate mixing ratio in
                                          ! layer k+1 (kg/kg)

REAL :: qclpk_c(n_dp), qclpk_c2(n_dp) ! Parcel liquid
                                      ! condensate mixing ratio in
                                      ! layer k (kg/kg)

REAL :: qclpkp1_c(n_dp), qclpkp1_c2(n_dp) ! Parcel liquid
                                          ! condensate mixing ratio in
                                          ! layer k+1 (kg/kg)

REAL :: qcfpk_c(n_dp), qcfpk_c2(n_dp) ! Parcel frozen
                                ! condensate mixing ratio in
                                ! layer k (kg/kg)

REAL :: qcfpkp1_c(n_dp), qcfpkp1_c2(n_dp) ! Parcel frozen
                                ! condensate mixing ratio in
                                ! layer k+1 (kg/kg)

REAL :: cflek_c2(n_dp),cflekp1_c2(n_dp) ! Environment liquid water
                                ! cloud volume ( )

REAL :: cffek_c2(n_dp),cffekp1_c2(n_dp) ! Environment frozen water
                                ! cloud volume ( )

REAL :: bcfek_c2(n_dp),bcfekp1_c2(n_dp) ! Environment bulk total
                                ! cloud volume ( )

REAL :: dqclek_c2(n_dp),dqclekp1_c2(n_dp) ! Environment increments
                                ! to liquid condensate mixing
                                ! ratio to convection (kg/kg/s)

REAL :: dqcfek_c2(n_dp),dqcfekp1_c2(n_dp) ! Environment increments
                                ! to frozen condensate mixing
                                ! ratio to convection (kg/kg/s)

REAL :: dcflek_c2(n_dp),dcflekp1_c2(n_dp) ! Environment increments
                                ! to liquid water cloud volume due
                                ! to convection (/s)

REAL :: dcffek_c2(n_dp),dcffekp1_c2(n_dp) ! Environment increments
                                ! to frozen water cloud volume due
                                ! to convection (/s)

REAL :: dbcfek_c2(n_dp),dbcfekp1_c2(n_dp) ! Environment increments
                                ! to bulk total cloud volume due
                                ! to convection (/s)

REAL :: amdetk_c2(n_dp)
LOGICAL :: bgmk_c2(n_dp)
LOGICAL :: bland_c2(n_dp)
LOGICAL :: blowst_c2(n_dp)
LOGICAL :: bterm_c2(n_dp)
REAL :: cape_c2(n_dp)
REAL :: cca_2d_c2(n_dp)
REAL :: cclwp_c2(n_dp)
REAL :: ccwkp1_c2(n_dp)
REAL :: dcpbydt_c2(n_dp)
REAL :: delpk_c2(n_dp)
REAL :: delpkp12_c2(n_dp)
REAL :: delpkp1_c2(n_dp)
REAL :: delp_uv_k_c2(n_dp)
REAL :: delp_uv_kp1_c2(n_dp)
REAL :: depth_c2(n_dp)
REAL :: dptot_c2(n_dp)
REAL :: eflux_u_ud_c2(n_dp)
REAL :: eflux_v_ud_c2(n_dp)
REAL :: ekp14_c(n_dp),ekp14_c2(n_dp)
REAL :: ekp34_c(n_dp),ekp34_c2(n_dp)
REAL :: exk_c(n_dp), exk_c2(n_dp)
REAL :: exkp1_c(n_dp),exkp1_c2(n_dp)
REAL :: expi_c2(n_dp)
INTEGER :: icct_c2(n_dp)
INTEGER :: iccb_c2(n_dp)
INTEGER :: lctop_c2(n_dp)
INTEGER :: lcbase_c2(n_dp)
REAL :: lcca_c2(n_dp)
REAL :: max_cfl_c2(n_dp)
REAL :: pk_c(n_dp),pk_c2(n_dp)
REAL :: pkp1_c(n_dp),pkp1_c2(n_dp)
REAL :: pstar_c2(n_dp)
REAL :: qpi_c2(n_dp)
REAL :: relh_c2(n_dp)
REAL :: tcw_c2(n_dp)
REAL :: thpi_c2(n_dp)
REAL :: xsbmin_v_c2(n_dp)
LOGICAL :: b_nodd(n_dp)   ! points with no downdraught
LOGICAL :: b_dd(n_dp)     ! points with downdraught on termination

! required by water conservation check

REAL ::                                                           &
 qminincolumn(n_dp)
                        ! Minimum value for q in column(kg/kg)
REAL :: temp1(n_dp)     ! work array
REAL :: temp2(n_dp)     ! work array

REAL ::      &
  rh_test    &   ! critical RH value for convective closure option
                 ! cldbase_opt_dp == 6
                 ! (RH-based CAPE timescale, timestep limited, reduced by w)
 ,rh_fac         ! factor for calculation in the above closure

! required by CMT

REAL ::   &
  wcld(n_dp)    &  ! Convective veloicty scale
 ,zlcl(n_dp)       ! lifting condensation level

REAL, PARAMETER :: qmin = 1.0e-8 ! Global minimum allowed Q


INTEGER :: warning   !  local integer (-ve) to pass to ereport

! Loop counters


INTEGER :: i,j,k,ktra,kt,kk

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DEEP_CONV_6A'


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise error_point

error_point=0


!-----------------------------------------------------------------------
! 2.0  Array Initialisation
!-----------------------------------------------------------------------

! Initialise logicals

DO i = 1,n_dp
  blowst(i)    = .TRUE.
  bterm(i)     = .FALSE.
  bconv(i)     = .FALSE.
  bcposs(i)    = .FALSE.
  b_nodd(i)    = .FALSE.
  b_dd(i)      = .FALSE.
  blatent(i)   = .FALSE.
  bgmk_term(i) = .FALSE.
END DO

DO i = 1,n_dp
  ind_cape_reduced(i) = 0.0
  cape_ts_used(i)     = 0.0
  cfl_limited(i)      = 0.0
  kterm(i)            = 0
  ind_deep(i)         = 0.0
  scale_f(i)          = 0.0
  start_lev(i)        = ntml(i)
END DO

!-----------------------------------------------------------------------
! 2.1  Initialise parcel properties and increment arrays
!-----------------------------------------------------------------------

!intialise parcel values over all levels
DO k = 1, nlev
  DO i = 1, n_dp
    qp(i,k)     = 0.0
    thp(i,k)    = 0.0
    qclp(i,k)   = 0.0
    qcfp(i,k)   = 0.0
    flx(i,k)    = 0.0
    precip(i,k) = 0.0
    ccw(i,k)    = 0.0
    area_ud(i,k) = 0.0
    area_dd(i,k) = 0.0
  END DO
END DO

IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k=1,nlev
      DO i = 1,n_dp
        trap(i,k,ktra) = 0.0
      END DO
    END DO
  END DO
END IF

DO k = 1,nlev
  DO i = 1,n_dp
    dthbydt(i,k)  = 0.0
    dqbydt(i,k)   = 0.0
    dqclbydt(i,k) = 0.0
    dqcfbydt(i,k) = 0.0
    dbcfbydt(i,k) = 0.0
    dcflbydt(i,k) = 0.0
    dcffbydt(i,k) = 0.0
    dt_dd(i,k)    = 0.0
    dq_dd(i,k)    = 0.0
    du_dd(i,k)    = 0.0
    dv_dd(i,k)    = 0.0
  END DO
END DO

IF (l_mom) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      dubydt(i,k) = 0.0
      dvbydt(i,k) = 0.0
    END DO
  END DO
END IF  ! L_mom

IF (l_tracer) THEN
  DO ktra = 1,ntra
    DO k = 1,nlev
      DO i = 1,n_dp
        dtrabydt(i,k,ktra) = 0.0
      END DO
    END DO
  END DO
END IF  ! L_tracer

!-----------------------------------------------------------------------
! 2.2  Initialise diagnostic arrays selected by STASH flags
!-----------------------------------------------------------------------
IF (flg_up_flx .OR. flg_mf_deep) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      up_flux(i,k)      = 0.0
    END DO
  END DO
END IF
IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      up_flux_half(i,k) = 0.0
    END DO
  END DO
END IF
IF (flg_dwn_flx) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      dwn_flux(i,k)     = 0.0
    END DO
  END DO
END IF
IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      entrain_up(i,k)   = 0.0
    END DO
  END DO
END IF
IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      detrain_up(i,k)   = 0.0
    END DO
  END DO
END IF
IF (flg_entr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      entrain_dwn(i,k)  = 0.0
    END DO
  END DO
END IF
IF (flg_detr_dwn) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      detrain_dwn(i,k)  = 0.0
    END DO
  END DO
END IF

IF (l_mom) THEN
  IF (flg_uw_dp) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        uw_deep(i,k)    = 0.0
      END DO
    END DO
  END IF
  IF (flg_vw_dp) THEN
    DO k = 1,nlev
      DO i = 1,n_dp
        vw_deep(i,k)    = 0.0
      END DO
    END DO
  END IF
END IF  ! L_mom

!-----------------------------------------------------------------------
! 2.3  Initialise radiation diagnostics
!-----------------------------------------------------------------------
DO i = 1,n_dp
  cca_2d(i) = 0.0
  cca_2d_term(i) = 0.0
  dcpbydt_term(i) = 0.0
  iccb(i)   = 0
  icct(i)   = 0
  tcw(i)    = 0.0
  cclwp(i)  = 0.0
  lcca(i)   = 0.0
  lctop(i)  = 0
  lcbase(i) = 0
END DO

DO k= 1,n_cca_lev
  DO i= 1,n_dp
    cca(i,k) = 0.0
  END DO
END DO


DO i = 1,n_dp
  !-----------------------------------------------------------------------
  ! 2.5  Initialise diagnostics for scaling and closure calculations
  !-----------------------------------------------------------------------
  flx_init(i)     = 0.0
  flx_init_new(i) = 0.0
  cape(i)         = 0.0
  cape_out(i)     = 0.0
  dcpbydt(i)      = 0.0
  max_cfl(i)      = 0.0
  det_lev(i)      = 0
  relh(i)         = 0.0
  dptot(i)        = 0.0
  wls_mean(i)     = 0.0
  mass_mean(i)    = 0.0

  !-----------------------------------------------------------------------
  ! 2.6  Initialise eddy flux arrays for updraught
  !-----------------------------------------------------------------------
  eflux_u_ud(i)   = 0.0
  eflux_v_ud(i)   = 0.0

  !-----------------------------------------------------------------------
  ! 2.7  Initialise surface precipitation arrays
  !-----------------------------------------------------------------------
  rain(i)         = 0.0
  snow(i)         = 0.0

  !-----------------------------------------------------------------------
  ! Initialise dummy variables that are not used by deep
  !-----------------------------------------------------------------------
  wsc_o_mb(i)     = 0.0
END DO

!Initialise adaptive entrainment variables
!initialise to level 2 'cos that's where parcel lift starts from
DO i = 1,n_dp
  thpk(i) = th(i,2)
  qpk(i)  = q(i,2)
END DO


!-----------------------------------------------------------------------
! Calculate parcel perturbations
!-----------------------------------------------------------------------

! Calculate XSBMIN and THPIXS constants based on layer thickness (Pa)
DO k = 1,nlev-1
  DO i = 1,n_dp
    xsbmin_v(i,k) = MIN( ((p_layer_centres(i,k) -               &
              p_layer_centres(i,k+1))/5000.0),1.0) *0.2

    thpixs_v(i,k) = MIN( ((p_layer_centres(i,k) -               &
              p_layer_centres(i,k+1))/5000.0),1.0) * thpixs_deep

    qpixs_v(i,k)  = qpixs_deep
  END DO
END DO  ! nlev

! Calculate cloud base mass flux
IF ( cldbase_opt_dp == 8 .OR. cldbase_opt_dp == 9 ) THEN
  ! Closure options that use the surface-flux based mass-flux at cloud-base

  ! Use value from CRM fits to deep simulations
  ! gives a value different to c_mass
  ! mb in units of m/s at this point.
  DO i = 1,n_dp
    mb(i) = b_cb*a_cb * wstar(i)
  END DO

ELSE ! original value
  DO i = 1,n_dp
    mb(i) = c_mass * wstar(i)
  END DO
END IF

! Define the LCL at the half level above ntml. Find environmental
! T at p_lcl by approximating theta there with
! th(i,k) + constant*(th(i,k+1)-th(i,k))  where constant is tunable.
! Similarly for q.

DO i = 1,n_dp
  k=ntml(i)
  p_lcl(i)  = p_layer_boundaries(i,k)
  th_lcl(i) = th(i,k) + 0.1 * (th(i,k+1) - th(i,k))
  t_lcl(i)  = th_lcl(i) * ((p_lcl(i) / pref)**kappa)
  q_lcl(i)  = q(i,k) + 0.1 * (q(i,k+1) - q(i,k))
END DO

! Calculate saturation mixing ratio at LCL

IF ( l_new_qsat_conv ) THEN
  CALL qsat_new(qse_lcl,t_lcl,p_lcl,n_dp)
ELSE
  ! DEPENDS ON: qsat_mix
  CALL qsat_mix(qse_lcl,t_lcl,p_lcl,n_dp,.FALSE.)
END IF

!-----------------------------------------------------------------------
! Initialize arrays required for Convective Momentum Transport(CMT)
!-----------------------------------------------------------------------
IF (l_mom) THEN
  SELECT CASE (deep_cmt_opt)

  CASE (2,6)         ! Gregory-Kershaw CMT

    ! need level near surface for initial parcel U & V values
    ! zsurf = 0.1*z_lcl

    DO i = 1,n_dp
      zsurf(i)  = 0.1*z_rho(i,ntml(i))
    END DO
    DO k=nlev-1,1,-1
      DO i = 1,n_dp
        IF (zsurf(i) <= z_theta(i,k)) THEN
          nstart(i) = k
        END IF
      END DO
    END DO
    l_mom_gk = .TRUE.

    ! Initialise winds for Gregory-Kershaw parcel calculation 
    DO k=1,nlev
      DO i = 1,n_dp
        up(i,k) = 0.0
        vp(i,k) = 0.0
      END DO
    END DO

    IF(deep_cmt_opt==6) THEN  ! Stabilized version
      l_mom_gk_stable = .TRUE.
    ELSE                      ! Original version 
      l_mom_gk_stable = .FALSE.
    END IF

  CASE DEFAULT    ! (0/1/5) Alan Grant's eddy viscosity based CMT

    ! Note: In terms of array indices p and phalf follow the convention
    !       used in the boundary layer scheme. phalf(k,*) refers to the
    !       lower boundary of uv layer k. This follows the convention for
    !       um UM4.5 and before

    !       Also note that p_layer_boundaries(0) and p_layer_centres(0)
    !       = pstar, so p_uv(k,1) and phalf_uv(k,1) will be equal.

    !       Because of the definition of nlcl, the pressure of the top of
    !       the mixed layer is phalf_uv(nlcl,*)


    k=1
    DO i = 1,n_dp
      p_uv(k,i)     = p_layer_boundaries(i,k-1)
      phalf_uv(k,i) = p_layer_centres(i,k-1)
      ue_p(k,i)     = u(i,k)
      ve_p(k,i)     = v(i,k)
      flxkp12(k,i)  = 0.0
      nlcl_uv(i)    = ntml(i) + 1
      n_0degc(i)    = freeze_lev(i)
    END DO

    DO i = 1,n_dp
      DO k = 2,nlev
        p_uv(k,i)     = p_layer_boundaries(i,k-1)
        phalf_uv(k,i) = p_layer_centres(i,k-1)
        ue_p(k,i)     = u(i,k)
        ve_p(k,i)     = v(i,k)
        flxkp12(k,i)  = 0.0
        exk_temp      = (p_uv(k,i)/pref)**kappa
        rho_uv(k,i)   = 2.0 * p_uv(k,i) / (r * exk_temp *        &
                        (th(i,k-1) + th(i,k)))
      END DO
      plcl_uv(i)      = phalf_uv(nlcl_uv(i),i)
      p_0degc_uv(i)   = phalf_uv(n_0degc(i),i)
      rho_uv(1,i)     = rho_uv(2,i)
    END DO

    l_mom_gk = .FALSE.
    l_mom_gk_stable = .FALSE.

  END SELECT      ! test on deep_cmt_opt

ELSE

  ! Initialise variable
  l_mom_gk = .FALSE.
  l_mom_gk_stable = .FALSE.

END IF     !L_mom

! Calculate theta and q perturbation (perturbation is based on
! environment buoyancy gradient)
! Reset thpixs and qpixs at ntml
DO i = 1,n_dp
  IF (t_lcl(i) >  tm) THEN
    dq_sat_env  = repsilon * lc * qse_lcl(i)                      &
                / (r * t_lcl(i) * t_lcl(i))
  ELSE
    dq_sat_env  = repsilon * (lc+lf) * qse_lcl(i)                 &
                / (r * t_lcl(i) * t_lcl(i))
  END IF

  b_calc      = t_lcl(i) * c_virtual * dq_sat_env + 1.0           &
              + c_virtual * qse_lcl(i)

  thv_pert(i) = -0.5 * (th(i,ntml(i)+1)                           &
              * (1.0+c_virtual * q(i,ntml(i)+1))                  &
              - th(i,ntml(i)) * (1.0 + c_virtual                  &
              * q(i,ntml(i)))) + 0.5

  c_calc    = th_lcl(i) * c_virtual * (qse_lcl(i)                 &
            - q_lcl(i)) - thv_pert(i)

  thpert(i) = MAX(MIN(-c_calc / b_calc, max_dp_thpert),           &
              min_dp_thpert) ! ignore term in thpert**2

  thpixs_v(i,ntml(i)) = thpert(i)

  qpert(i)  = MAX(MIN(qse_lcl(i) + ((p_lcl(i) / pref)             &
            **kappa) * thpert(i) * dq_sat_env                     &
            - q_lcl(i),                                           &
              max_dp_qpert_fac * qse_lcl(i)),0.0)

  qpixs_v(i,ntml(i))  = qpert(i)

END DO ! n_dp


! Set bwater=.true. on points where water will condense rather than
! ice.
CALL flag_wet(n_dp,n_dp,nlev,th,exner_layer_centres,bwater)

!-----------------------------------------------------------------------
! 3.0  Main loop over all levels
!-----------------------------------------------------------------------

DO k = 2,nlev-1
  !-----------------------------------------------------------------------
  ! Initialise environment variables
  ! NB These variable are only used by layer_cn.
  !-----------------------------------------------------------------------
  DO i = 1,n_dp
    thek(i)   = th(i,k)
    qek(i)    = q(i,k)
    thekp1(i) = th(i,k+1)
    qekp1(i)  = q(i,k+1)
    !Note that unlike p_layer_boundaries, where k indexing is offset
    !by one compared to the dynamics numbering, z retains the numbering
    !convention for dynamics variables i.e. for theta levels, k->k
    !and for rho levels k+1/2 -> k+1
    zkm1(i)   = z_theta(i,k-1)
    zk(i)     = z_theta(i,k)
    zkp12(i)  = z_rho(i,k+1)
    zkp1(i)   = z_theta(i,k+1)
    rhum(i)   = q(i,k) / qse(i,k)
  END DO

  !-----------------------------------------------------------------------
  ! Initialise parcel properties (theta,q,tracer,momentum) if convection
  ! is not occurring at level k and has not convected in column before
  !-----------------------------------------------------------------------
  DO i = 1,n_dp
    IF ( .NOT. bconv(i) .AND. det_lev(i) == 0) THEN
      expi(i)     = exner_layer_centres(i,k)
      bgmk(i)     = .FALSE.
      depth(i)    = 0.0
      thpi(i)     = th(i,k) + thpixs_v(i,k)
      thp(i,k)    = thpi(i)
      thpk(i)     = thp(i,k)
      qpi(i)      = q(i,k)  + qpixs_v(i,k)
      qp(i,k)     = qpi(i)
      qpk(i)      = qp(i,k)
      IF (l_q_interact) THEN
        qclp(i,k) = qcl(i,k)
        qcfp(i,k) = qcf(i,k)
      ELSE
        qclp(i,k) = 0.0
        qcfp(i,k) = 0.0
      END IF
      IF (l_mom_gk) THEN  ! Gregory Kershaw CMT
        ! Set initial parcel values at cloud base to values of near surface winds
        up(i,k)   = u(i,nstart(i))
        vp(i,k)   = v(i,nstart(i))
      END IF
    END IF
  END DO  ! n_dp
  IF (l_tracer) THEN
    DO ktra=1,ntra
      DO i = 1,n_dp
        IF ( .NOT. bconv(i)) THEN
          trap(i,k,ktra)  = tracer(i,k,ktra)
        END IF  !not bconv
      END DO
    END DO
  END IF

  !-----------------------------------------------------------------------
  ! 3.1  Calculate layer dependent constants (pressure,
  !      layer thickness, entrainment coefficients, detrainment
  !      coefficients)
  !-----------------------------------------------------------------------

  CALL layer_cn_6a(k, n_dp, nlev,                                   &
                   mdet_dp_on,                                      &
                   ntml, ntpar, start_lev,                          &
                   exner_layer_boundaries, exner_layer_centres,     &
                   p_layer_boundaries, p_layer_centres,             &
                   z_rho,                                           &
                   conv_prog_precip,                                &
                   recip_pstar, entrain_coef, rhum,                 &
                   zk, zkp12, zkp1,                                 &
                   thek, qek, thekp1, qekp1,                        &
                   thpk, qpk,                                       &
                   wsc_o_mb, qsat_lcl, w_max,                       &
                   deep,                                            &
                   bconv,                                           &
                   ! Out
                   pk, pkp1, exk, exkp1,                            &
                   delpk, delpkp12, delpkp1,                        &
                   delp_uv_k, delp_uv_kp1,                          &
                   ekp14, ekp34, amdetk                             &
                   )

  ! Maximum initial convective mass flux
  DO i = 1,n_dp
    flxmax(i) = delpk(i) / ((1.0 + ekp14(i)) * timestep)
  END DO  ! n_dp

  !-----------------------------------------------------------------------
  ! Initial test to check if convection is possible in layer k
  !-----------------------------------------------------------------------
  ! Convection is possible if
  ! - the point was convecting (bconv = .T.) and did not terminate
  !   in the previous layer
  ! - or if at the top level of the surface mixed layer (k = ntml)
  DO i = 1,n_dp
    bcposs(i) = bconv(i) .OR. k  ==  ntml(i)
  END DO  ! n_dp

  ! Calculate number of points which may convect (ncposs) and
  ! set compression indices (index1)
  ncposs = 0
  DO i = 1,n_dp
    IF (bcposs(i)) THEN
      ncposs          = ncposs + 1
      index1(ncposs)  = i
    END IF
  END DO

  !-----------------------------------------------------------------------
  ! Compress points where convection may occur
  ! NB This process is used to update some single level variables that are
  ! defined on level k and kp1 by using the full field variables.
  ! NB The order in which the variables are compressed
  ! is the same as the argument list for LIFT_PAR
  !-----------------------------------------------------------------------
  IF (ncposs  >   0) THEN
    DO i = 1,ncposs
      !INTENT(IN) for lift_par/water_loading
      thek_c(i)     = th(index1(i),k)
      thekp1_c(i)   = th(index1(i),k+1)
      qek_c(i)      = q(index1(i),k)
      qekp1_c(i)    = q(index1(i),k+1)
      qclek_c(i)    = qcl(index1(i),k)
      qcfek_c(i)    = qcf(index1(i),k)
      qclekp1_c(i)  = qcl(index1(i),k+1)
      qcfekp1_c(i)  = qcf(index1(i),k+1)
      pk_c(i)       = pk(index1(i))
      pkp1_c(i)     = pkp1(index1(i))
      exkp1_c(i)    = exkp1(index1(i))
      thpk_c(i)     = thp(index1(i),k)
      qpk_c(i)      = qp(index1(i),k)
      qclpk_c(i)    = qclp(index1(i),k)
      qcfpk_c(i)    = qcfp(index1(i),k)
      ekp14_c(i)    = ekp14(index1(i))
      ekp34_c(i)    = ekp34(index1(i))
      bwk_c(i)      = bwater(index1(i),k)
      bwkp1_c(i)    = bwater(index1(i),k+1)
      !INTENT(IN) for water_loading only
      exk_c(i)      = exk(index1(i))
    END DO
  END IF  ! ncposs>0

  !-----------------------------------------------------------------------
  ! 3.2  Lift parcel from layer k to layer k+1
  !-----------------------------------------------------------------------
  IF (ncposs > 0) THEN

    CALL lift_par_6a(ncposs, thek_c, thekp1_c,                      &
                qek_c, qekp1_c, qclek_c, qcfek_c,                   &
                qclekp1_c, qcfekp1_c,                               &
                pk_c, pkp1_c, exkp1_c,                              &
                thpk_c, qpk_c, qclpk_c, qcfpk_c,                    &
                ekp14_c, ekp34_c,                                   &
                l_q_interact, bwk_c, bwkp1_c,                       &
                !Out
                bgmkp1_c, thpkp1_c, qpkp1_c,                        &
                qclpkp1_c, qcfpkp1_c,                               &
                Qlkp1_c, Qfkp1_c, Frezkp1_c)

    !-----------------------------------------------------------------------
    ! Calculate the water loading for level k and k+1
    !-----------------------------------------------------------------------
    CALL water_loading(ncposs, pk_c, exk_c, thek_c, thpk_c,         &
                     qclek_c, qcfek_c, qclpk_c, qcfpk_c,            &
                     watldek_c, watldpk_c)

    CALL water_loading(ncposs, pkp1_c, exkp1_c, thekp1_c, thpkp1_c, &
                     qclekp1_c, qcfekp1_c, qclpkp1_c, qcfpkp1_c,    &
                     watldekp1_c, watldpkp1_c)

    ! NEC compiler directive
    !CDIR NODEP

    !-----------------------------------------------------------------------
    ! Test if convection is starting from layer k
    !-----------------------------------------------------------------------
    DO i = 1,ncposs ! Loop over points which may convect

      ! Calculate buoyancy (virt. potential temp.) of parcel in layer k and k+1
      rbuoyk_c(i)   = thpk_c(i) * (1.0 + c_virtual *qpk_c(i)        &
                    - watldpk_c(i))                                 &
                    - thek_c(i) * (1.0 + c_virtual *qek_c(i)        &
                    - watldek_c(i))

      rbuoykp1_c(i) = thpkp1_c(i) * (1.0 + c_virtual *qpkp1_c(i)    &
                    - watldpkp1_c(i))                               &
                    - thekp1_c(i) * (1.0 + c_virtual *qekp1_c(i)    &
                    - watldekp1_c(i))

      ! Allow parcel to convect from ntml.
      IF (k  ==  ntml(index1(i))) THEN
        bconv(index1(i))  = .TRUE.  ! convection active
        blowst(index1(i)) = .TRUE.  ! convection initialised in layer

        ! Set parcel mass flux
        flxk_c(i)         = mb(index1(i)) * g                       &
                          * p_layer_centres(index1(i),k)            &
                          / ( r * thpk_c(i)                         &
                          * (p_layer_centres(index1(i),k)/ pref)**kappa )

        ! If mass flux out of the initial layer is greater than the mass flux
        ! of the layer over the timestep then limit mass flux to mass of layer.
        IF (flxk_c(i)  >   flxmax(index1(i))) THEN
          flxk_c(i) = flxmax(index1(i))
        END IF

        ! Write compressed mass flux back to full array
        flx(index1(i),k) = flxk_c(i)

        ! Store diagnostics linked to initial convective mass flux for
        ! calculation of final closure.
        flx_init(index1(i))    = flxk_c(i)
        flxmax_init(index1(i)) = flxmax(index1(i))

      ELSE     ! k=ntml test
        blowst(index1(i)) = .FALSE. ! convection not initialise in layer
      END IF   ! k=ntml test


      ! Reset threshold for forced detrainment to the initial
      ! (positive or negative) buoyancy (limit positive buoy.
      ! threshold to XSBMIN fn(delta P)), ONLY for first 5 levels of lift

      IF (k  >=  ntml(index1(i)) .AND.                              &
          k  <=  ntml(index1(i)) + 4) THEN

        xsbmin_v(index1(i),k) = MIN ( xsbmin_v(index1(i),k),        &
                                      thv_pert(index1(i)))
      END IF

    END DO  !ncposs
  END IF    !ncposs>0

  ! Calculate number of points which are convecting  (nconv)
  ! set compression indices (index2).
  nconv = 0
  DO i = 1,ncposs
    IF (bconv(index1(i))) THEN
      nconv         = nconv + 1
      index2(nconv) = i
    END IF
  END DO

  !-----------------------------------------------------------------------
  ! Second compression to form arrays of length nconv to be passed
  ! into CONVEC2.
  ! NB This process is used to update some single level variables that are
  ! defined on level k and kp1 by using the full field variables.
  ! NB The order in which the variables are compressed
  ! is the same as the argument list for CONVEC2.
  !-----------------------------------------------------------------------
  IF (nconv  >   0) THEN
    !Compression for INTENT(IN)
    DO i = 1,nconv
      start_lev_c2(i)   = ntml(index1(index2(i)))
      pstar_c2(i)       = pstar(index1(index2(i)))
      pk_c2(i)          = pk_c(index2(i))
      pkp1_c2(i)        = pkp1_c(index2(i))
      delpk_c2(i)       = delpk(index1(index2(i)))
      delpkp1_c2(i)     = delpkp1(index1(index2(i)))
      delpkp12_c2(i)    = delpkp12(index1(index2(i)))
      delp_uv_k_c2(i)   = delp_uv_k(index1(index2(i)))
      delp_uv_kp1_c2(i) = delp_uv_kp1(index1(index2(i)))
      exk_c2(i)         = exk_c(index2(i))
      exkp1_c2(i)       = exkp1_c(index2(i))
      thek_c2(i)        = thek_c(index2(i))
      thekp1_c2(i)      = thekp1_c(index2(i))
      qek_c2(i)         = qek_c(index2(i))
      qekp1_c2(i)       = qekp1_c(index2(i))
      qclek_c2(i)       = qclek_c(index2(i))
      qclekp1_c2(i)     = qclekp1_c(index2(i))
      qcfek_c2(i)       = qcfek_c(index2(i))
      qcfekp1_c2(i)     = qcfekp1_c(index2(i))
      qsek_c2(i)        = qse(index1(index2(i)),k)
      qsekp1_c2(i)      = qse(index1(index2(i)),k+1)
      cflek_c2(i)       = cf_liquid(index1(index2(i)),k)
      cflekp1_c2(i)     = cf_liquid(index1(index2(i)),k+1)
      cffek_c2(i)       = cf_frozen(index1(index2(i)),k)
      cffekp1_c2(i)     = cf_frozen(index1(index2(i)),k+1)
      bcfek_c2(i)       = bulk_cf(index1(index2(i)),k)
      bcfekp1_c2(i)     = bulk_cf(index1(index2(i)),k+1)
      thpk_c2(i)        = thpk_c(index2(i))
      qpk_c2(i)         = qpk_c(index2(i))
      qclpk_c2(i)       = qclpk_c(index2(i))
      qcfpk_c2(i)       = qcfpk_c(index2(i))
      thpi_c2(i)        = thpi(index1(index2(i)))
      qpi_c2(i)         = qpi(index1(index2(i)))
      expi_c2(i)        = expi(index1(index2(i)))
      rbuoyk_c2(i)      = rbuoyk_c(index2(i))
      rbuoykp1_c2(i)    = rbuoykp1_c(index2(i))
      xsbmin_v_c2(i)    = xsbmin_v(index1(index2(i)),k)
      watldek_c2(i)     = watldekp1_c(index2(i))
      watldekp1_c2(i)   = watldekp1_c(index2(i))
      watldpk_c2(i)     = watldpkp1_c(index2(i))
      watldpkp1_c2(i)   = watldpkp1_c(index2(i))
      Qlkp1_c2(i)       = Qlkp1_c(index2(i))
      Qfkp1_c2(i)       = Qfkp1_c(index2(i))
      Frezkp1_c2(i)     = Frezkp1_c(index2(i))
      ekp14_c2(i)       = ekp14_c(index2(i))
      ekp34_c2(i)       = ekp34_c(index2(i))
      amdetk_c2(i)      = amdetk(index1(index2(i)))
      flxk_c2(i)        = flx(index1(index2(i)),k)
      flx_init_c2(i)    = flx_init(index1(index2(i)))
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        uek_c2(i)       = u(index1(index2(i)),k)
        uekp1_c2(i)     = u(index1(index2(i)),k+1)
        vek_c2(i)       = v(index1(index2(i)),k)
        vekp1_c2(i)     = v(index1(index2(i)),k+1)
        upk_c2(i)       = up(index1(index2(i)),k)
        vpk_c2(i)       = vp(index1(index2(i)),k)
      END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,nconv
          traek_c2(i,ktra)   = tracer(index1(index2(i)),k,ktra)
          traekp1_c2(i,ktra) = tracer(index1(index2(i)),k+1,ktra)
          trapk_c2(i,ktra)   = trap(index1(index2(i)),k,ktra)
        END DO
      END DO
    END IF
    IF (flg_w_eqn .OR. l_new_dd) THEN
      DO i=1, nconv
        zkm1_c2(i)      = zkm1  (index1(index2(i)))
        zk_c2(i)        = zk    (index1(index2(i)))
        zkp12_c2(i)     = zkp12 (index1(index2(i)))
        zkp1_c2(i)      = zkp1  (index1(index2(i)))
        w2p_km1_c2(i)   = w2p(index1(index2(i)),k-1)
      END DO
    END IF
    DO i = 1,nconv
      bgmk_c2(i)        = bgmk(index1(index2(i)))
      bgmkp1_c2(i)      = bgmkp1_c(index2(i))
      bwk_c2(i)         = bwk_c(index2(i))
      bwkp1_c2(i)       = bwkp1_c(index2(i))
      blowst_c2(i)      = blowst(index1(index2(i)))
      bland_c2(i)       = bland(index1(index2(i)))
    END DO
    !Compression for INTENT(INOUT)
    DO i = 1,nconv
      lcbase_c2(i)      = lcbase(index1(index2(i)))
      lctop_c2(i)       = lctop(index1(index2(i)))
      thpkp1_c2(i)      = thpkp1_c(index2(i))
      qpkp1_c2(i)       = qpkp1_c(index2(i))
      qclpkp1_c2(i)     = qclpkp1_c(index2(i))
      qcfpkp1_c2(i)     = qcfpkp1_c(index2(i))
      dthek_c2(i)       = dthbydt(index1(index2(i)),k)
      dqek_c2(i)        = dqbydt(index1(index2(i)),k)
      dqclek_c2(i)      = dqclbydt(index1(index2(i)),k)
      dqcfek_c2(i)      = dqcfbydt(index1(index2(i)),k)
      tcw_c2(i)         = tcw(index1(index2(i)))
      depth_c2(i)       = depth(index1(index2(i)))
      cclwp_c2(i)       = cclwp(index1(index2(i)))
      lcca_c2(i)        = lcca(index1(index2(i)))
      cape_c2(i)        = cape(index1(index2(i)))
      dcpbydt_c2(i)     = dcpbydt(index1(index2(i)))
      relh_c2(i)        = relh(index1(index2(i)))
      dptot_c2(i)       = dptot(index1(index2(i)))
      max_cfl_c2(i)     = max_cfl(index1(index2(i)))
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        eflux_u_ud_c2(i)= eflux_u_ud(index1(index2(i)))
        eflux_v_ud_c2(i)= eflux_v_ud(index1(index2(i)))
        duek_c2(i)      = dubydt(index1(index2(i)),k)
        dvek_c2(i)      = dvbydt(index1(index2(i)),k)
      END DO
    END IF
    IF (l_tracer) THEN
      DO ktra = 1,ntra
        DO i = 1,nconv
          dtraek_c2(i,ktra) = dtrabydt(index1(index2(i)),k,ktra)
        END DO
      END DO
    END IF
    IF (flg_w_eqn .OR. l_new_dd) THEN
      DO i=1, nconv
        w2p_k_c2(i)     = w2p(index1(index2(i)),k)
      END DO
    END IF
    DO i = 1,nconv
      bterm_c2(i)       = .FALSE.
      blatent_c2(i)     = blatent(index1(index2(i)))
    END DO
    !Compression for INTENT(OUT)
    ! Most INTENT(OUT) variables do not need to be initialised because they
    ! are always set in CONVEC2 or if they are not set in CONVEC2 then they
    ! are not used. However, several of the cloud variables do need to be
    ! initialised.
    DO i = 1,nconv
      iccb_c2(i)        = iccb(index1(index2(i)))
      icct_c2(i)        = icct(index1(index2(i)))
      cca_2d_c2(i)      = cca_2d(index1(index2(i)))
    END DO

    !-----------------------------------------------------------------------
    ! 3.3  Calculate the rest of the parcel ascent  and the effect of
    !      convection on the large-scale atmosphere.

    !      Subroutine CONVEC2

    !      UM Documentation paper 27, sections (5),(6),(7),(8),(9),(10)
    !-----------------------------------------------------------------------

    CALL convec2_6a  (k, nconv, n_dp, nlev, ntra, dp_on, dp_new_termc,       &
                      start_lev_c2,                                          &
                      pstar_c2, pk_c2, pkp1_c2, delpk_c2,                    &
                      delpkp1_c2, delpkp12_c2, delp_uv_k_c2, delp_uv_kp1_c2, &
                      exk_c2, exkp1_c2,                                      &
                      thek_c2, thekp1_c2, qek_c2, qekp1_c2,                  &
                      qclek_c2, qclekp1_c2, qcfek_c2, qcfekp1_c2,            &
                      qsek_c2, qsekp1_c2,                                    &
                      cflek_c2, cflekp1_c2,  cffek_c2,  cffekp1_c2,          &
                      bcfek_c2,  bcfekp1_c2,                                 &
                      thpk_c2, qpk_c2, qclpk_c2, qcfpk_c2,                   &
                      thpi_c2, qpi_c2, expi_c2,                              &
                      rbuoyk_c2, rbuoykp1_c2,                                &
                      watldek_c2, watldekp1_c2, watldpk_c2, watldpkp1_c2,    &
                      Qlkp1_c2, Qfkp1_c2, Frezkp1_c2,                        &
                      ekp14_c2, ekp34_c2, amdetk_c2, flxk_c2, flx_init_c2,   &
                      uek_c2, uekp1_c2, vek_c2, vekp1_c2,                    &
                      upk_c2, vpk_c2,                                        &
                      traek_c2, traekp1_c2, trapk_c2,                        &
                      zkm1_c2, zk_c2, zkp12_c2, zkp1_c2, w2p_km1_c2,         &
                      l_q_interact, l_mom_gk, l_mom_gk_stable, l_tracer,     &
                      bgmk_c2, bgmkp1_c2, bwk_c2,                            &
                      bwkp1_c2, blowst_c2, bland_c2,                         &

                      ! In/out
                      lcbase_c2, lctop_c2,                                   &
                      thpkp1_c2, qpkp1_c2, qclpkp1_c2, qcfpkp1_c2,           &
                      dthek_c2, dqek_c2, dqclek_c2, dqcfek_c2,               &
                      tcw_c2, depth_c2, cclwp_c2, lcca_c2,                   &
                      cape_c2, dcpbydt_c2, relh_c2, dptot_c2, max_cfl_c2,    &
                      eflux_u_ud_c2, eflux_v_ud_c2,                          &
                      duek_c2, dvek_c2,                                      &
                      dtraek_c2, w2p_k_c2,                                   &
                      bterm_c2, blatent_c2, xsbmin_v_c2,                     &

                      ! Out
                      iccb_c2, icct_c2,                                      &
                      dcflek_c2, dcffek_c2, dbcfek_c2,                       &
                      dthekp1_c2, dqekp1_c2, dqclekp1_c2, dqcfekp1_c2,       &
                      dcflekp1_c2, dcffekp1_c2, dbcfekp1_c2,                 &
                      prekp1_c2, thrk_c2, qrk_c2, deltak_c2,                 &
                      flxkp12_c2, flxkp1_c2,                                 &
                      cca_2d_c2, ccwkp1_c2,                                  &
                      upkp1_c2, vpkp1_c2,                                    &
                      duekp1_c2, dvekp1_c2,                                  &
                      trapkp1_c2,                                            &
                      dtraekp1_c2, w2p_kp1_c2 )

  END IF ! nconv > 0

  !-----------------------------------------------------------------------
  ! Decompression of compressed variables coming out of
  ! of CONVEC2.
  ! NB The order in which the variables are decompressed
  ! is the same as the the argument list for CONVEC2.
  !-----------------------------------------------------------------------
  DO i = 1,n_dp
    depth(i)      = 0.0
    bgmk(i)       = .FALSE.
    bterm(i)      = .FALSE.
  END DO

  IF (nconv  >   0) THEN
    !Decompression for INTENT(IN)
    DO i = 1,nconv
      bgmk(index1(index2(i)))         = bgmkp1_c2(i)
    END DO
    !Decompression for INTENT(INOUT)
    DO i = 1,nconv
      lcbase(index1(index2(i)))       = lcbase_c2(i)
      lctop(index1(index2(i)))        = lctop_c2(i)
      thp(index1(index2(i)),k+1)      = thpkp1_c2(i)
      qp(index1(index2(i)),k+1)       = qpkp1_c2(i)
      qclp(index1(index2(i)),k+1)     = qclpkp1_c2(i)
      qcfp(index1(index2(i)),k+1)     = qcfpkp1_c2(i)
      dthbydt(index1(index2(i)),k)    = dthek_c2(i)
      dqbydt(index1(index2(i)),k)     = dqek_c2(i)
      dqclbydt(index1(index2(i)),k)   = dqclek_c2(i)
      dqcfbydt(index1(index2(i)),k)   = dqcfek_c2(i)
      tcw(index1(index2(i)))          = tcw_c2(i)
      depth(index1(index2(i)))        = depth_c2(i)
      cclwp(index1(index2(i)))        = cclwp_c2(i)
      lcca(index1(index2(i)))         = lcca_c2(i)
      cape(index1(index2(i)))         = cape_c2(i)
      dcpbydt(index1(index2(i)))      = dcpbydt_c2(i)
      relh(index1(index2(i)))         = relh_c2(i)
      dptot(index1(index2(i)))        = dptot_c2(i)
      max_cfl(index1(index2(i)))      = max_cfl_c2(i)
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        eflux_u_ud(index1(index2(i))) = eflux_u_ud_c2(i)
        eflux_v_ud(index1(index2(i))) = eflux_v_ud_c2(i)
        dubydt(index1(index2(i)),k)   = duek_c2(i)
        dvbydt(index1(index2(i)),k)   = dvek_c2(i)
      END DO
    END IF
    IF (l_tracer) THEN
      DO i = 1,nconv
        DO ktra = 1,ntra
          dtrabydt(index1(index2(i)),k,ktra)    = dtraek_c2(i,ktra)
        END DO
      END DO
    END IF
    IF (flg_w_eqn .OR. l_new_dd) THEN
      DO i=1, nconv
        w2p(index1(index2(i)),k)      = w2p_k_c2(i)
      END DO
    END IF
    DO i = 1,nconv
      bterm(index1(index2(i)))        = bterm_c2(i)
      blatent(index1(index2(i)))      = blatent_c2(i)
    END DO
    !Decompression for INTENT(OUT)
    DO i = 1,nconv
      iccb(index1(index2(i)))         = iccb_c2(i)
      icct(index1(index2(i)))         = icct_c2(i)
      dcflbydt(index1(index2(i)),k)   = dcflek_c2(i)
      dcffbydt(index1(index2(i)),k)   = dcffek_c2(i)
      dbcfbydt(index1(index2(i)),k)   = dbcfek_c2(i)
      dthbydt(index1(index2(i)),k+1)  = dthekp1_c2(i)
      dqbydt(index1(index2(i)),k+1)   = dqekp1_c2(i)
      dqclbydt(index1(index2(i)),k+1) = dqclekp1_c2(i)
      dqcfbydt(index1(index2(i)),k+1) = dqcfekp1_c2(i)
      dcflbydt(index1(index2(i)),k+1) = dcflekp1_c2(i)
      dcffbydt(index1(index2(i)),k+1) = dcffekp1_c2(i)
      dbcfbydt(index1(index2(i)),k+1) = dbcfekp1_c2(i)
      precip(index1(index2(i)),k+1)   = prekp1_c2(i)
      !dummy                          = deltak_c2(i)
    END DO
    IF (l_mom) THEN    ! needed for all versions
      DO i = 1,nconv
        flxkp12(k,index1(index2(i)))  = flxkp12_c2(i)
      END DO
    END IF  ! L_mom
    DO i = 1,nconv
      flx(index1(index2(i)),k+1)      = flxkp1_c2(i)
      cca_2d(index1(index2(i)))       = cca_2d_c2(i)
      ccw(index1(index2(i)),k+1)      = ccwkp1_c2(i)
    END DO
    IF (l_mom_gk) THEN
      DO i = 1,nconv
        up(index1(index2(i)),k+1)     = upkp1_c2(i)
        vp(index1(index2(i)),k+1)     = vpkp1_c2(i)
        dubydt(index1(index2(i)),k+1) = duekp1_c2(i)
        dvbydt(index1(index2(i)),k+1) = dvekp1_c2(i)
      END DO
    END IF
    IF (l_tracer) THEN
      DO i = 1,nconv
        DO ktra = 1,ntra
          trap(index1(index2(i)),k+1,ktra)      = trapkp1_c2(i,ktra)
          dtrabydt(index1(index2(i)),k+1,ktra)  = dtraekp1_c2(i,ktra)
        END DO
      END DO
    END IF
    IF (flg_w_eqn .OR. l_new_dd) THEN
      DO i=1, nconv
        w2p(index1(index2(i)),k+1)    = w2p_kp1_c2(i)
      END DO
    END IF

  END IF      ! nconv > 0

  !-----------------------------------------------------------------------
  ! 3.4  Cape and CFL scaling - adjust initial mass flux so that cape is
  !      removed by convection over timescale cape_timescale.
  !-----------------------------------------------------------------------

  ! Set up integer nterm which is the total number of points where
  ! convection has terminated.


  nterm = 0
  DO i = 1,n_dp
    IF (bterm(i)) THEN
      nterm       = nterm + 1
      bgmk_term(i)= bgmk(i)
      rh_mean(i)  = relh(i)/dptot(i)
      kterm(i)    = k
      ! If convection has terminated write cape to diagnostic output
      ! variable (cape_out).
      cape_out(i) = cape(i)
      cape(i)     = 0.0
      bconv(i)    = .FALSE.
      det_lev(i)  = k+1
      DO kk=ntml(i),k
        wls_mean(i)  = wls_mean(i)  + w(i,kk)*r2rho_th(i,kk)*dr_across_th(i,kk)
        mass_mean(i) = mass_mean(i) + r2rho_th(i,kk)*dr_across_th(i,kk)
      END DO
      wls_mean(i) =  wls_mean(i)/mass_mean(i)
    END IF
  END DO

  !-----------------------------------------------------------------------
  ! Write out entrainment, detrainment and half-level mass flux diagnostics.
  ! They will be scaled by the full level mass flux outside
  ! of the level loop
  !-----------------------------------------------------------------------
  ! Calculate fractional entrainment rate for level k.
  IF (flg_entr_up) THEN
    DO i = 1,nconv
      entrain_up(index1(index2(i)),k) = (1.0 - deltak_c2(i))        &
               * (1.0 - amdetk_c2(i)) * (ekp14_c2(i) + ekp34_c2(i)  &
               * (1.0 + ekp14_c2(i)))
    END DO
  END IF

  ! Calculate fractional detrainment rate for level k
  IF (flg_detr_up) THEN
    DO i = 1,nconv
      detrain_up(index1(index2(i)),k) = -(amdetk_c2(i)              &
                      + deltak_c2(i) * (1.0 - amdetk_c2(i)))
    END DO
  END IF

  ! Calculate the half level mass flux for level k
  ! Only the scaling factor between full level and half levels is calculated
  ! here. This is scaled by the full level mass flux outside the level loop
  IF (flg_up_flx_half) THEN
    DO i = 1,nconv
      up_flux_half(index1(index2(i)),k) = (1.0 - deltak_c2(i))      &
                * (1.0 - amdetk_c2(i)) * (1.0 + ekp14_c2(i))
    END DO
  END IF

  ! Copy SCM profile diagnostics for which profiles are not stored elsewhere
  IF (l_scm_convss_dg) THEN
    DO i = 1, nconv
      j = index1(index2(i))

      ! Save buoyancy profiles at level k
      scm_convss_dg(j) % par_thetav_excess(k) = rbuoyk_c2(i)
      scm_convss_dg(j) % par_thetav(k)                                    &
           = thpk_c2(i) * ( 1.0 + c_virtual*qpk_c2(i) - watldpk_c2(i) )
      scm_convss_dg(j) % env_thetav(k)                                    &
           = thek_c2(i) * ( 1.0 + c_virtual*qek_c2(i) - watldek_c2(i) )

      ! If this is the termination level, also save these at k+1
      IF ( bterm_c2(i) ) THEN
        scm_convss_dg(j) % par_thetav_excess(k+1) = rbuoykp1_c2(i)
        scm_convss_dg(j) % par_thetav(k+1)                                &
           = thpkp1_c2(i) * ( 1.0 + c_virtual*qpkp1_c2(i) - watldpkp1_c2(i) )
        scm_convss_dg(j) % env_thetav(k+1)                                &
           = thekp1_c2(i) * ( 1.0 + c_virtual*qekp1_c2(i) - watldekp1_c2(i) )
      END IF

      ! Entrainment and mixing detrainment profiles
      scm_convss_dg(j) % ekp14(k)        = ekp14_c2(i)
      scm_convss_dg(j) % ekp34(k)        = ekp34_c2(i)
      scm_convss_dg(j) % amdetk(k)       = amdetk_c2(i)

      ! Save adaptive detrainment diagnostics at level k / k+1
      scm_convss_dg(j) % deltak(k)       = deltak_c2(i)
      scm_convss_dg(j) % rbuoy_star(k+1) = rbuoykp1_c2(i)
      scm_convss_dg(j) % xsbmin(k+1)     = xsbmin_v_c2(i)
      IF (k==ntml(i)) THEN
        scm_convss_dg(j) % rbuoy_star(k) = rbuoyk_c2(i)
        scm_convss_dg(j) % xsbmin(k)     = rbuoyk_c2(i)
      END IF
      scm_convss_dg(j) % thrk(k)         = thrk_c2(i)
      scm_convss_dg(j) % qrk(k)          = qrk_c2(i)
      ! buoyancy excess of detrained air, assuming same waterloading as parcel
      IF ( thrk_c2(i) > 0.0 ) THEN
        scm_convss_dg(j) % thvrk_excess(k)                                &
            = thrk_c2(i) * ( 1.0 + c_virtual*qrk_c2(i) - watldpk_c2(i) )  &
            - thek_c2(i) * ( 1.0 + c_virtual*qek_c2(i) - watldek_c2(i) )
      END IF

      ! Save the updraft mass-flux before closure scaling
      scm_convss_dg(j) % up_flx_guess(k) = flx(j,k)

    END DO
  END IF

  !-----------------------------------------------------------------------
  ! 3.6  End of main loop over levels
  !-----------------------------------------------------------------------
END DO


! If used, copy convection profile diagnostics
IF (l_scm_convss_dg) THEN

  ! Copy profile diagnostics
  DO k = 1, nlev
    DO i = 1, n_dp
      IF ( k>=ntml(i) .AND. k<=det_lev(i) ) THEN
        scm_convss_dg(i) % par_theta(k) = thp(i,k)
        scm_convss_dg(i) % par_q(k)     = qp(i,k)
        scm_convss_dg(i) % par_qcl(k)   = qclp(i,k)
        scm_convss_dg(i) % par_qcf(k)   = qcfp(i,k)
      END IF
      scm_convss_dg(i) % env_theta(k) = th(i,k)
      scm_convss_dg(i) % env_q(k)     = q(i,k)
      scm_convss_dg(i) % env_qcl(k)   = qcl(i,k)
      scm_convss_dg(i) % env_qcf(k)   = qcf(i,k)
    END DO
  END DO

END IF


!-----------------------------------------------------------------------
! 4.0 Choice of cloud base closure option
!-----------------------------------------------------------------------
SELECT CASE ( cldbase_opt_dp )

  ! Default 4a convection scheme - RH-based CAPE closure
CASE ( 0 )

  DO i = 1,n_dp
    IF (dcpbydt(i)  >   0.0) THEN
      cape_ts_new(i) =                                                    &
                MIN(MAX(900.0*(1.0 - rh_mean(i))/0.1,60.0),cape_timescale)
      IF (cape_ts_new(i) < cape_timescale) THEN
        ind_cape_reduced(i) = 1.0
      END IF
    END IF
  END DO

  ! Modified 4a convection scheme - RH-based CAPE closure
  ! timescale limited to timestep
CASE ( 1 )

  DO i = 1,n_dp
    IF (dcpbydt(i)  >   0.0) THEN
      cape_ts_new(i) =                                                   &
                  MIN(MAX(cape_timescale*(1.0-rh_mean(i))/0.4,timestep)  &
                ,cape_timescale)
      IF (cape_ts_new(i) < cape_timescale) THEN
        ind_cape_reduced(i) = 1.0
      END IF
    END IF  ! dcpbydt > 0
  END DO

  ! Fixed cape timescale
CASE ( 2 )

  DO i = 1,n_dp
    IF (dcpbydt(i)  >   0.0) THEN
      cape_ts_new(i) =  cape_timescale
    END IF
  END DO

  ! w based cape closure; if w_cape_limit > 1000. reverts to cape_timescale
CASE ( 3 )
  ! Initialise array to cape timescale and then alter as required.
  cape_ts_new(:) =  cape_timescale

  IF ( w_cape_limit < 1000.0 ) THEN
    !  This section includes test on w_max
    DO i = 1,n_dp
      IF ( dcpbydt(i) > 0.0 ) THEN
        ! new denominator introduced at vn6.6
        IF ( w_max(i) > w_cape_limit ) THEN
          cape_ts_new(i) =   cape_timescale * w_cape_limit/           &
                    (w_cape_limit+ (w_max(i)-w_cape_limit)*wcape_fac)
          ! set indicator that CAPE reduced
          ind_cape_reduced(i) = 1.0
        END IF !  w_max(i) > w_cape_limit
      END IF  ! dcpbydt > 0
    END DO
  END IF  ! w_cape_limit

  ! w based cape closure with grid-box area scaling
CASE ( 4 )

  IF ( w_cape_limit < 1000.0 ) THEN
    !  This section includes test on w_max
    DO i = 1,n_dp
      IF ( dcpbydt(i) > 0.0 ) THEN
        IF ( w_max(i) > w_cape_limit ) THEN
          cape_ts_new(i) = cape_timescale * w_cape_limit/w_max(i)
        ELSE
          cape_ts_new(i) = cape_timescale * cape_out(i) / cape_min +    &
                             cape_timescale * EXP(-cape_out(i) / cape_min)
        END IF !  w_max(i) > w_cape_limit
      END IF  ! dcpbydt > 0
    END DO
  ELSE
    DO i = 1,n_dp
      IF ( dcpbydt(i) > 0.0 ) THEN
        cape_ts_new(i) = cape_timescale * cape_out(i) / cape_min +      &
                         cape_timescale * EXP( - cape_out(i) /cape_min)
      END IF  ! dcpbydt > 0
    END DO
  END IF  ! w_cape_limit

  ! w based cape closure (experimental option; not used); 
CASE ( 5 )

  IF ( w_cape_limit < 1000.0 ) THEN
    !  This section includes test on w_max
    DO i = 1,n_dp
      IF ( dcpbydt(i) > 0.0 ) THEN
        IF ( w_max(i) > w_cape_limit ) THEN
          cape_ts_new(i) =   cape_timescale * w_cape_limit/w_max(i)
        ELSE
          IF ( rh_mean(i) >= 0.75 ) THEN
            cape_ts_new(i) = cape_timescale *( 0.2373 / (rh_mean(i))**5)
            ind_cape_reduced(i) = 1.0
          ELSE
            cape_ts_new(i) = cape_timescale
          END IF ! rh_mean(i) >= 0.75
        END IF   !  w_max(i) > w_cape_limit
      END IF  ! dcpbydt > 0
    END DO
  ELSE
    DO i = 1,n_dp
      IF ( dcpbydt(i) > 0.0 ) THEN
        IF ( rh_mean(i) >= 0.75 ) THEN
          cape_ts_new(i) = cape_timescale *( 0.2373 / (rh_mean(i))**5)
        ELSE
          cape_ts_new(i) = cape_timescale
        END IF ! rh_mean(i) >= 0.75
      END IF   ! dcpbydt > 0
    END DO
  END IF  ! w_cape_limit

  ! RH and w based CAPE option
  ! Expects a sensible w_cape_limit or will do nothing
CASE ( 6 )

  rh_test = 0.60         ! critical RH value
  rh_fac  = 1.0/ (1.0 - rh_test)

  IF ( w_cape_limit < 1000.0 ) THEN
    !  This section includes test on w_max
    DO i = 1,n_dp
      IF ( dcpbydt(i) > 0.0 ) THEN
        ! work out any reduction in cape_timescale due to RH
        ! linearly falls to 1/2 given cape_timescale for RH above rh_test
        IF ( rh_mean(i) >= rh_test ) THEN
          cape_ts_new(i) = cape_timescale*0.5*                          &
                             (1.0+(1.0-rh_mean(i))*rh_fac)
          ind_cape_reduced(i) =1.0
        ELSE
          cape_ts_new(i) = cape_timescale
        END IF
        ! Further reduction if w_max above critical value
        IF ( w_max(i) > w_cape_limit ) THEN
          cape_ts_new(i) = cape_ts_new(i) * w_cape_limit/               &
                          (w_cape_limit + (w_max(i)-w_cape_limit)*wcape_fac)
          ind_cape_reduced(i) =1.0
        END IF
        ! Limit CAPE timescale to convective timestep
        cape_ts_new(i)  = MAX(cape_ts_new(i), timestep)
      END IF  ! dcpbydt > 0
    END DO
  END IF  ! w_cape_limit

  ! CAPE timescale dependent on large-scale w,
  ! or
  ! CAPE timescale dependent on large-scale w but with a lower limit of a 
  !   cloud base mass flux of a_cb*b_cb*wstar (applied later)
CASE ( 7, 8 )

  DO i = 1,n_dp
    IF ( dcpbydt(i) > 0.0 ) THEN
      IF (wls_mean(i) > 0.0) THEN
        cape_ts_new(i) = a_cape *(wls_mean(i)**b_cape)
      ELSE    ! set to maximum value
        cape_ts_new(i) = max_cape
      END IF

      ! Limit CAPE timescale to a maximum of 4 hours
      cape_ts_new(i)  = MIN(cape_ts_new(i), max_cape)

      ! Limit CAPE timescale to a minimum of convective timestep
      cape_ts_new(i)  = MAX(cape_ts_new(i), timestep)
    END IF  ! dcpbydt > 0
  END DO

  ! Boundary layer and large-scale vertical velocity closure
  ! NO CAPE element to this closure
  ! Works out required new cloud base mass flux directly
CASE ( 9 )

  DO i = 1,n_dp
    IF ( dcpbydt(i) > 0.0 ) THEN
      IF (wls_mean(i) > 0.0) THEN
        ! Closure based on large-scale vertical velocity and wstar
        ! Currently flx_init(i) = a_cb*b_cb*wstar*rho_cb
        ! Instead of a fixed sigma_cb = b_cb we want to use
        ! sigma_cb = b_cb+c_cb*wls_mean
        flx_init_new(i) = flx_init(i)*(1.0+c_cb*wls_mean(i)/b_cb)
      ELSE    ! Assume a fixed fractional area
        flx_init_new(i) = flx_init(i)
      END IF
    END IF  ! dcpbydt > 0
  END DO

END SELECT        ! cldbase_opt_dp


!---------------------------------------------------------------------------
! Calculate closure scaling to apply to convective tendencies etc
!---------------------------------------------------------------------------
! If using some sort of CAPE-closure (all closure options except the bl / w one
IF (cldbase_opt_dp /= 9) THEN
  DO i = 1,n_dp
    IF (dcpbydt(i) > 0.0) THEN
      ! Calculate new mass-flux at cloud-base, so-as to remove
      ! CAPE over time cape_ts_new(i)
      flx_init_new(i)   = flx_init(i) * cape_out(i) /                 &
                          (cape_ts_new(i)*dcpbydt(i))
      ! Apply limit for the mass-flux out of the originating layer
      IF (flx_init_new(i) > flxmax_init(i)) THEN
        flx_init_new(i) = flxmax_init(i)
      END IF
      ! Scale max_cfl with cape scale
      max_cfl(i)        = max_cfl(i) * flx_init_new(i) / flx_init(i)
    ELSE
      ! If convective tendencies don't erode CAPE, set mass-flux to zero
      flx_init_new(i) = 0.0
    END IF  ! dcpbydt > 0
  END DO
END IF

! If using surface-limited w-based CAPE timescale:
IF (cldbase_opt_dp == 8) THEN
  DO i = 1,n_dp
    IF ( flx_init_new(i) > 0.0 .AND. flx_init(i) > flx_init_new(i) ) THEN

      ! Use original closure i.e. surface based, if the surface-based
      ! closure implies a larger mass-flux than the CAPE-based closure

      ! First Scale max_cfl with cape scale
      max_cfl(i)        = max_cfl(i) * flx_init(i) / flx_init_new(i)

      ! And then apply the surface based closure
      flx_init_new(i)   = flx_init(i)

      ! Apply limit for the mass-flux out of the originating layer
      IF (flx_init_new(i) > flxmax_init(i)) THEN
        flx_init_new(i) = flxmax_init(i)
      END IF

    END IF
  END DO
END IF

! Work out scaled mass flux needed to keep cfl ratio below limit.
! This applies whether CAPE or surface based closure.
DO i = 1,n_dp
  ! Scale max_cfl by the timestep
  max_cfl(i)        = max_cfl(i) * timestep

  ! If the max CFL in the profile exceeds the limit, reduce the mass-flux
  ! accordingly
  IF (max_cfl(i) > cfl_limit) THEN
    flx_init_new(i) = flx_init_new(i) * cfl_limit / max_cfl(i)
    cfl_limited(i)  = 1.0       ! flag that deep convection is CFL limited
  END IF

  ! Also re-apply separately calculated CFL-limit for the initiating layer.
  IF (flx_init_new(i)  >   flxmax_init(i)) THEN
    flx_init_new(i) = flxmax_init(i)
  END IF

  max_cfl(i) = 0.0
END DO

! Compute the final closure scaling and the effective CAPE timescale diagnostic
DO i = 1,n_dp   
  ! Calculate the final scaling.
  scale_f(i)      = flx_init_new(i) / flx_init(i)

  ! set flx_init to the new value to provide the real initial mass
  ! flux in all conditions
  flx_init(i)     = flx_init_new(i)

  IF (flx_init_new(i) > 0.0) THEN
    ! Work out effective CAPE timescale being used taking into account
    ! all restrictions being applied to cloud base mass flux (diagnostic).
    cape_ts_used(i) = cape_out(i)/(scale_f(i) * dcpbydt(i))
  END IF
END DO


!---------------------------------------------------------------------------
! Test to catch "failed convection" events
!---------------------------------------------------------------------------

! First, save SCM diagnostic of whether the convection failed, and if so, why:
IF ( l_scm_convss_dg ) THEN
  DO i = 1,n_dp

    ! Convection failed due to the ascent being too shallow or cloud-free
    IF ( icct(i)-iccb(i) <= 3 .OR. ( .NOT. blatent(i) ) ) THEN
      scm_convss_dg(i) % status_deep = 1

    ! Convection failed because the closure set the mass-flux to zero
    ELSE IF ( flx_init_new(i)  <=  minflx ) THEN
      scm_convss_dg(i) % status_deep = 2

    ! Real deep convection occurred!
    ELSE
      scm_convss_dg(i) % status_deep = 3

    END IF

  END DO
END IF

! Do the actual test for failed convection, and reset outputs to zero
! where convection has failed
DO i = 1,n_dp
  IF ( (flx_init_new(i)  <=  minflx)                        &
                  .OR. ( (icct(i)-iccb(i)) <= 3 )           &
                  .OR. ( .NOT. blatent(i) ) ) THEN
    ! False convection diagnosed if:
    ! - the new initial mass flux is less than zero
    ! - or the convecting layer it too thin
    ! - or the parcel has never released latent heat
    ! 3d variables are reset below by setting scale_f to zero.
    flx_init(i)     = 0.0
    flx_init_new(i) = 0.0
    mb(i)           = 0.0
    scale_f(i)      = 0.0
    cca_2d(i)       = 0.0
    iccb(i)         = 0
    icct(i)         = 0
    tcw(i)          = 0.0
    cclwp(i)        = 0.0
    lcca(i)         = 0.0
    lctop(i)        = 0
    lcbase(i)       = 0
    ! kterm should be reset to zero if convection fails but because it is used
    ! for reasons of efficiency in 3d variable scaling it is reset later
    ! instead here
    !   kterm(i)        = 0
    cape_ts_used(i) = 0.0     ! Not real deep so set to zero
    ind_cape_reduced(i) = 0   ! Ensure not set to 1 as not real deep
  ELSE
    ! True convection
    ind_deep(i)     = 1.0  ! real deep event indicator
  END IF
END DO


! Maximum termination level
kmax_term = 2
DO i = 1,n_dp
  IF (kterm(i)+1 >  kmax_term) THEN
    kmax_term = kterm(i)+1
  END IF
END DO

IF (kmax_term > nlev) THEN
  kmax_term = nlev
END IF
!-----------------------------------------------------------------------
! Apply cape and cfl scaling
!-----------------------------------------------------------------------
DO kt = 2, kmax_term
  DO i = 1,n_dp
    IF (kt  >=  ntml(i) .AND. kt <= kterm(i)+1 ) THEN

      dthbydt(i,kt)   = dthbydt(i,kt)  * scale_f(i)
      dqbydt(i,kt)    = dqbydt(i,kt)   * scale_f(i)
      dqclbydt(i,kt)  = dqclbydt(i,kt) * scale_f(i)
      dqcfbydt(i,kt)  = dqcfbydt(i,kt) * scale_f(i)
      dcflbydt(i,kt)  = dcflbydt(i,kt) * scale_f(i)
      dcffbydt(i,kt)  = dcffbydt(i,kt) * scale_f(i)
      dbcfbydt(i,kt)  = dbcfbydt(i,kt) * scale_f(i)

      IF (l_mom_gk) THEN
        dubydt(i,kt)  = dubydt(i,kt) * scale_f(i)
        dvbydt(i,kt)  = dvbydt(i,kt) * scale_f(i)
      END IF
      IF (l_mom) THEN     ! required for all versions
        flxkp12(kt,i) = flxkp12(kt,i) * scale_f(i)
      END IF
      IF (l_tracer) THEN
        DO ktra = 1,ntra
          dtrabydt(i,kt,ktra) = dtrabydt(i,kt,ktra)*scale_f(i)
        END DO
      END IF
      flx(i,kt)    = flx(i,kt) * scale_f(i)
      precip(i,kt) = precip(i,kt) * scale_f(i)

    END IF !kt>ntml and flx_init_new>0
  END DO  ! i loop
END DO  ! kt loop

!-----------------------------------------------------------------------
! Scale cloud fraction
! Additional check on scale_f needed because of logarithmic dependence
!-----------------------------------------------------------------------
DO i = 1,n_dp
  IF (scale_f(i) > 0.0) THEN
    cca_2d(i)  = cca_2d(i) + 0.06 * LOG(scale_f(i))

    ! Check scaled cloud fraction not smaller than minimum value
    ! (2.0E-5) or greater than unity.
    cca_2d(i) = MAX(2.0e-5, cca_2d(i))
    cca_2d(i) = MIN(1.0e+0, cca_2d(i))
  END IF
END DO      ! i

!-----------------------------------------------------------------------
! Reset kterm value now that increments etc for all levels have been set
! to zero in the cases of failed convection
!-----------------------------------------------------------------------
DO i = 1,n_dp
  IF (scale_f(i) == 0.0) THEN
    kterm(i) = 0
  END IF
END DO      ! i

!-----------------------------------------------------------------------
! Write out updraught massflux diagnostics and scale the
! entrainment and detrainment diagnostics by the mass flux.
!-----------------------------------------------------------------------
IF (flg_up_flx .OR. flg_mf_deep) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      up_flux(i,k) = flx(i,k)
    END DO
  END DO
END IF

IF (flg_up_flx_half) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      up_flux_half(i,k) = up_flux_half(i,k) * flx(i,k)
    END DO
  END DO
END IF

IF (flg_entr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      entrain_up(i,k) = entrain_up(i,k) * flx(i,k)
    END DO
  END DO
END IF

IF (flg_detr_up) THEN
  DO k = 1,nlev
    DO i = 1,n_dp
      detrain_up(i,k) = detrain_up(i,k) * flx(i,k)
    END DO
  END DO
END IF

IF (flg_area_ud .OR. l_new_dd) THEN  ! updraught area
  DO k = 1,nlev
    DO i = 1,n_dp
      IF (flx(i,k) > 0.0) THEN
        ! expect to get a valid wup
        ! updraught core area from wup and flx

        IF (w2p(i,k) > 0.0) THEN
          area_ud(i,k) =  flx(i,k)/(g*SQRT(w2p(i,k))* rho_theta(i,k))
        ELSE
          ! The calculation of w2p is not producing a sensible value
          ! Assume wcb =1m/s and calculate a value for cloud base as value
          ! equired for new DD scheme.
          IF (k == iccb(i)+1) THEN 
            area_ud(i,k) =  flx(i,k)/(g* rho_theta(i,k))
          END IF  
        END IF
      END IF
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
! 5.0  Mixes the increments from the initial parcel perturbation throughout
!      the subcloud layer
!-----------------------------------------------------------------------

CALL mix_ipert_6a(n_dp, nlev, nbl, ntml, p_layer_boundaries,    &
             exner_layer_centres, dthbydt, dqbydt, flx_init,    &
             thpert, qpert)


!-----------------------------------------------------------------------
! 6.0 Down draughts  - now 2 options
!                     (a) Emanuel downdraught scheme
!                     (b) Original mass flux code
!                     (c) New mass flux scheme 
! Note the level at which deep convection terminates has been stored
! in the above updraught loop as Kterm.
!-----------------------------------------------------------------------

IF (l_eman_dd) THEN

  !-----------------------------------------------------------------------
  ! (a) Emanuel downdraught scheme
  !-----------------------------------------------------------------------

  ! Work out maximum termination level
  kmax_term = 2
  DO i = 1,n_dp
    IF (kterm(i) >  kmax_term) THEN
      kmax_term = kterm(i)
    END IF
  END DO

  IF (l_snow_rain) THEN   ! revised Emanuel DD

    CALL eman_dd_rev (n_dp,kmax_term,nlev,trlev,ntra,                 &
                  kterm,l_tracer,                                     &
                  exner_layer_centres,                                &
                  p_layer_centres, p_layer_boundaries,                &
                  timestep, scale_f, th, q, qse, tracer, precip,      &
                  dthbydt, dqbydt, dtrabydt,                          &
                  rain, snow ,dwn_flux, dt_dd, dq_dd)

  ELSE

    CALL eman_dd (n_dp,kmax_term,nlev,trlev,ntra                  &
,                      kterm,l_tracer                             &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries        &
,                      timestep, th, q, qse, tracer, precip       &
,                      dthbydt, dqbydt, dtrabydt                  &
,                      rain, snow ,dwn_flux, dt_dd, dq_dd         &
                     )
  END IF

ELSE
   
  IF (l_new_dd) THEN
    !-----------------------------------------------------------------------
    ! (c) New Downdraught & evaporation scheme
    !-----------------------------------------------------------------------
    ! Work out maximum termination level and whether there is any
    ! precipitation from the updraught. If no precipitation is produced
    ! then there is no point in calling the scheme.
    
    DO i = 1,n_dp
      tot_precip(i) = 0.0
    END DO
    ! Sum precipitation in the column
    DO k = 1, nlev
      DO i = 1,n_dp
        tot_precip(i) = tot_precip(i) + precip(i,k)
      END DO
    END DO

    kmax_term = 2
    npossdd = 0
    DO i = 1,n_dp
      IF(kterm(i) >  kmax_term) THEN
        kmax_term = kterm(i)
      END IF
      ! Real deep convection with precipitation
      IF(kterm(i) > ntml(i) .AND. flx_init_new(i) >0.0 .AND.       &
                                               tot_precip(i) > 0.0 ) THEN
        npossdd = npossdd +1
        index_possdd(npossdd) = i
      END IF
    END DO
    IF (npossdd > 0) THEN
      ! Define precip_rain & precip_snow from precip
      ! May eventually come direct from updraught
      ! bwater(i,k) not defined for k=1
      DO k = 2,nlev
        DO i = 1,n_dp
          IF (bwater(i,k)) THEN                ! assumed all rain
            precip_rain(i,k) = precip(i,k)
            precip_snow(i,k) = 0.0
          ELSE                                 ! assumed all snow
            precip_snow(i,k) = precip(i,k)
            precip_rain(i,k) = 0.0
          END IF
        END DO
      END DO
      ! no precipitation from level 1 so  
      DO i = 1,n_dp
        precip_snow(i,1) = 0.0
        precip_rain(i,1) = 0.0
      END DO

      CALL dd_evap_ctl(n_dp, npossdd, nlev, trlev, ntra, kmax_term,        &
                       index_possdd, iccb, kterm, freeze_lev,              &
                       l_tracer, timestep, r_theta, r_rho,                 &
                       z_theta, z_rho, dr_across_th,dr_across_rh,          &
                       exner_layer_centres,exner_layer_boundaries,         &
                       p_layer_centres, p_layer_boundaries,                &
                       rho, rho_dry, rho_theta, rho_dry_theta,             &
                       th, q, u, v, qse,                                   &
                       precip_rain, precip_snow, flx, area_ud,             &
                       thp, qp, up, vp, trap, tracer,                      &
                       rain, snow, rain_3d, snow_3d,                       &
                       dthbydt, dqbydt, dubydt, dvbydt, dtrabydt,          & 
                       dwn_flux, area_dd, dt_dd, dq_dd, du_dd, dv_dd)
    END IF 
  ELSE 
    !-----------------------------------------------------------------------
    ! (b) Original mass flux downdraughts & evaporation scheme
    !-----------------------------------------------------------------------
    ! 6.1  Downdraught calculation - on all points where convection is
    !      terminating.
    !      UM Documentation Paper 27
    !-----------------------------------------------------------------------

    DO i = 1,n_dp
      IF (kterm(i) >= ntml(i) .AND. flx_init_new(i) >0.0) THEN
        k = kterm(i)
        tempnum = 0.0
        IF (iccb(i)  >   0) THEN
          deltap_cld = p_layer_centres(i,iccb(i))                     &
                           - p_layer_centres(i,k)
          DO kt = iccb(i), k+1
            tempnum = tempnum + precip(i,kt)
          END DO
        ELSE
          deltap_cld = 0.0
        END IF


        ! Downdraughts possible if pressure thickness of convective
        ! cloud (deltap_cld) is greater than 15000Pa, the point is saturated
        ! and the precip. in the layer is greater than a threshold
        ! value (1E-12).
        ! Set logical for use later
        IF (deltap_cld  >   15000.0 .AND. bgmk_term(i) .AND.             &
                               tempnum  >   1e-12) THEN
          b_dd(i) = .TRUE.
        ELSE
          b_nodd(i) = .TRUE.
        END IF

      END IF   ! test on whether deep really happened
    END DO
    !-----------------------------------------------------------------------
    ! 6.2  Downdraught calculation - on all points where convection is
    !      terminating.
    !      UM Documentation Paper 27
    !-----------------------------------------------------------------------

    npossdd = 0
    DO i = 1,n_dp
      IF (b_dd(i)) THEN
        npossdd = npossdd +1
        index_possdd(npossdd) = i
      END IF
    END DO

    IF (npossdd  >   0) THEN

      ! Work out maximum termination level
      kmax_term = 2
      DO i = 1,npossdd
        IF(kterm(index_possdd(i)) >  kmax_term) THEN
          kmax_term = kterm(index_possdd(i))
        END IF
      END DO

      CALL dd_all_call_6a (n_dp,npossdd,kmax_term,nlev,trlev,ntra &
,                      kterm, iccb, icct, index_possdd, l_tracer  &
,                      bwater(1,2)                                &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries,pstar  &
,                      recip_pstar,timestep , cca_2d              &
,                      thp, qp, th, q, qse, trap,tracer, flx,precip &
,                      dthbydt, dqbydt, dtrabydt                  &
,                      rain, snow , rain_3d, snow_3d, dwn_flux    &
,                      entrain_dwn, detrain_dwn, dt_dd, dq_dd)


    END IF

    !-----------------------------------------------------------------------
    ! 6.3 Surface precipitation calculation for terminating points with
    !     no downdraught (moved outside level loop) ie do this calculation
    !     on all points at the end.
    !-----------------------------------------------------------------------
    ! Points where no downdraught possible
    nnodd = 0
    DO i = 1,n_dp

      IF (b_nodd(i)) THEN
        nnodd = nnodd +1
        index_nodd(nnodd) = i
      END IF
    END DO

    IF (nnodd  >   0) THEN

      ! Work out maximum termination level
      kmax_term = 2
      DO i = 1,nnodd
        IF(kterm(index_nodd(i)) >  kmax_term) THEN
          kmax_term = kterm(index_nodd(i))
        END IF
      END DO
      ! Only add 1 if kmax_term is less than model levels (which should be
      ! true).
      IF (kmax_term  <  nlev ) THEN
        kmax_term = kmax_term + 1
      END IF

      ! Surface precipitation calculation


      CALL evap_bcb_nodd_all(n_dp,nnodd,kmax_term,kterm           &
,                      iccb, index_nodd, bwater(1,2)              &
,                      exner_layer_centres,exner_layer_boundaries &
,                      p_layer_centres, p_layer_boundaries,pstar  &
,                      timestep , cca_2d, th, q, qse, precip      &
,                      dthbydt, dqbydt                            &
,                      rain, snow, rain_3d, snow_3d, dt_dd, dq_dd )

    END IF

  END IF        ! choice of mass flux scheme
END IF        ! test on down draught type

!-----------------------------------------------------------------------
! 7.0  Convective Momentum Transport (if L_mom = .true.)
!-----------------------------------------------------------------------

IF (l_mom) THEN

  SELECT CASE (deep_cmt_opt)
  CASE (2,6)     ! Gregory-Kershaw deep CMT

    ! Do nothing here as calculated in parcel ascent earlier

  CASE DEFAULT   ! (0/1/5) Alan Grant's Eddy viscosity CMT

    ! altered to use kterm instead of ntpar
    DO i = 1,n_dp
      IF (kterm(i) >=  nlcl_uv(i)) THEN
        ntop_uv(i)    = kterm(i) + 1
      ELSE     ! case where deep convection fails
        ! I think in this case the cloud base mass flux will be zero so
        ! there will be no CMT. (The value will not matter)
        ntop_uv(i)    = ntpar(i) + 1
      END IF

      ptop_uv(i)    = phalf_uv(ntop_uv(i),i)
    END DO

    nterm = 0

    ! Set cloud base mass flux equal to mass flux at half level below
    ! the LCL. mb now completely reset and units changed to Pa/s

    DO i = 1,n_dp
      IF (mb(i) > 0.0 .AND. kterm(i) < nlev -1) THEN
        nterm = nterm + 1
        cu_term(nterm) = i
        cu_tend(nterm) = i
        mb(i) = flxkp12(nlcl_uv(i),i)
        DO j = 1,nlev
          flxkp12(j,i) = 0.0
        END DO
      ELSE IF (kterm(i) == nlev - 1) THEN
        ! Problem deep convection has gone to the top of the model
        ! Return location of deep problem point plus profiles
        IF (printstatus >= prstatus_normal) THEN
          warning = -1
          WRITE(umMessage,'(A)')                                         &
            "  PROBLEM: deep convection has gone to the top of the model."
          CALL umPrint(umMessage,src='deep_conv_6a')
          WRITE(umMessage,'(a23,i6,a7,i6,a7,i6)')                        &
              ' Deep convection point ',i,                               &
              ' kterm ',kterm(i),' nlcl  ',nlcl_uv(i)
          CALL umPrint(umMessage,src='deep_conv_6a')
          WRITE(umMessage,'(A6,4A26)') 'k','Theta','q','qcl','qcf'
          CALL umPrint(umMessage,src='deep_conv_6a')
          DO k=1,nlev
            WRITE(umMessage,'(I6,4G26.18)') k,th(i,k),q(i,k),qcl(i,k),qcf(i,k)
            CALL umPrint(umMessage,src='deep_conv_6a')
          END DO
          CALL ereport( 'deep_conv_6a', warning,                         &
            "Convection has gone to the top of the model."    //newline//&
            "This usually means there is junk in the inputs"  //newline//&
            "to the convection scheme.  Check the profiles at"//newline//&
            "the grid-point where this has occured (these are"//newline//&
            "printed in the run output file)." )
        END IF
        error_point = i
        IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, &
                                zhook_handle)
        RETURN

      END IF

      ! initialise output arrays as lower level subroutines don't set all
      ! values

      DO j = 1,nlev
        uw(j,i)       =0.0
        vw(j,i)       =0.0
        uw_base(j,i)  =0.0
        vw_base(j,i)  =0.0
        visc(j,i)     =0.0
        mass_dwn(j,i) =0.0
      END DO

    END DO

    IF (nterm  >   0) THEN
      CALL cmt_mass(n_dp, n_dp, nlev, nterm, cu_term,                 &
              kterm, cu_tend, n_0degc, nlcl_uv, ntop_uv,              &
              mb, p_0degc_uv, plcl_uv, ptop_uv, phalf_uv, p_uv,       &
            ! Output arguments
              flxkp12 ,mass_dwn, visc)

      CALL deep_grad_stress(n_dp,n_dp,n_dp,nlev,nlcl_uv,ntop_uv,         &
                            nterm,cu_term,cu_tend,                       &
                            ue_p,ve_p,visc,phalf_uv,p_uv,rho_uv,timestep,&
                            ! Output
                            uw,vw)


      CALL deep_ngrad_stress(n_dp,n_dp,n_dp,nterm,nlev,               &
                             nlcl_uv,ntop_uv,cu_term,cu_tend,cu_tend, &
                             pstar,uw0,vw0,zlcl_uv,ue_p,ve_p,visc,    &
                             flxkp12,p_uv,phalf_uv,rho_uv,timestep,   &
                             ! Input/output
                             uw,vw,                                   &
                             ! Output
                             uw_base,vw_base,uw_deep,vw_deep)


      CALL deep_cmt_incr(n_dp,n_dp,n_dp,nlev,nterm,                  &
                         nlcl_uv,ntop_uv,cu_term,cu_tend,            &
                         zlcl_uv,phalf_uv,p_uv,rho_uv,               &
                         uw_base,vw_base,uw,vw,                      &
                         ! Output
                         dubydt,dvbydt)

    END IF  ! nterm > 0

  CASE (3,4)        ! New Turbulence scheme using heights

    nterm = 0   ! count of number of deep points which actually convected
    DO i = 1,n_dp
      IF (kterm(i)  >=  nlcl_uv(i)) THEN
        nterm = nterm + 1
        cu_term(nterm) = i

        ! Use CAPE scaled mass flux as initial mass flux rather than CRM derived value
        ! to be consistent with thermodynamic part of convection.
        ! mb now completely reset, units now Pa/s
        mb(i) = flxkp12(nlcl_uv(i),i)
        zlcl(i) = z_rho(i,ntml(i))

        ! Cloud velocity scale - derived from CRM simulations
        !             wcld = (C_mass*wstar*CAPE)**(1/3)

        wcld(i) = (delthvu(i) * c_mass * wstar(i) * g / (th(i,ntml(i)) &
              * (1.0 + c_virtual * q(i,ntml(i)))))**0.3333
      END IF
    END DO

    IF (nterm > 0) THEN
      CALL deep_turb_cmt (n_dp, nterm, nlev, deep_cmt_opt,         &
                    ntml, kterm,cu_term,freeze_lev,                &
                    timestep,                                      &
                    uw0, vw0, mb, wcld, wstar ,zlcl,               &
                    flx,                                           &
                    r_rho, r_theta, z_rho, z_theta,rho,rho_theta,  &
                    r2rho, r2rho_th, dr_across_th, dr_across_rh,   &
                    u, v,                                          &
                    dubydt, dvbydt, uw_deep, vw_deep)

    END IF           ! nterm > 0

  END SELECT           ! deep_cmt_opt

END IF  ! L_mom

!-----------------------------------------------------------------------
! 7.1  Add the dissipative heating from the CMT to the theta increment
!-----------------------------------------------------------------------
IF (l_mom .AND. l_cmt_heating) THEN
  CALL cmt_heating(n_dp, nlev,                                      &
                   z_theta, z_rho, exner_layer_centres,             &
                   u, v, dubydt, dvbydt,                            &
                   ! Out
                   dthbydt)
END IF


!-----------------------------------------------------------------------
! 8.0  Energy (and optionally water) correction calculation
!-----------------------------------------------------------------------
DO i = 1,n_dp
  index1(i) = i
END DO

IF (l_cv_conserve_check) THEN
  CALL cor_engy_6a(n_dp, n_dp, nlev, index1,                        &
                   p_layer_boundaries, exner_layer_centres,         &
                   exner_rho, r_theta, r_rho, rho,                  &
                   r2rho_th, r2rho,                                 &
                   dr_across_th, dr_across_rh,                      &
                   dubydt, dvbydt, dqclbydt, dqcfbydt,              &
                   rain, snow, th, q, qcl, qcf, u, v,               &
                   !In/Out
                   dthbydt, dqbydt)
END IF

!-----------------------------------------------------------------------
! 9.0  Correct negative/very small humidities
!-----------------------------------------------------------------------

! only check columns where convection has occurred.

DO i = 1,n_dp
  qminincolumn(i) = q(i,nlev)
END DO
DO k = 1,nlev-1
  DO i = 1,n_dp
    IF (q(i,k)  <   qminincolumn(i)) THEN
      qminincolumn(i) = q(i,k)
    END IF
  END DO
END DO

! Ensure Q does not go below global allowed minimum (QMIN)

DO i = 1,n_dp
  qminincolumn(i)=MAX(qmin,qminincolumn(i))
END DO

! Apply an artificial upwards flux from k-1 level to ensure Q
! remains above minimum value in the column.

DO k = nlev,2,-1
  DO i = 1,n_dp
    IF (dqbydt(i,k) /= 0.0) THEN
      temp1(i)=q(i,k) + dqbydt(i,k) * timestep
      IF (temp1(i)  <   qminincolumn(i)) THEN

        dqbydt(i,k-1) = dqbydt(i,k-1) -                         &
            ((qminincolumn(i) - q(i,k)) / timestep-dqbydt(i,k)) &
             * (r2rho_th(i,k)*dr_across_th(i,k))                &
             / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

        dqbydt(i,k) = (qminincolumn(i) - q(i,k)) / timestep
      END IF
    END IF
  END DO ! n_dp loop
END DO  ! nlev

! check negative q

k=1
DO i = 1,n_dp
  temp1(i)=q(i,k) + dqbydt(i,k) * timestep
  IF (temp1(i)  <   qminincolumn(i) .AND.                                    &
      printstatus >= prstatus_normal) THEN

    warning = -1
    WRITE(umMessage,'(A)')                                                   &
      '  PROBLEM: negative q after deep convection.'
    CALL umPrint(umMessage,src='deep_conv_6a')
    WRITE(umMessage,'(2(a4,i6),a9,g26.18,a7,g26.18)')                        &
      '  i:',i,' k: ',k, ' q after ',temp1(i),' dq/dt ',dqbydt(i,k)
    CALL umPrint(umMessage,src='deep_conv_6a')

    WRITE(umMessage,'(A6,5A16)') 'k','q inc','qcl inc','qcf inc','mf','dp/dt'
    CALL umPrint(umMessage,src='deep_conv_6a')

    DO k=1,nlev
      WRITE(umMessage,'(I6,5G16.8)') k, dqbydt(i,k), dqclbydt(i,k),          &
                                     dqcfbydt(i,k), flx(i,k),                &
                (p_layer_boundaries(i,k)-p_layer_boundaries(i,k-1))/timestep
      CALL umPrint(umMessage,src='deep_conv_6a')
    END DO

    CALL ereport( 'deep_conv_6a', warning,                               &
            "Negative q after deep convection."               //newline//&
            "Check the profiles at the grid-point where this" //newline//&
            "has occured (these are printed in the run"       //newline//&
            "output file)." )

  END IF
END DO ! n_dp loop

!-----------------------------------------------------------------------
! 10.0  3D - Convective cloud amount assumed 3d required ie L_3d_cca
!      is true in old code
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------------
!      CCRad - Calculate CCA for Deep levels only
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------
! 10.1 Calculate CCA_2D of Deep Cloud
!-----------------------------------------------------------------
SELECT CASE (cca2d_dp_opt)
CASE (srf_precip)
  DO i=1, n_dp
    IF (iccb(i) > 0) THEN ! Deep convection was successful

      ! Use determination of CCA_2D based on surface precip. rate
      ! from deep cloud.

      ! NOTE: at present the a_ and b_ parameters for LAND and SEA
      !       are equal, so there will be no difference between
      !       land and sea points.

      IF ((rain(i) + snow(i)) > 0.0) THEN
        IF (bland(i)) THEN
          ! Land point
          tempnum = a_land + b_land                                 &
                  * LOG(rsec_per_day * (rain(i)+snow(i)))
        ELSE
          ! Sea point

          tempnum = a_sea + b_sea                                   &
                  * LOG(rsec_per_day * (rain(i)+snow(i)))
        END IF


        cca_2d(i) = MAX(2.0e-5, tempnum)

        ! Grab lowest cca value before any tuning occurs
        ! This will overwrite lcca in ni_conv_ctl only if neither
        ! shallow or deep have occured.  This is under a switch in
        ! the 4a scheme.
        !
        ! NOTE: Downdraughts & Evaporation still being fed cca_2d
        !       derived from TCW, This issue may require further
        !       investigation.
        lcca(i) = cca_2d(i)

      END IF

    END IF      ! iccb
  END DO      ! i (n_dp)

CASE (total_condensed_water)
  ! cca_2d_dp left unchanged from code, which is based on
  ! TCW (Total Condensed Water) (This is a rate)

END SELECT

! l_dcpl_cld4pc2 is set to true for 5a scheme, so
! CCRad tuning knobs applied in glue_conv

!---------------------------------------------------------------------
! 10.2 Apply CCA_2D to 3d cloud profile
!---------------------------------------------------------------------

IF (l_anvil) THEN

  ! Apply anvil scheme to deep cloud
  CALL calc_3d_cca                                                    &
    ( n_dp, n_dp, nlev, n_cca_lev, nbl, iccb, icct                    &
    , p_layer_boundaries, freeze_lev, cca_2d, cca, z_theta, z_rho )

  ! NOTE: iccb, icct are layer centres (theta levels) at this
  !        point.

ELSE

  ! Apply cca_2d to all levels from deep base to deep top
  DO i=1, n_dp
    ! Only copy across to cca if both iccb(i) and icct(i) have
    ! been set and are within cca array bounds
    IF ((iccb(i) > 0) .AND. (icct(i) > 0) .AND.                       &
        (iccb(i) <= n_cca_lev)) THEN
      DO k=iccb(i), MIN(icct(i), n_cca_lev)
        cca(i,k) = cca_2d(i)
      END DO
    END IF
  END DO
END IF      ! l_anvil

! l_dcpl_cld4pc2 is set to true for 5a scheme, so
! CCRad tuning knobs applied in glue_conv


!-----------------------------------------------------------------------
! Final SCM convection sub-step diagnostics
!-----------------------------------------------------------------------
IF ( l_scm_convss_dg ) THEN

  ! 3-D diagnostics
  DO k = 1, nlev
    DO i = 1, n_dp
      IF ( k>=ntml(i) .AND. k<=det_lev(i) ) THEN
        ! Save the final updraft mass-flux
        scm_convss_dg(i) % up_flx(k) = flx(i,k)
      END IF
    END DO
  END DO

  ! 2-D diagnostics
  DO i = 1, n_dp
    scm_convss_dg(i) % precip_deep = rain(i) + snow(i)
  END DO

END IF


!-----------------------------------------------------------------------
! 11.0  End Subroutine
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE deep_conv_6a
END MODULE deep_conv_6a_mod
