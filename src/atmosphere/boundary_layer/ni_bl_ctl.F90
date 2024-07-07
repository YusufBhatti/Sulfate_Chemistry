! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! purpose: Interface to boundary layer mixing coefficients calculation
!   language: fortran 95
!   this code is written to umdp3 programming standards.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary layer
!---------------------------------------------------------------------
MODULE ni_bl_ctl_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'NI_BL_CTL_MOD'
CONTAINS

SUBROUTINE ni_bl_ctl (                                                  &
! IN parameters for iterative SISL scheme
 cycleno,l_jules_call,                                                  &
! IN time stepping information
 val_year, val_day_number, val_hour, val_minute, val_second,            &
! IN model dimensions.
 land_points, ntiles, bl_levels,                                        &
! IN switches
 l_scrn, l_aero_classic,                                                &
! IN data fields.
 p, p_layer_centres, rho_wet_rsq, rho_wet, rho_dry, u_p, v_p,           &
 u_px, v_px, u_0_px, v_0_px,                                            &
 land_sea_mask, q, qcl, qcf, p_star, theta, exner_theta_levels, rad_hr, &
 micro_tends, soil_layer_moisture, rho_wet_tq, z_uv, z_tq,              &
! IN ancillary fields and fields needed to be kept from tstep to tstep
 hcon, smvccl, smvcwt, smvcst, sthf, sthu, sil_orog_land,               &
 ho2r2_orog, sd_orog, ice_fract_cat, k_sice,                            &
 land_index, photosynth_act_rad,                                        &
! IN variables required for mineral dust scheme
 soil_clay,soil_sand,dust_mrel1,dust_mrel2,dust_mrel3,                  &
 dust_mrel4,dust_mrel5,dust_mrel6,                                      &
! IN additional variables for JULES
 canopy ,catch, catch_snow, snow_tile, z0_tile, z0h_tile_bare, z0m_soil,&
 lw_down,sw_tile,tstar_tile,tsurf_elev_surft,                           &
 co2_3d,asteps_since_triffid,                                           &
 cs,frac,canht_ft,lai_ft,fland,flandg,albsoil,cos_zenith_angle,         &
! IN everything not covered so far
 t_soil, ti, ti_cat, tstar, zh_prev,ddmfx, cf_bulk, zhpar, zlcl,        &
! IN SCM namelist data
 l_spec_z0, z0m_scm, z0h_scm,flux_e, flux_h, ustar_in,                  &
! SCM Diagnostics (dummy values in full UM) and STASH
 nSCMDpkgs,L_SCMDiags, BL_diag, sf_diag,                                &
! INOUT data
 Gs,z0msea,w,etadot,tstar_sea,tstar_sice_cat,zh,dzh,cumulus,ntml,ntpar, &
 l_shallow,error_code,                                                  &
! INOUT additional variables for JULES
 g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,                        &
! INOUT variables for TKE based turbulence schemes
 e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                          &
! INOUT variables from bdy_expl1 needed elsewhere
 bq_gb, bt_gb, dtrdz_charney_grid, rdz_charney_grid,                    &
 dtrdz_u, dtrdz_v, rdz_u, rdz_v, k_blend_tq, k_blend_uv,                &
! INOUT variables from Jules needed elsewhere
 flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m, cdr10m_n, cd10m_n,   &
 fqw, ftl, rib_gb, vshr, z0m_eff_gb, z0h_eff_gb, r_b_dust,              &
 rho_aresist,aresist,resist_b, rhokm,rhokh,                             &
! INOUT diagnostics required for soil moisture nudging scheme :
 wt_ext,                                                                &
! INOUT variables required in IMP_SOLVER
 alpha1_sea, alpha1_sice, ashtf_sea, ashtf,  uStarGBM,                  &
! INOUT additional variables for JULES
 ftl_tile,radnet_sea,radnet_sice,rib_tile,rho_aresist_tile,             &
 aresist_tile,resist_b_tile,alpha1,ashtf_tile,fqw_tile,epot_tile,       &
 fqw_ice,ftl_ice,fraca,resfs,resft,rhokh_tile,rhokh_sice,rhokh_sea,     &
 z0hssi,z0h_tile,z0m_gb,z0mssi,z0m_tile,chr1p5m,chr1p5m_sice,smc,       &
 gpp,npp,resp_p,g_leaf,gpp_ft,npp_ft,resp_p_ft,resp_s,resp_s_tot,       &
 resp_w_ft,gc,canhc_tile,wt_ext_tile,flake,tile_index,tile_pts,         &
 tile_frac,fsmc,vshr_land,vshr_ssi,tstar_land,tstar_ssi,dtstar_tile,    &
 dtstar_sea,dtstar_sice,hcons,emis_tile,emis_soil, t1_sd, q1_sd,fb_surf,&
! OUT variables for message passing
 tau_fd_x, tau_fd_y, rhogamu, rhogamv, f_ngstress,                      &
! OUT diagnostics (done after implicit solver)
 zht, zhnl, shallowc, cu_over_orog, bl_type_1,bl_type_2,bl_type_3,      &
 bl_type_4,bl_type_5,bl_type_6,bl_type_7, bl_w_var,                     &
! OUT data required for tracer mixing & dust:
 dust_flux,dust_emiss_frac,                                             &
 u_s_t_tile,u_s_t_dry_tile,u_s_std_tile,kent, we_lim, t_frac, zrzi,     &
 kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                      &
! OUT fields
 nbdsc, ntdsc,wstar, wthvs,uw0,vw0,taux_p,tauy_p,rhcpt, rib_ssi         &
 )

USE atm_fields_bounds_mod, ONLY:                                        &
                                udims, vdims, tdims, tdims_s, tdims_l,  &
                                pdims, pdims_s, wdims
USE atm_step_local, ONLY: dim_cs1, dim_cs2, land_pts_trif, npft_trif,   &
       co2_dim_len,co2_dim_row
USE bdy_layr_mod, ONLY: bdy_layr
USE bl_diags_mod, ONLY: strnewbldiag
USE bl_option_mod, ONLY: l_quick_ap2
USE c_elevate, ONLY: l_elev_absolute_height
USE cloud_inputs_mod, ONLY: i_rhcpt
USE cv_run_mod, ONLY: l_jules_flux
USE dynamics_input_mod, ONLY: numcycles
USE dust_parameters_mod, ONLY: ndiv, ndivh
USE jules_sea_seaice_mod, ONLY: l_ctile, nice, nice_use
USE jules_soil_mod, ONLY: sm_levels
USE jules_surface_types_mod, ONLY: ntype, npft
USE level_heights_mod, ONLY:   r_theta_levels
USE model_domain_mod, ONLY: model_type, mt_single_column
USE mphys_inputs_mod, ONLY: l_subgrid_qcl_mp
USE pc2_constants_mod, ONLY: rhcpt_tke_based
USE planet_constants_mod, ONLY: planet_radius
USE sf_diags_mod, ONLY: strnewsfdiag
USE stash_array_mod, ONLY: sf
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE
!---------------------------------------------------------------------
! arguments with intent in. ie: input variables.
!---------------------------------------------------------------------
INTEGER, INTENT(IN) ::                                                  &
 cycleno    ! Iteration no

INTEGER, INTENT(IN) ::                                                  &
 !time information for current timestep
 val_year, val_day_number, val_hour, val_minute, val_second

LOGICAL, INTENT(IN) ::                                                  &
  l_aero_classic,                                                       &
       !switch for CLASSIC aerosol scheme
  l_jules_call  ! switch for whether this call is for jules only or the
                ! main scheme

! Model dimensions
INTEGER, INTENT(IN) ::                                                  &
  land_points,                                                          &
            ! IN No.of land points being processed, can be 0.
  ntiles,                                                               &
            ! IN No. of land-surface tiles ( JULES )
  bl_levels

LOGICAL, INTENT(IN) ::                                                  &
  l_scrn
                         ! Logical to control output
                         !    of screen level T,Q,QCL,QCF

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
 nSCMDpkgs                ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
 L_SCMDiags(nSCMDpkgs)    ! Logicals for SCM diagnostics packages

! Data arrays
REAL, INTENT(IN) ::                                                     &
  u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      pdims%k_start:pdims%k_end),                                       &
  v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,              &
      pdims%k_start:pdims%k_end),                                       &
  rho_wet_rsq(pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end),                           &
!                       ! wet density times r^2 on rho levels (kg/m3)
    rho_wet(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            1:tdims%k_end),                                             &
!                       ! wet density on rho levels (kg/m3)
    rho_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            pdims%k_start:pdims%k_end),                                 &
!                       ! dry density on rho levels (kg/m3)
    z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
           bl_levels+1),                                                &
                              ! Z_uv(*,K) is height of half
!                                   ! level k-1/2.
    z_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
           bl_levels),                                                  &
    p(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,      &
      pdims_s%k_start:pdims_s%k_end),                                   &
    p_layer_centres(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end,0:tdims%k_end),           &
    p_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
    theta(tdims_s%i_start:tdims_s%i_end,                                &
          tdims_s%j_start:tdims_s%j_end,                                &
          tdims_s%k_start:tdims_s%k_end),                               &
    exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                       tdims_s%j_start:tdims_s%j_end,                   &
                       tdims_s%k_start:tdims_s%k_end),                  &
    q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,      &
      tdims_l%k_start:tdims_l%k_end),                                   &
    qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,    &
        tdims_l%k_start:tdims_l%k_end),                                 &
    qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,    &
        tdims_l%k_start:tdims_l%k_end),                                 &
    rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
                          2,bl_levels),                                 &
                        ! IN (LW,SW) radiative heating rate (K/s)
    micro_tends(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                       2, bl_levels),                                   &
!                          ! Tendencies from microphys within BL levels
!                          ! (TL, K/s; QW, kg/kg/s)
    u_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,   &
         bl_levels),                                                    &
    v_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,   &
         bl_levels),                                                    &
    u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),&
    v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

!     Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag
REAL, INTENT(IN) ::                                                     &
 soil_layer_moisture(land_points,sm_levels)!IN soil moisture
!                 ! per layer (kg m-2)
LOGICAL, INTENT(IN) ::                                                  &
  land_sea_mask(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

! ancillary arrays and fields required to be saved from tstep to tstep
INTEGER, INTENT(IN) ::                                                  &
  land_index (land_points),                                             &
                           ! set from land_sea_mask
 asteps_since_triffid

REAL, INTENT(IN) ::                                                     &
  k_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         nice_use),                                                     &
               ! sea ice surface layer effective conductivity (W/m2/K)
  hcon (land_points),                                                   &
                     ! soil/qrparm.soil.hcond
  smvccl (land_points,sm_levels),                                       &
                     ! soil/qrparm.soil.crit
  smvcwt (land_points,sm_levels),                                       &
                     ! soil/qrparm.soil.wilt
  smvcst (land_points,sm_levels),                                       &
                     ! soil/qrparm.soil.satn
  sthf(land_points,sm_levels),                                          &
                        ! IN Frozen soil moisture content of
                        !    each layer as a fraction of
                        !    saturation.
  sthu(land_points,sm_levels),                                          &
                        ! IN Unfrozen soil moisture content
                        !    of each layer as a fraction of
                        !    saturation.
  sil_orog_land (land_points),                                          &
                         ! orog/qrparm.orog.as
  ho2r2_orog (land_points),                                             &
                         ! orog/qrparm.orog.h2root2
  sd_orog (land_points),                                                &
                         ! orog/qrparm.orog.stdev
  ice_fract_cat (pdims%i_start:pdims%i_end,                             &
                 pdims%j_start:pdims%j_end,nice_use),                   &
                    ! ice/qrclim.ice.(month)
                    ! If nice_use=1, this is the sum of the categories
  photosynth_act_rad(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end)
                                     ! Net downward
!                                 shortwave radiation in band 1 (w/m2).

REAL, INTENT(IN) ::                                                     &
   ! mineral dust fields
  soil_clay (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  soil_sand (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  dust_mrel1 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel2 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel3 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel4 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel5 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel6 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

REAL, INTENT(IN) ::                                                     &
 rho_wet_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
!                               ! RHO_WET_TQ(*,K) is the density at half
!                               ! level k+1/2.
    z0m_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                           ! Fixed Sea surface roughness
                           ! length(m) for momentum (SCM)
    z0h_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                           ! Fixed Sea surface roughness
                           ! length(m) for heat (SCM)
    t_soil(land_points,sm_levels),                                      &
                                 ! slt/qrclim.slt_pm(lev).(month)
    ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                        ! set equal to tstar

    ti_cat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),   &
    cf_bulk(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,tdims%k_end),                     &
    zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! IN boundary layer height from
!                                 !    previous timestep
    ddmfx(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
!                                 ! IN Convective downdraught
!                                 !    mass-flux at cloud base
    flux_e(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
    flux_h(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
    ustar_in(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
    zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                               ! height of ntpar
    zlcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! height of lcl accurate value not
                               ! a model level height (m)
LOGICAL, INTENT(IN) ::                                                  &
 l_spec_z0                  ! T if roughness lengths have been specified

REAL, INTENT(IN) ::                                                     &
 canopy(land_points,ntiles),                                            &
                              ! IN Surface/canopy water for
!                                  !    snow-free land tiles (kg/m2)
   catch(land_points,ntiles),                                           &
                                ! IN Surface/canopy water capacity
!                                  !    of snow-free land tiles (kg/m2).
   catch_snow(land_points,ntiles),                                      &
                             ! IN Snow interception capacity of
!                                  !    tiles (kg/m2).
   snow_tile(land_points,ntiles),                                       &
                                ! IN Lying snow on tiles (kg/m2)
   z0_tile(land_points,ntiles),                                         &
                                ! IN Tile roughness lengths (m).
   z0h_tile_bare(land_points,ntiles),                                   &
                                ! IN Tile thermal roughness lengths (m)
                                ! without snow
   z0m_soil(land_points),                                               &
                                ! IN bare soil momentum roughness length (m)
   lw_down(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                             ! IN Surface downward LW radiation
!                                  !    (W/m2).
   sw_tile(land_points,ntiles),                                         &
                                ! IN Surface net SW radiation on
!                                  !    land tiles (W/m2).
   tstar_tile(land_points,ntiles),                                      &
                                ! IN Surface tile temperatures
   tsurf_elev_surft(land_points,ntiles),                                &
                                ! IN
   co2_3d(co2_dim_len,co2_dim_row),                                     &
!                                  ! IN 3D CO2 field if required.
   cs(land_points,dim_cs1),                                             &
                           ! IN Soil carbon (kg C/m2).
   frac(land_points,ntype),                                             &
                                ! IN Fractions of surface types.
   canht_ft(land_points,npft),                                          &
                                ! IN Canopy height (m)
   lai_ft(land_points,npft),                                            &
                                ! IN Leaf area index
   fland(land_points),                                                  &
                             ! IN Land fraction on land points.
   flandg(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end),                               &
                             ! IN Land fraction on all points.
   albsoil(land_points),                                                &
                             ! Soil albedo.
   cos_zenith_angle(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end)
!                                  ! Cosine of the zenith angle
!---------------------------------------------------------------------
! arguments with intent INOUT. ie: input variables changed on output.
!---------------------------------------------------------------------
REAL, INTENT(INOUT) ::                                                  &
  z0msea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                         ! Sea surface roughness length(m)
                         ! for momentum
  zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
  dzh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                         ! INOUT inversion thickness
  w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,                &
    0:wdims%k_end),                                                     &
  etadot(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,           &
         0:wdims%k_end),                                                &
 tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
 tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                           ! INOUT Open sea sfc temperature (K).
 tstar_sice_cat(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,nice_use),                    &
                           ! INOUT Sea-ice sfc temperature (K).
 Gs(land_points),                                                       &
                              ! INOUT "Stomatal" conductance to
!                                  !        evaporation (m/s).
   g_leaf_acc(land_points,npft),                                        &
                                ! INOUT Accumulated G_LEAF
   npp_ft_acc(land_pts_trif,npft_trif),                                 &
!                                  ! INOUT Accumulated NPP_FT
   resp_w_ft_acc(land_pts_trif,npft_trif),                              &
!                                  ! INOUT Accum RESP_W_FT
   resp_s_acc(land_pts_trif,dim_cs2)
                                   ! INOUT Accumulated RESP_S

! INOUT variables for TKE based turbulence schemes
REAL, INTENT(INOUT) ::                                                  &
  e_trb(tdims%i_start:tdims%i_end,                                      &
        tdims%j_start:tdims%j_end,bl_levels),                           &
!                   ! TKE defined on theta levels K-1
    tsq_trb(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,bl_levels),                       &
!                   ! Self covariance of liquid potential temperature
!                   ! (thetal'**2) defined on theta levels K-1
    qsq_trb(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,bl_levels),                       &
!                   ! Self covariance of total water
!                   ! (qw'**2) defined on theta levels K-1
    cov_trb(tdims%i_start:tdims%i_end,                                  &
            tdims%j_start:tdims%j_end,bl_levels),                       &
!                   ! Correlation between thetal and qw
!                   ! (thetal'qw') defined on theta levels K-1
    zhpar_shcu(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!                   ! Height of mixed layer used to evaluate
!                   ! the non-gradient buoyancy flux

LOGICAL, INTENT(INOUT) ::                                               &
  cumulus (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                           ! *APL bl convection flag
  l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!                            ! Logical indicator of shallow convection

INTEGER, INTENT(INOUT) ::                                               &
  ntml (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                               ! INOUT Number of model layers in the
                               !      turbulently mixed layer
  ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! IN/OUT top of diagnostic parcel
!                                !        ascent
INTEGER, INTENT(INOUT) ::                                               &
  error_code

! Variables which are OUT but need to be declared as INOUT because this
! routine is called twice, i.e. they are OUT from its first call, but
! are either needed for the 2nd call or elsewhere in the UM after the
! 2nd call

! ... ones from bdy_expl1
REAL, INTENT(INOUT) ::                                                  &
  dtrdz_charney_grid(pdims%i_start:pdims%i_end,                         &
                     pdims%j_start:pdims%j_end,bl_levels),              &
  rdz_charney_grid(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,bl_levels),                &
  dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,          &
                          bl_levels),                                   &
  dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,          &
                            bl_levels),                                 &
  bq_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
  bt_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
  rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
                        2:bl_levels),                                   &
  rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
                          2:bl_levels)

INTEGER, INTENT(INOUT) ::                                               &
 k_blend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
 k_blend_uv(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

! ... ones from surf_couple_explicit
REAL, INTENT(INOUT) ::                                                  &
  uStarGBM(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
!       ! GBM surface friction velocity
  t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                        ! set to zero initially
  q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                        ! set to zero initially
  fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                        ! Surface flux buoyancy over density (m^2/s^3)
  z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  z0h_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                        ! Effective grid-box roughness
!                                 lengths for momentum and for
!                                 heat, moisture
  rhokm(pdims_s%i_start:pdims_s%i_end,                                  &
        pdims_s%j_start:pdims_s%j_end,bl_levels),                       &
  alpha1_sice(pdims%i_start:pdims%i_end,                                &
              pdims%j_start:pdims%j_end,nice_use),                      &
  ashtf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),  &
  alpha1_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  ashtf_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
  flandfac(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end),                              &
  fseafac(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end),                               &
  rhokm_land(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end),                            &
  rhokm_ssi(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end),                             &
  cdr10m(pdims_s%i_start:pdims_s%i_end,                                 &
         pdims_s%j_start:pdims_s%j_end),                                &
  cdr10m_n(pdims_s%i_start:pdims_s%i_end,                               &
         pdims_s%j_start:pdims_s%j_end),                                &
  cd10m_n(pdims_s%i_start:pdims_s%i_end,                                &
         pdims_s%j_start:pdims_s%j_end),                                &
  rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                        ! INOUT RHOSTAR*CD_STD*VSHR
                        !     for CLASSIC aerosol scheme
  aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                        ! INOUT 1/(CD_STD*VSHR)
                        !     for CLASSIC aerosol scheme
  resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                        ! INOUT (1/CH-1/(CD_STD)/VSHR
                        !     for CLASSIC aerosol scheme
  r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           ndiv),                                                       &
                              !INOUT surface layer resist for dust
    vshr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
    ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
                                       ! needed as diagnostic
    fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
                                       ! needed as diagnostic ?
    rib_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                          ! Mean bulk Richardson number for
                               !  lowest layer.
    rhokh (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
    wt_ext(land_points,sm_levels),                                      &
                                !INOUT cumulative fract of transp'n
   tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                             ! INOUT Land mean sfc temperature (K)

   tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),      &
                             ! INOUT Sea mean sfc temperature (K).
 ftl_tile(land_points,ntiles),                                          &
                              ! INOUT Surface FTL for land tiles
 radnet_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                           ! INOUT Surface net radiation on open sea (W/m2)
 radnet_sice(pdims%i_start:pdims%i_end,                                 &
             pdims%j_start:pdims%j_end,nice_use),                       &
                           ! INOUT Surface net radiation on sea-ice (W/m2)
   rib_tile(land_points,ntiles),                                        &
                                ! INOUT RIB for land tiles.
   rho_aresist_tile(land_points,ntiles),                                &
!                                  ! INOUT RHOSTAR*CD_STD*VSHR on land
!                                  !     tiles for CLASSIC aerosol scheme
   aresist_tile(land_points,ntiles),                                    &
!                                  ! INOUT 1/(CD_STD*VSHR) on land tiles
                                   !     for CLASSIC aerosol scheme
   resist_b_tile(land_points,ntiles),                                   &
!                                  ! INOUT (1/CH-1/CD_STD)/VSHR on land
!                                  !     tiles for CLASSIC aerosol scheme
   alpha1(land_points,ntiles),                                          &
                                ! INOUT Mean gradient of saturated
!                                  !     specific humidity with respect
!                                  !     to temperature between the
!                                  !     bottom model layer and tile
!                                  !     surfaces
   ashtf_tile(land_points,ntiles),                                      &
                                !INOUT Coefficient to calculate
!                                  !     surface heat flux into land
!                                  !     tiles.
   fqw_tile(land_points,ntiles),                                        &
                                ! INOUT Surface FQW for land tiles
   epot_tile(land_points,ntiles),                                       &
                                ! INOUT Local EPOT for land tiles.
   fqw_ice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           nice_use),                                                   &
                             ! INOUT Surface FQW for sea-ice
   ftl_ice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           nice_use),                                                   &
                             ! INOUT Surface FTL for sea-ice
   fraca(land_points,ntiles),                                           &
                                ! INOUT Fraction of surface moisture
!                                  !     flux with only aerodynamic
!                                  !     resistance for snow-free land
!                                  !     tiles.
   resfs(land_points,ntiles),                                           &
                                ! INOUT Combined soil, stomatal
!                                  !     and aerodynamic resistance
!                                  !     factor for fraction (1-FRACA)
!                                  !     of snow-free land tiles.
   resft(land_points,ntiles),                                           &
                                ! INOUT Total resistance factor.
!                                  !     FRACA+(1-FRACA)*RESFS for
!                                  !     snow-free land, 1 for snow.
   rhokh_tile(land_points,ntiles),                                      &
                                ! INOUT Surface exchange coefficient
!                                  !     for land tiles
   rhokh_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
                                                             nice_use), &
                             ! INOUT Surface exchange coefficients
!                                  !     for sea-ice.
   rhokh_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                             ! INOUT Surface exchange coefficients
!                                  !     for sea
   z0hssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   z0mssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                             ! INOUT Roughness lengths over sea (m).
   z0h_tile(land_points,ntiles),                                        &
                             ! INOUT Tile roughness lengths for h
!                                  !     and moisture (m).
   z0m_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                             ! INOUT Grid-box mean roughness length
!                                  !      for momentum (m).
   z0m_tile(land_points,ntiles),                                        &
                             ! INOUT Tile roughness lengths for
!                                  !     momentum.
   chr1p5m(land_points,ntiles),                                         &
                             ! INOUT Ratio of coefffs for
!                                  !     calculation of 1.5m temp for
!                                  !     land tiles.
   chr1p5m_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
!                                  ! INOUT CHR1P5M for sea and sea-ice
!                                  !     (leads ignored).
   smc(land_points),                                                    &
                                ! INOUT Available moisture in the
!                                  !     soil profile (mm).
   gpp(land_points),                                                    &
                                ! INOUT Gross primary productivity
!                                  !     (kg C/m2/s).
   npp(land_points),                                                    &
                                ! INOUT Net primary productivity
!                                  !     (kg C/m2/s).
   resp_p(land_points),                                                 &
                                ! INOUT Plant respiration (kg C/m2/s
   g_leaf(land_points,npft),                                            &
                                ! INOUT Leaf turnover rate (/360days
   gpp_ft(land_points,npft),                                            &
                                ! INOUT Gross primary productivity
!                                  !     on PFTs (kg C/m2/s).
   npp_ft(land_points,npft),                                            &
                                ! INOUT Net primary productivity
!                                  !     (kg C/m2/s).
   resp_p_ft(land_points,npft),                                         &
                                ! INOUT Plant respiration on PFTs
!                                  !     (kg C/m2/s).
   resp_s(land_points,dim_cs1),                                         &
                            ! INOUT Soil respiration (kg C/m2/s)
   resp_s_tot(dim_cs2),                                                 &
                           ! INOUT total soil respiration
   resp_w_ft(land_points,npft),                                         &
                                ! INOUT Wood maintenance respiration
!                                  !     (kg C/m2/s).
   gc(land_points,ntiles),                                              &
                                ! INOUT "Stomatal" conductance to
!                                  !      evaporation for land tiles
!                                  !      (m/s).
   canhc_tile(land_points,ntiles),                                      &
                                ! INOUT Areal heat capacity of canopy
!                                  !    for land tiles (J/K/m2).
   wt_ext_tile(land_points,sm_levels,ntiles),                           &
!                                  ! INOUT Fraction of evapotranspiration
!                                  !    which is extracted from each
!                                  !    soil layer by each tile.
   flake(land_points,ntiles),                                           &
                                ! INOUT Lake fraction.
   tile_frac(land_points,ntiles),                                       &
                                ! INOUT Tile fractions including
!                                  !     snow cover in the ice tile.
   fsmc(land_points,npft),                                              &
                                ! INOUT Moisture availability
!                                  !     factor.
    vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
    vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
 dtstar_tile(land_points,ntiles),                                       &
                           ! Change in TSTAR over timestep
!                                  ! for land tiles
   dtstar_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                             ! Change is TSTAR over timestep
                             ! for open sea
   dtstar_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),&
                             ! Change is TSTAR over timestep
                             ! for sea-ice
   hcons(land_points),                                                  &
                             ! Soil thermal conductivity
                             ! including water and ice
   emis_tile(land_points,ntiles),                                       &
                             ! Emissivity for land tiles
   emis_soil(land_points)
                             ! Emissivity of underlying soil


INTEGER, INTENT(INOUT) ::                                               &
 tile_index(land_points,ntype),                                         &
                              ! INOUT Index of tile points
 tile_pts(ntype)             ! INOUT Number of tile points

!---------------------------------------------------------------------
! arguments with intent OUT. ie: output variables.
!---------------------------------------------------------------------
REAL, INTENT(OUT) ::                                                    &
    uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
!                       ! U-component of surface wind stress (P-grid)
    vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
!                       ! V-component of surface wind stress (P-grid)
    wstar(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end),         &
                               ! surface-based mixed layer
!                                    ! velocity scale
    wthvs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                               ! surface buoyancy flux
REAL, INTENT(OUT) :: taux_p( pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,                 &
                             bl_levels )
REAL, INTENT(OUT) :: tauy_p( pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,                 &
                             bl_levels )
         ! Wind stresses on theta-levels, on p-grid (for convection)

INTEGER, INTENT(OUT) ::                                                 &
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
 kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                            ! OUT grid-level of SML inversion
 kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! OUT grid-level of DSC inversion

REAL, INTENT(OUT) :: rhcpt(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,1:tdims%k_end)
! variables passed from BDY_LAYR to IMP_SOLVER
REAL, INTENT(OUT) ::                                                    &
  tau_fd_x(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end,bl_levels),                    &
  tau_fd_y(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end,bl_levels),                    &
  rhogamu(pdims_s%i_start:pdims_s%i_end,                                &
           pdims_s%j_start:pdims_s%j_end ,bl_levels),                   &
  rhogamv(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
  f_ngstress(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,2:bl_levels)

! Diagnostics needed in NI_imp_ctl
REAL, INTENT(OUT) ::                                                    &
  dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            ndiv),                                                      &
                               !OUT dust emissions (kg m-2 s-1)
  dust_emiss_frac(land_points,ntiles),                                  &
                        ! OUT fraction of tile can emit dust
  u_s_t_tile(land_points,ntiles,ndivh),                                 &
                                   !OUT threshold frict. vel
  u_s_t_dry_tile(land_points,ntiles,ndivh),                             &
                                       !OUT dry soil value
  u_s_std_tile(land_points,ntiles),                                     &
                                !OUT friction velocity
  we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                            ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
    zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                              ! OUT (z-z_base)/(z_i-z_base)
    t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                              ! OUT a fraction of the timestep
    we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
!                                   ! OUT rho*entrainment rate implied b
!                                   !     placing of subsidence
    zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                              ! OUT (z-z_base)/(z_i-z_base)
    t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
!                                   ! OUT a fraction of the timestep
    zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! OUT Top of decoupled layer
    zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                ! Max height of turb mixing
    zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                ! non-local BL depth
    shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                             !OUT Indicator set to 1.0 if shallow,
!                                  !   0.0 if not shallow or not cumulus
    cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
!                                  ! OUT Indicator for cumulus
!                                  !     over steep orography
!                                  !   Indicator set to 1.0 if true,
!                                  !   0.0 if false. Exclusive.
    bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT Indicator set to 1.0 if stable
!                                  !     b.l. diagnosed, 0.0 otherwise.
    bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT Indicator set to 1.0 if Sc over
!                                  !     stable surface layer diagnosed,
!                                  !     0.0 otherwise.
    bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT Indicator set to 1.0 if well
!                                  !     mixed b.l. diagnosed,
!                                  !     0.0 otherwise.
    bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer (not over
!                                  !     cumulus) diagnosed,
!                                  !     0.0 otherwise.
    bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT Indicator set to 1.0 if
!                                  !     decoupled Sc layer over cumulus
!                                  !     diagnosed, 0.0 otherwise.
    bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT Indicator set to 1.0 if a
!                                  !     cumulus capped b.l. diagnosed,
!                                  !     0.0 otherwise.
    bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!                                  !     shear-dominated  b.l.
!                                  !     diagnosed, 0.0 otherwise.

REAL, INTENT(OUT) ::                                                    &
   rib_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                             ! OUT Sea mean bulk Richardson number
                             !     for lowest layer.

! Boundary layer TKE diagnostic variable
REAL, INTENT(OUT) :: bl_w_var( tdims%i_start : tdims%i_end,             &
                               tdims%j_start : tdims%j_end,             &
                                           2 : tdims%k_end+1 )
                          ! BL w-variance diagnostic (m/s)^2
                          ! Note: not defined at surface (k=1)

!---------------------------------------------------------------------
! local variables.
!---------------------------------------------------------------------
! loop counters
INTEGER ::                                                              &
  i, j, k ,l

! Diagnostic switches
LOGICAL :: l_apply_diag

REAL ::                                                                 &
  t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)

REAL ::                                                                 &
 z_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
!                       !land height over fractional land points

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='NI_BL_CTL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! error information
IF ( error_code  ==  0) THEN
  ! ----------------------------------------------------------------------
  ! Section BL.0 Initialisation of variables.
  ! ----------------------------------------------------------------------
  IF (l_jules_call .OR. .NOT. l_jules_flux) THEN
    SELECT CASE (model_type)
    CASE DEFAULT
      ! Apply diags at last cycle only or if l_quick_ap2 is .true.
      l_apply_diag = (cycleno == numcycles .OR. l_quick_ap2)

      ! Set diagnostic flags required for boundary layer diagnostics from
      ! STASHflags.
      ! --------------------------------------------------------
      ! Note that an equivalent block of code exists in routine
      ! ni_imp_ctl, and needs to be kept consistent.
      ! --------------------------------------------------------
      ! Windspeed (227, 230) and u, v at 10m on 'B' or 'C' grid
      sf_diag%su10    = (l_apply_diag .AND.                             &
                  ( sf(209,3) .OR. sf(225,3) .OR. sf(227,3) .OR.        &
                    sf(230,3) .OR. sf(463,3) ))

      sf_diag%sv10    = (l_apply_diag .AND.                             &
                  ( sf(210,3) .OR. sf(226,3) .OR. sf(227,3) .OR.        &
                    sf(230,3) .OR. sf(463,3) ))

      sf_diag%sq_t1p5 = (l_apply_diag .AND.                             &
                  ( sf(236,3) .OR. sf(237,3) .OR. sf(245,3) .OR.        &
                    sf(247,3) .OR. sf(248,3) .OR. sf(250,3) .OR.        &
                    sf(341,3) .OR. sf(342,3) .OR. sf(253,3) .OR.        &
                    sf(328,3) .OR. sf(329,3) .OR. l_scrn) )

      sf_diag%sq1p5   = (l_apply_diag .AND. sf_diag%sq_t1p5)
      sf_diag%st1p5   = (l_apply_diag .AND. sf_diag%sq_t1p5)
      sf_diag%l_t10m  = (l_apply_diag .AND. sf(344,3))
      sf_diag%l_q10m  = (l_apply_diag .AND. sf(345,3))
      sf_diag%sfme    = (l_apply_diag .AND. sf(224,3))
      sf_diag%l_cd_ssi = (l_apply_diag .AND. sf(538,3))
      sf_diag%l_ch_ssi = (l_apply_diag .AND. sf(541,3))
      sf_diag%sz0heff = (l_apply_diag .AND. sf(027,3))
      sf_diag%l_ra    = (l_apply_diag .AND. sf(054,3)) ! aerodynamic resistance

      sf_diag%l_et_stom = (l_apply_diag .AND. sf(539,3))
      sf_diag%l_et_stom_surft = (l_apply_diag .AND. sf(540,3))
      sf_diag%l_ftemp = (l_apply_diag .AND. sf(485,3))
      sf_diag%l_fsth  = (l_apply_diag .AND. sf(486,3))
      sf_diag%l_fprf  = (l_apply_diag .AND. sf(489,3))

      !       Set switches for BL diagnostic arrays

      ! counter gradient term for u
      BL_diag%l_rhogamu      = (l_apply_diag .AND. sf(130,3))
      ! counter gradient term for v
      BL_diag%l_rhogamv      = (l_apply_diag .AND. sf(131,3))
      ! counter gradient term for t
      BL_diag%l_rhogamt      = (l_apply_diag .AND. sf(132,3))
      ! counter gradient term for q
      BL_diag%l_rhogamq      = (l_apply_diag .AND. sf(133,3))
      ! mixing length
      BL_diag%l_elm          = (l_apply_diag .AND. sf(134,3))
      ! production rate of TKE by shear
      BL_diag%l_tke_shr_prod = (l_apply_diag .AND. sf(135,3))
      ! production rate of TKE by buoyancy
      BL_diag%l_tke_boy_prod = (l_apply_diag .AND. sf(136,3))
      ! dissipation rate of TKE
      BL_diag%l_tke_dissp    = (l_apply_diag .AND. sf(137,3))
      ! non-dimensional diffusion coef. for u, v
      BL_diag%l_sm           = (l_apply_diag .AND. sf(138,3))
      ! non-dimensional diffusion coef. for t, q
      BL_diag%l_sh           = (l_apply_diag .AND. sf(139,3))
      ! non-gradient buoyancy flux
      BL_diag%l_wb_ng        = (l_apply_diag .AND. sf(140,3))
      ! cloud fraction used in the TKE schemes
      BL_diag%l_cf_trb       = (l_apply_diag .AND. sf(141,3))
      ! condensed water used in the TKE schemes
      BL_diag%l_ql_trb       = (l_apply_diag .AND. sf(142,3))
      ! standard deviation of the distribution function in the TKE schemes
      BL_diag%l_sgm_trb      = (l_apply_diag .AND. sf(143,3))
      ! Heating increment from turbulence dissipation
      BL_diag%l_dtfric       = (l_apply_diag .AND. sf(188,3))
      ! Top of surface mixed layer (Ksurf profile)
      BL_diag%l_smltop       = (l_apply_diag .AND. sf(356,3))
      ! Top of decoupled stratocu layer
      BL_diag%l_dsctop       = (l_apply_diag .AND. sf(357,3))
      ! BL depth diagnosed from Ri>RiCrit
      BL_diag%l_zhlocal      = (l_apply_diag .AND. sf(358,3))
      ! Height of diagnosis parcel top
      BL_diag%l_zhpar        = (l_apply_diag .AND. sf(359,3))
      ! Decoupled stratocu base height
      BL_diag%l_dscbase      = (l_apply_diag .AND. sf(360,3))
      ! BL cloud base height
      BL_diag%l_cldbase      = (l_apply_diag .AND. sf(361,3))
      ! Entrainment rate
      BL_diag%l_weparm       = (l_apply_diag .AND. sf(362,3))
      ! Entrainment rate for decoupled stratocu
      BL_diag%l_weparm_dsc   = (l_apply_diag .AND. sf(363,3))
      ! inversion thickness
      BL_diag%l_dzh          = (l_apply_diag .AND. sf(364,3))
      ! l_apply_diag is actually implicit in suv10m_n, but is
      ! reapplied for consistency with other components of the structure.
      ! x-cpt of 10 m neutral wind
      sf_diag%l_u10m_n       = (l_apply_diag .AND. sf_diag%suv10m_n)
      ! y-cpt of 10 m neutral wind
      sf_diag%l_v10m_n       = (l_apply_diag .AND. sf_diag%suv10m_n)
      ! x-cpt of 10 m pseudostress
      sf_diag%l_mu10m_n      = (l_apply_diag .AND. sf_diag%suv10m_n)
      ! y-cpt of 10 m pseudostress
      sf_diag%l_mv10m_n      = (l_apply_diag .AND. sf_diag%suv10m_n)
      ! Obukhov length, also required for gustiness diagnostic (463)
      BL_diag%l_oblen        = (l_apply_diag .AND. (sf(464,3) .OR. sf(463,3)))
      ! Friction velocity, also required for gustiness diagnostic (463)
      BL_diag%l_ustar        = (l_apply_diag .AND. (sf(465,3) .OR. sf(463,3)))
      ! Surface buoyancy flux
      BL_diag%l_wbsurf       = (l_apply_diag .AND. sf(467,3))
      ! Gradient Richardson number
      BL_diag%l_gradrich     = (l_apply_diag .AND. sf(468,3))
      ! Convective velocity scale
      BL_diag%l_wstar        = (l_apply_diag .AND. sf(466,3))
      ! Stratification
      BL_diag%l_dbdz         = (l_apply_diag .AND. sf(469,3))
      ! Modulus of shear
      BL_diag%l_dvdzm        = (l_apply_diag .AND. sf(470,3))
      ! Momentum diffusivity
      BL_diag%l_rhokm        = (l_apply_diag .AND. sf(471,3))
      ! Thermal diffusivity
      BL_diag%l_rhokh        = (l_apply_diag .AND. sf(472,3))
      ! Turbulent kinetic energy
      BL_diag%l_tke        = sf(473,3) .AND. l_apply_diag               &
                             .OR. i_rhcpt == rhcpt_tke_based            &
                             .OR. l_subgrid_qcl_mp
      ! x component of orographic stress
      BL_diag%l_ostressx     = (l_apply_diag .AND. sf(474,3))
      ! y component of orographic stress
      BL_diag%l_ostressy     = (l_apply_diag .AND. sf(475,3))
      ! local mixing length for momentum
      BL_diag%l_elm3d        = (l_apply_diag .AND. sf(501,3))
      ! local mixing length for scalars
      BL_diag%l_elh3d        = (l_apply_diag .AND. sf(502,3))
      ! local momentum diffusion coefficient
      BL_diag%l_rhokmloc     = (l_apply_diag .AND. sf(503,3))
      ! local scalar diffusion coefficient
      BL_diag%l_rhokhloc     = (l_apply_diag .AND. sf(504,3))
      ! surface driven momentum diffusion coefficient
      BL_diag%l_rhokmsurf    = (l_apply_diag .AND. sf(505,3))
      ! surface driven scalar diffusion coefficient
      BL_diag%l_rhokhsurf    = (l_apply_diag .AND. sf(506,3))
      ! stratocu-top-driven momentum diffusion coefficient
      BL_diag%l_rhokmsc      = (l_apply_diag .AND. sf(507,3))
      ! stratocu-top-driven scalar diffusion coefficient
      BL_diag%l_rhokhsc      = (l_apply_diag .AND. sf(508,3))
      ! stability function for scalars
      BL_diag%l_fh           = (l_apply_diag .AND. sf(511,3))
      ! stability function for momentum
      BL_diag%l_fm           = (l_apply_diag .AND. sf(512,3))
      ! weighting applied to 1D BL scheme in Smag blending
      BL_diag%l_weight1d     = (l_apply_diag .AND. sf(513,3))

    CASE (mt_single_column)
      ! Set the SCM diagnostic flags below to mostly true,
      ! consistent with the settings in NI_imp_ctl
      sf_diag%su10    = .TRUE.
      sf_diag%sv10    = .TRUE.
      sf_diag%sq1p5   = .TRUE.
      sf_diag%st1p5   = .TRUE.
      sf_diag%l_t10m  = .FALSE.
      sf_diag%l_q10m  = .FALSE.
      sf_diag%sfme    = .TRUE.
      sf_diag%sz0heff = .TRUE.

      ! Set diagnostic flags for Single Column Model use
      ! (mostly false as most diags are output separately)
      BL_diag%l_tke      = .TRUE.
      BL_diag%l_dscbase  = .TRUE.
      BL_diag%l_dtfric   = .TRUE.
      BL_diag%l_oblen    = .FALSE.
      BL_diag%l_ustar    = .FALSE.
      BL_diag%l_wbsurf   = .FALSE.
      BL_diag%l_gradrich = .FALSE.
      BL_diag%l_wstar    = .FALSE.
      BL_diag%l_dbdz     = .FALSE.
      BL_diag%l_dvdzm    = .FALSE.
      BL_diag%l_rhokm    = .FALSE.
      BL_diag%l_rhokh    = .FALSE.
      BL_diag%l_ostressx = .FALSE.
      BL_diag%l_ostressy = .FALSE.
      BL_diag%l_cldbase  = .FALSE.
      BL_diag%l_weparm   = .FALSE.
      BL_diag%l_weparm_dsc   = .FALSE.
      BL_diag%l_smltop       = .FALSE.
      BL_diag%l_dsctop       = .FALSE.
      BL_diag%l_zhlocal      = .FALSE.
      BL_diag%l_zhpar        = .FALSE.
      BL_diag%l_dzh          = .FALSE.
      sf_diag%l_u10m_n       = .FALSE.
      sf_diag%l_v10m_n       = .FALSE.
      sf_diag%l_mu10m_n      = .FALSE.
      sf_diag%l_mv10m_n      = .FALSE.
      BL_diag%l_rhogamu      = .FALSE.
      BL_diag%l_rhogamv      = .FALSE.
      BL_diag%l_rhogamt      = .FALSE.
      BL_diag%l_rhogamq      = .FALSE.
      BL_diag%l_elm          = .FALSE.
      BL_diag%l_tke_shr_prod = .TRUE.
      BL_diag%l_tke_boy_prod = .TRUE.
      BL_diag%l_tke_dissp    = .TRUE.
      BL_diag%l_sm           = .FALSE.
      BL_diag%l_sh           = .FALSE.
      BL_diag%l_wb_ng        = .FALSE.
      BL_diag%l_cf_trb       = .FALSE.
      BL_diag%l_ql_trb       = .FALSE.
      BL_diag%l_sgm_trb      = .FALSE.
      BL_diag%l_elm3d        = .FALSE.
      BL_diag%l_elh3d        = .FALSE.
      BL_diag%l_rhokmloc     = .FALSE.
      BL_diag%l_rhokhloc     = .FALSE.
      BL_diag%l_rhokmsurf    = .FALSE.
      BL_diag%l_rhokmsc      = .FALSE.
      BL_diag%l_rhokhsurf    = .FALSE.
      BL_diag%l_rhokhsc      = .FALSE.
      BL_diag%l_weight1d     = .FALSE.
      BL_diag%l_fm           = .FALSE.
      BL_diag%l_fh           = .FALSE.
      sf_diag%l_et_stom   = .FALSE.
      sf_diag%l_et_stom_surft  = .FALSE.   
      sf_diag%l_ftemp        = .FALSE.
      sf_diag%l_fsth         = .FALSE.
      sf_diag%l_fprf         = .FALSE.

    END SELECT ! model_type

    !       Allocate space for those BL diagnostic arrays required and zero
    !       the elements explicitly
    IF (sf_diag%sfme) THEN
      ALLOCATE(sf_diag%fme(                                             &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      sf_diag%fme(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%fme(1,1))
    END IF
    IF (sf_diag%l_ra) THEN
      ALLOCATE(sf_diag%ra(land_points))
      sf_diag%ra(:) = 0.0
    ELSE
      ALLOCATE(sf_diag%ra(1))
    END IF
    IF (sf_diag%l_et_stom .OR. sf_diag%l_et_stom_surft) THEN
      ALLOCATE(sf_diag%et_stom_ij(                                      &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      sf_diag%et_stom_ij(:,:) = 0.0
      ALLOCATE(sf_diag%et_stom_surft(land_points,ntiles))
      sf_diag%et_stom_surft(:,:) = 0.0
      ALLOCATE(sf_diag%resfs_stom(land_points,ntiles))
      sf_diag%resfs_stom(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%et_stom_ij(1,1))
      ALLOCATE(sf_diag%et_stom_surft(1,1))
      ALLOCATE(sf_diag%resfs_stom(1,1))
    END IF

    ! Soil respiration rate modifiers due to temperature, moisture and
    ! vegetation cover: ftemp, fsth and fprf
    IF (sf_diag%l_ftemp) THEN
      ALLOCATE(sf_diag%ftemp(land_points,sm_levels))
      sf_diag%ftemp(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%ftemp(1,1))
    END IF

    IF (sf_diag%l_fsth) THEN
      ALLOCATE(sf_diag%fsth(land_points,sm_levels))
      sf_diag%fsth(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%fsth(1,1))
    END IF

    IF (sf_diag%l_fprf) THEN
      ALLOCATE(sf_diag%fprf(land_points))
      sf_diag%fprf(:) = 0.0
    ELSE
      ALLOCATE(sf_diag%fprf(1))
    END IF


    IF (BL_diag%l_oblen) THEN
      ALLOCATE(BL_diag%oblen(                                           &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%oblen(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%oblen(1,1))
    END IF
    IF (BL_diag%l_ustar) THEN
      ALLOCATE(BL_diag%ustar(                                           &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%ustar(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%ustar(1,1))
    END IF
    IF (BL_diag%l_wbsurf) THEN
      ALLOCATE(BL_diag%wbsurf(                                          &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%wbsurf(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%wbsurf(1,1))
    END IF
    IF (BL_diag%l_gradrich) THEN
      ALLOCATE(BL_diag%gradrich(pdims%i_start:pdims%i_end,              &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%gradrich(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%gradrich(1,1,1))
    END IF
    IF (BL_diag%l_wstar) THEN
      ALLOCATE(BL_diag%wstar(                                           &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%wstar(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%wstar(1,1))
    END IF
    IF (BL_diag%l_dbdz) THEN
      ALLOCATE(BL_diag%dbdz(pdims%i_start:pdims%i_end,                  &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%dbdz(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%dbdz(1,1,1))
    END IF
    IF (BL_diag%l_dvdzm) THEN
      ALLOCATE(BL_diag%dvdzm(pdims%i_start:pdims%i_end,                 &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%dvdzm(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%dvdzm(1,1,1))
    END IF
    IF (BL_diag%l_rhokm) THEN
      ALLOCATE(BL_diag%rhokm(pdims%i_start:pdims%i_end,                 &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokm(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokm(1,1,1))
    END IF
    IF (BL_diag%l_rhokh) THEN
      ALLOCATE(BL_diag%rhokh(pdims%i_start:pdims%i_end,                 &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokh(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokh(1,1,1))
    END IF
    IF (BL_diag%l_tke) THEN
      ALLOCATE(BL_diag%tke(pdims%i_start:pdims%i_end,                   &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%tke(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%tke(1,1,1))
    END IF
    IF (BL_diag%l_ostressx) THEN
      ALLOCATE(BL_diag%ostressx(pdims%i_start:pdims%i_end,              &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%ostressx(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%ostressx(1,1,1))
    END IF
    IF (BL_diag%l_ostressy) THEN
      ALLOCATE(BL_diag%ostressy(pdims%i_start:pdims%i_end,              &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%ostressy(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%ostressy(1,1,1))
    END IF
    IF (BL_diag%l_smltop) THEN
      ALLOCATE(BL_diag%smltop(                                          &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%smltop(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%smltop(1,1))
    END IF
    IF (BL_diag%l_dsctop) THEN
      ALLOCATE(BL_diag%dsctop(                                          &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%dsctop(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%dsctop(1,1))
    END IF
    IF (BL_diag%l_zhlocal) THEN
      ALLOCATE(BL_diag%zhlocal(                                         &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%zhlocal(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%zhlocal(1,1))
    END IF
    IF (BL_diag%l_zhpar) THEN
      ALLOCATE(BL_diag%zhpar(                                           &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%zhpar(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%zhpar(1,1))
    END IF
    IF (BL_diag%l_dscbase) THEN
      ALLOCATE(BL_diag%dscbase(                                         &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%dscbase(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%dscbase(1,1))
    END IF
    IF (BL_diag%l_cldbase) THEN
      ALLOCATE(BL_diag%cldbase(                                         &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%cldbase(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%cldbase(1,1))
    END IF
    IF (BL_diag%l_weparm) THEN
      ALLOCATE(BL_diag%weparm(                                          &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%weparm(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%weparm(1,1))
    END IF
    IF (BL_diag%l_weparm_dsc) THEN
      ALLOCATE(BL_diag%weparm_dsc(                                      &
               pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%weparm_dsc(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%weparm_dsc(1,1))
    END IF
    IF (BL_diag%l_dzh) THEN
      ALLOCATE(BL_diag%dzh(                                             &
             pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
      BL_diag%dzh(:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%dzh(1,1))
    END IF
    IF (sf_diag%l_u10m_n) THEN
      ALLOCATE(sf_diag%u10m_n(udims%i_start:udims%i_end,                &
                              udims%j_start:udims%j_end))
      sf_diag%u10m_n(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%u10m_n(1,1))
    END IF
    IF (sf_diag%l_v10m_n) THEN
      ALLOCATE(sf_diag%v10m_n(vdims%i_start:vdims%i_end,                &
                                vdims%j_start:vdims%j_end))
      sf_diag%v10m_n(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%v10m_n(1,1))
    END IF
    IF (sf_diag%l_mu10m_n) THEN
      ALLOCATE(sf_diag%mu10m_n(udims%i_start:udims%i_end,               &
                                udims%j_start:udims%j_end))
      sf_diag%mu10m_n(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%mu10m_n(1,1))
    END IF
    IF (sf_diag%l_mv10m_n) THEN
      ALLOCATE(sf_diag%mv10m_n(vdims%i_start:vdims%i_end,               &
                                vdims%j_start:vdims%j_end))
      sf_diag%mv10m_n(:,:) = 0.0
    ELSE
      ALLOCATE(sf_diag%mv10m_n(1,1))
    END IF
    ! 10m t and q diagnostics over sea/sea-ice
    IF (sf_diag%l_t10m .OR. sf_diag%l_q10m) THEN
      ALLOCATE( sf_diag%chr10m(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end))
      sf_diag%chr10m = 0.0
    ELSE
      ALLOCATE( sf_diag%chr10m(1,1))
    END IF
    ! Sea and sea ice drag coefficient (momentum)
    IF (sf_diag%l_cd_ssi) THEN
      ALLOCATE( sf_diag%cd_ssi(tdims%i_start:tdims%i_end,               &
                              tdims%j_start:tdims%j_end))
      sf_diag%cd_ssi = 0.0
    ELSE
      ALLOCATE( sf_diag%cd_ssi(1,1))
    END IF
    ! Sea and sea ice drag coefficient (heat)
    IF (sf_diag%l_ch_ssi) THEN
      ALLOCATE( sf_diag%ch_ssi(tdims%i_start:tdims%i_end,               &
                              tdims%j_start:tdims%j_end))
      sf_diag%ch_ssi = 0.0
    ELSE
      ALLOCATE( sf_diag%ch_ssi(1,1))
    END IF
    !
    IF (BL_diag%l_dtfric) THEN
      ALLOCATE(BL_diag%dTfric(pdims%i_start:pdims%i_end,                &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%dTfric(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%dTfric(1,1,1))
    END IF
    IF (BL_diag%l_elm3d) THEN
      ALLOCATE(BL_diag%elm3d(pdims%i_start:pdims%i_end,                 &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%elm3d(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%elm3d(1,1,1))
    END IF
    IF (BL_diag%l_elh3d) THEN
      ALLOCATE(BL_diag%elh3d(pdims%i_start:pdims%i_end,                 &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%elh3d(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%elh3d(1,1,1))
    END IF
    IF (BL_diag%l_rhokmloc) THEN
      ALLOCATE(BL_diag%rhokmloc(pdims%i_start:pdims%i_end,              &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokmloc(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokmloc(1,1,1))
    END IF
    IF (BL_diag%l_rhokhloc) THEN
      ALLOCATE(BL_diag%rhokhloc(pdims%i_start:pdims%i_end,              &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokhloc(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokhloc(1,1,1))
    END IF
    IF (BL_diag%l_rhokmsurf) THEN
      ALLOCATE(BL_diag%rhokmsurf(pdims%i_start:pdims%i_end,             &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokmsurf(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokmsurf(1,1,1))
    END IF
    IF (BL_diag%l_rhokhsurf) THEN
      ALLOCATE(BL_diag%rhokhsurf(pdims%i_start:pdims%i_end,             &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokhsurf(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokhsurf(1,1,1))
    END IF
    IF (BL_diag%l_rhokmsc) THEN
      ALLOCATE(BL_diag%rhokmsc(pdims%i_start:pdims%i_end,               &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokmsc(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokmsc(1,1,1))
    END IF
    IF (BL_diag%l_rhokhsc) THEN
      ALLOCATE(BL_diag%rhokhsc(pdims%i_start:pdims%i_end,               &
                               pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhokhsc(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhokhsc(1,1,1))
    END IF
    IF (BL_diag%l_weight1d) THEN
      ALLOCATE(BL_diag%weight1d(pdims%i_start:pdims%i_end,              &
                                pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%weight1d(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%weight1d(1,1,1))
    END IF
    IF (BL_diag%l_fm) THEN
      ALLOCATE(BL_diag%fm(pdims%i_start:pdims%i_end,                    &
                        pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%fm(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%fm(1,1,1))
    END IF
    IF (BL_diag%l_fh) THEN
      ALLOCATE(BL_diag%fh(pdims%i_start:pdims%i_end,                    &
                        pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%fh(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%fh(1,1,1))
    END IF

    IF (BL_diag%l_rhogamu) THEN
      ALLOCATE(BL_diag%rhogamu(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhogamu(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhogamu(1,1,1))
    END IF

    IF (BL_diag%l_rhogamv) THEN
      ALLOCATE(BL_diag%rhogamv(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhogamv(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhogamv(1,1,1))
    END IF

    IF (BL_diag%l_rhogamt) THEN
      ALLOCATE(BL_diag%rhogamt(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhogamt(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhogamt(1,1,1))
    END IF

    IF (BL_diag%l_rhogamq) THEN
      ALLOCATE(BL_diag%rhogamq(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%rhogamq(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%rhogamq(1,1,1))
    END IF

    IF (BL_diag%l_elm) THEN
      ALLOCATE(BL_diag%elm(pdims%i_start:pdims%i_end,                   &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%elm(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%elm(1,1,1))
    END IF

    IF (BL_diag%l_tke_shr_prod) THEN
      ALLOCATE(BL_diag%tke_shr_prod(pdims%i_start:pdims%i_end,          &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%tke_shr_prod(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%tke_shr_prod(1,1,1))
    END IF

    IF (BL_diag%l_tke_boy_prod) THEN
      ALLOCATE(BL_diag%tke_boy_prod(pdims%i_start:pdims%i_end,          &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%tke_boy_prod(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%tke_boy_prod(1,1,1))
    END IF

    IF (BL_diag%l_tke_dissp) THEN
      ALLOCATE(BL_diag%tke_dissp(pdims%i_start:pdims%i_end,             &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%tke_dissp(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%tke_dissp(1,1,1))
    END IF

    IF (BL_diag%l_sm) THEN
      ALLOCATE(BL_diag%sm(pdims%i_start:pdims%i_end,                    &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%sm(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%sm(1,1,1))
    END IF

    IF (BL_diag%l_sh) THEN
      ALLOCATE(BL_diag%sh(pdims%i_start:pdims%i_end,                    &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%sh(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%sh(1,1,1))
    END IF

    IF (BL_diag%l_wb_ng) THEN
      ALLOCATE(BL_diag%wb_ng(pdims%i_start:pdims%i_end,                 &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%wb_ng(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%wb_ng(1,1,1))
    END IF

    IF (BL_diag%l_cf_trb) THEN
      ALLOCATE(BL_diag%cf_trb(pdims%i_start:pdims%i_end,                &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%cf_trb(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%cf_trb(1,1,1))
    END IF

    IF (BL_diag%l_ql_trb) THEN
      ALLOCATE(BL_diag%ql_trb(pdims%i_start:pdims%i_end,                &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%ql_trb(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%ql_trb(1,1,1))
    END IF

    IF (BL_diag%l_sgm_trb) THEN
      ALLOCATE(BL_diag%sgm_trb(pdims%i_start:pdims%i_end,               &
                            pdims%j_start:pdims%j_end,bl_levels))
      BL_diag%sgm_trb(:,:,:) = 0.0
    ELSE
      ALLOCATE(BL_diag%sgm_trb(1,1,1))
    END IF
  END IF !l_jules_call .OR. .NOT. l_jules_flux
  ! ----------------------------------------------------------------------
  ! Section BL.1 Calculate T at old time level.
  ! Modified to use latest values to avoid time-level inconsistencies
  ! with cloud data.
  ! ---------------------------------------------------------------------
  DO k = 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        t(i,j,k) = theta(i,j,k) * exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  ! ----------------------------------------------------------------------
  ! Calculate the land height over fractional land points:
  ! ----------------------------------------------------------------------
  IF (l_jules_call .OR. .NOT. l_jules_flux) THEN
    l=0
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        z_land(i,j) = 0.0

        IF (land_sea_mask(i,j)) THEN
          l=l+1

          IF ( (l_ctile .AND. fland(l) >  0.0 .AND. fland(l) <  1.0)    &
               .OR. ANY(l_elev_absolute_height)                         &
             ) THEN

            z_land(i,j) = r_theta_levels(i,j,0) - planet_radius
            IF (z_land(i,j) <  0.0)z_land(i,j) = 0.0
          END IF

        END IF
      END DO
    END DO
  END IF !l_jules_call .OR. .NOT. l_jules_flux
  ! ----------------------------------------------------------------------
  ! Section BL.2a Call boundary_layer scheme.
  ! ----------------------------------------------------------------------
  CALL bdy_layr(                                                        &
  ! IN  parameters for iterative SISL scheme
     cycleno,l_jules_call,                                              &
  ! IN values defining field dimensions and subset to be processed :
     ntiles,land_points,bl_levels,                                      &
  ! IN time stepping information
     val_year, val_day_number, val_hour, val_minute, val_second,        &
  ! IN values defining vertical grid of model atmosphere :
     p,p_layer_centres,rho_wet_rsq,rho_wet,rho_dry,rho_wet_tq,          &
     z_uv, z_tq,                                                        &
  ! IN U and V momentum fields.
     u_p, v_p, u_px, v_px, u_0_px, v_0_px,                              &
  ! IN soil/vegetation/land surface data :
     land_index,canopy,catch,catch_snow,hcon,smvccl,                    &
     smvcst,smvcwt,sthf,sthu,sil_orog_land,ho2r2_orog,sd_orog,          &
  ! IN for dust scheme
     soil_layer_moisture,                                               &
  ! IN sea/sea-ice data :
     ice_fract_cat, k_sice,                                             &
  ! IN cloud data :
     cf_bulk, q, qcf, qcl, t,                                           &
  ! IN everything not covered so far :
     photosynth_act_rad, p_star,                                        &
     rad_hr,micro_tends,zh_prev,ddmfx,zhpar,zlcl,                       &
  ! IN Variables for: prescribed surface flux forcing
     flux_e, flux_h, ustar_in, l_spec_z0, z0m_scm, z0h_scm,             &
  ! IN variables required for CLASSIC/mineral dust schemes
     l_aero_classic,soil_clay,soil_sand,                                &
     dust_mrel1,dust_mrel2,dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6, &
  ! IN additional variables for JULES
     snow_tile, z0_tile, z0h_tile_bare, z0m_soil, lw_down,              &
     sw_tile,t_soil,ti,ti_cat,tstar_tile, tsurf_elev_surft,             &
     co2_3d,asteps_since_triffid,                                       &
     cs,frac,canht_ft,lai_ft,fland,flandg,z_land,albsoil,               &
     cos_zenith_angle,                                                  &
  ! SCM Diagnostics (dummy values in full UM) and STASH
     nSCMDpkgs,L_SCMDiags, BL_diag, sf_diag,                            &
  ! INOUT data
     Gs,z0msea,w, etadot,tstar,tstar_sea,tstar_sice_cat,                &
     zh,dzh,ntml,ntpar,l_shallow,cumulus,                               &
  ! INOUT additional variables for JULES
     g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,                    &
  ! INOUT variables for TKE based turbulence schemes
     e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                      &
  ! INOUT variables from bdy_expl1 needed elsewhere
     bq_gb,bt_gb,dtrdz_charney_grid,rdz_charney_grid,                   &
     dtrdz_u,dtrdz_v,rdz_u,rdz_v,k_blend_tq, k_blend_uv,                &
  ! INOUT variables from Jules needed elsewhere
     flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m,cdr10m_n, cd10m_n,&
     fqw, fqw_tile, epot_tile, ftl,ftl_tile, rhokh, rhokm, r_b_dust,    &
     rib_gb, vshr, z0m_eff_gb,z0h_eff_gb, rho_aresist,aresist,resist_b, &
  ! INOUT diagnostics required for soil moisture nudging scheme :
     wt_ext,                                                            &
  ! INOUT variables required in IMP_SOLVER
     alpha1,ashtf,fraca,rhokh_tile,smc,chr1p5m,resfs,z0hssi,z0mssi,     &
     uStarGBM,                                                          &
  ! INOUT additional variables for JULES
     radnet_sea,radnet_sice,rib_tile,rho_aresist_tile,aresist_tile,     &
     resist_b_tile,alpha1_sea,alpha1_sice,ashtf_sea,ashtf_tile,fqw_ice, &
     ftl_ice,resft,rhokh_sice,rhokh_sea,                                &
     z0h_tile,z0m_gb,z0m_tile,chr1p5m_sice,                             &
     g_leaf,gpp_ft,npp_ft, resp_p_ft,resp_s,resp_s_tot,resp_w_ft,       &
     gc,canhc_tile,wt_ext_tile,flake,tile_index,tile_pts,tile_frac,fsmc,&
     vshr_land,vshr_ssi,tstar_land,tstar_ssi,                           &
     dtstar_tile,dtstar_sea,dtstar_sice,hcons,emis_tile,emis_soil,      &
     gpp,npp,resp_p,t1_sd,q1_sd,fb_surf,                                &
  ! OUT variables for message passing
     tau_fd_x, tau_fd_y, rhogamu, rhogamv, f_ngstress,                  &
  ! OUT  diagnostic not requiring stash flags :
     zht, zhnl, shallowc,cu_over_orog,bl_type_1,bl_type_2,bl_type_3,    &
     bl_type_4,bl_type_5,bl_type_6, bl_type_7, bl_w_var,                &
  ! OUT data required for tracer mixing :
     kent, we_lim, t_frac, zrzi,                                        &
     kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                  &
  !OUT variables required for mineral dust scheme
     dust_flux,dust_emiss_frac,u_s_t_tile,u_s_t_dry_tile,u_s_std_tile,  &
  ! OUT data required elsewhere in um system :
     ntdsc, nbdsc, wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt,rib_ssi      &
     )

END IF                    ! on error code = 0

! end of routine NI_bl_ctl
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE ni_bl_ctl
END MODULE ni_bl_ctl_mod
