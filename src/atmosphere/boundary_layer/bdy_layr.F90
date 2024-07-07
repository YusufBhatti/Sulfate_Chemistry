! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Calculate turbulent fluxes of heat, moisture and momentum
!           between (a) surface and atmosphere, (b) atmospheric levels
!           within the boundary layer, and/or the effects of these
!           fluxes on the primary model variables.  The flux of heat
!           into and through the soil is also modelled.  Numerous
!           related diagnostics are also calculated.

!  Programming standard : UMDP 3

!  Documentation: UMDP 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
MODULE bdy_layr_mod

IMPLICIT NONE

SAVE

REAL, ALLOCATABLE ::                                                    &
   recip_l_mo_sea(:,:),                                                 &
                                ! Reciprocal of the surface
                                ! Obukhov length over the sea
                                ! (m-1).
   h_blend_orog(:,:),                                                   &
                                ! Blending height used as part of
                                !  effective roughness scheme
   rhostar(:,:),                                                        &
                                ! Surface air density
   rho_mix(:,:,:),                                                      &
                                ! density on UV (ie. rho) levels;
                                ! used in RHOKH so dry density if
                                ! l_mr_physics is true
   rho_mix_tq(:,:,:),                                                   &
                                ! density on TQ (ie. theta) levels;
                                ! used in non-turb flux integration
                                ! so dry density if l_mr_physics is true
   dzl_charney(:,:,:),                                                  &
                                ! DZL(,K) is depth in m of theta level
                                !  K, i.e. distance from boundary
                                !  K-1/2 to boundary K+1/2.
   rdz(:,:,:),                                                          &
                                ! RDZ(,1) is the reciprocal of the
                                !  height of level 1, i.e. of the
                                !  middle of layer 1.  For K > 1,
                                !  RDZ(,K) is the reciprocal
                                !  of the vertical distance
                                !  from level K-1 to level K.
   qw(:,:,:),                                                           &
   tl(:,:,:),                                                           &
   bt(:,:,:),                                                           &
                                ! A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
   bq(:,:,:),                                                           &
                                ! A buoyancy parameter for clear air
                                ! on p,T,q-levels (full levels).
   bt_cld(:,:,:),                                                       &
                                ! A buoyancy parameter for cloudy air
                                ! on p,T,q-levels (full levels).
   bq_cld(:,:,:),                                                       &
                                ! A buoyancy parameter for cloudy air
                                ! on p,T,q-levels (full levels).
   a_qs(:,:,:),                                                         &
                                ! Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   a_dqsdt(:,:,:),                                                      &
                                ! Saturated lapse rate factor
                                ! on p,T,q-levels (full levels).
   dqsdt(:,:,:)
                                ! Derivative of q_SAT w.r.t. T

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BDY_LAYR_MOD'
CONTAINS

SUBROUTINE bdy_layr (                                                   &
! IN  parameters for iterative SISL scheme
 cycleno,l_jules_call,                                                  &
! IN values defining field dimensions and subset to be processed :
 ntiles,land_pts,bl_levels,                                             &
! IN time stepping information
 val_year, val_day_number, val_hour, val_minute, val_second,            &
! IN values defining vertical grid of model atmosphere :
 p,p_theta_levels, rho_wet_rsq, rho_wet, rho_dry, rho_wet_tq, z_uv,z_tq,&
! IN U, V and W momentum fields.
 u_p, v_p, u_px, v_px, u_0_px, v_0_px,                                  &
! IN soil/vegetation/land surface data :
 land_index,canopy,catch,catch_snow,hcon,smvccl,                        &
 smvcst,smvcwt,sthf,sthu,sil_orog_land,ho2r2_orog,sd_orog,              &
! IN for dust scheme
 soil_layer_moisture,                                                   &
! IN sea/sea-ice data :
 ice_fract_cat, k_sice,                                                 &
! IN cloud data :
 cf_bulk,q,qcf,qcl,t,                                                   &
! IN everything not covered so far :
 photosynth_act_rad,pstar,                                              &
 rad_hr,micro_tends,zh_prev,ddmfx,zhpar,z_lcl,                          &
! IN SCM variables
 flux_e, flux_h, ustar_in, l_spec_z0, z0m_scm, z0h_scm,                 &
! IN variables required for CLASSIC aerosol/mineral dust schemes
 l_aero_classic,soil_clay,soil_sand,                                    &
 dust_mrel1,dust_mrel2,dust_mrel3,dust_mrel4,dust_mrel5,dust_mrel6,     &
! IN additional variables for JULES
 snow_tile,z0_tile,z0h_tile_bare,z0m_soil,lw_down,                      &
 sw_tile,t_soil,ti,ti_cat,tstar_tile,tsurf_elev_surft,                  &
 co2_3d,asteps_since_triffid,                                           &
 cs,frac,canht_ft,lai_ft,fland,flandg,z_land,albsoil,cos_zenith_angle,  &
! SCM Diagnostics (dummy values in full UM) and stash diags
 nSCMDpkgs,L_SCMDiags, BL_diag, sf_diag,                                &
! INOUT data :
 Gs,z0msea,w, etadot,tstar,tstar_sea,tstar_sice_cat,                    &
 zh,dzh,ntml,ntpar,l_shallow,cumulus,                                   &
! INOUT additional variables for JULES
 g_leaf_acc,npp_ft_acc,resp_w_ft_acc,resp_s_acc,                        &
! INOUT variables on TKE based turbulence schemes
 e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                          &
! INOUT variables from bdy_expl1 needed elsewhere
 bq_gb,bt_gb, dtrdz_charney_grid,rdz_charney_grid,                      &
 dtrdz_u,dtrdz_v,rdz_u,rdz_v, k_blend_tq,k_blend_uv,                    &
! INOUT variables from Jules needed elsewhere
 flandfac, fseafac, rhokm_land, rhokm_ssi, cdr10m, cdr10m_n, cd10m_n,   &
 fqw,fqw_tile,epot_tile, ftl,ftl_tile,rhokh,rhokm, r_b_dust,            &
 rib_gb,vshr,z0m_eff_gb,z0h_eff_gb, rho_aresist,aresist,resist_b,       &
! INOUT diagnostics required for soil moisture nudging scheme :
 wt_ext,                                                                &
! INOUT variables required in IMP_SOLVER
 alpha1,ashtf,fraca,rhokh_tile,smc,chr1p5m,resfs,z0hssi,z0mssi,uStarGBM,&
! INOUT additional variables for JULES
 radnet_sea,radnet_sice,rib_tile,rho_aresist_tile,aresist_tile,         &
 resist_b_tile,alpha1_sea,alpha1_sice,ashtf_sea,ashtf_tile,fqw_ice,     &
 ftl_ice,resft,rhokh_sice,rhokh_sea,                                    &
 z0h_tile,z0m_gb,z0m_tile,chr1p5m_sice,                                 &
 g_leaf,gpp_ft,npp_ft, resp_p_ft,resp_s,resp_s_tot,resp_w_ft,           &
 gc,canhc_tile,wt_ext_tile,flake, tile_index,tile_pts,tile_frac,fsmc,   &
 vshr_land,vshr_ssi,tstar_land,tstar_ssi,                               &
 dtstar_tile,dtstar_sea,dtstar_sice,hcons,emis_tile,emis_soil,          &
 gpp,npp,resp_p,t1_sd,q1_sd,fb_surf,                                    &
! OUT variables for message passing
 tau_fd_x, tau_fd_y, rhogamu, rhogamv, f_ngstress,                      &
! OUT Diagnostic not requiring STASH flags :
 zht,zhnl,shallowc,cu_over_orog, bl_type_1,bl_type_2,bl_type_3,         &
 bl_type_4,bl_type_5,bl_type_6, bl_type_7, bl_w_var,                    &
! OUT data required for tracer mixing :
 kent, we_lim, t_frac, zrzi,                                            &
 kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc, zhsc,                      &
!OUT variables required for mineral dust scheme
 dust_flux,dust_emiss_frac, u_s_t_tile,u_s_t_dry_tile,u_s_std_tile,     &
! OUT data required elsewhere in UM system :
 ntdsc,nbdsc,wstar,wthvs,uw0,vw0,taux_p,tauy_p,rhcpt,rib_ssi            &
 )

USE atm_fields_bounds_mod, ONLY:                                        &
                                 udims, vdims, tdims, pdims, pdims_s,   &
                                 wdims, tdims_l
USE atm_step_local, ONLY: dim_cs1, dim_cs2, land_pts_trif, npft_trif,   &
     co2_dim_len,co2_dim_row
USE bl_diags_mod, ONLY: strnewbldiag
USE bl_option_mod, ONLY: i_bl_vn,                                       &
                         i_bl_vn_1a,                                    &
                         i_bl_vn_9b,                                    &
                         i_bl_vn_9c,                                    &
                         on, one_third, l_calc_tau_at_p,                &
                         flux_bc_opt,interactive_fluxes,                &
                         specified_fluxes_only, specified_fluxes_cd
USE cv_run_mod, ONLY: l_jules_flux
USE dust_parameters_mod, ONLY: ndiv, ndivh, ndivl, l_dust, l_dust_diag
USE dynamics_input_mod, ONLY: numcycles
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE jules_soil_mod, ONLY: sm_levels
USE jules_surface_mod, ONLY: i_modiscopt
USE jules_surface_types_mod, ONLY: ntype, npft
USE planet_constants_mod, ONLY: cp, g, grcp, c_virtual, p_zero, kappa,  &
     repsilon, r
USE qsat_mod, ONLY: qsat, qsat_mix
USE sf_diags_mod, ONLY: strnewsfdiag
USE surf_couple_explicit_mod, ONLY: surf_couple_explicit
USE bdy_expl3_mod, ONLY: bdy_expl3
USE level_heights_mod, ONLY: r_rho_levels
USE water_constants_mod, ONLY: lc
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE bdy_expl2_mod, ONLY: bdy_expl2
USE bdy_expl2_1a_mod, ONLY: bdy_expl2_1a
USE dust_calc_emiss_frac_mod, ONLY: dust_calc_emiss_frac
USE dust_srce_mod, ONLY: dust_srce
IMPLICIT NONE

!---------------------------------------------------------------------
!  Inputs :-
!---------------------------------------------------------------------
! (a) Defining horizontal grid and subset thereof to be processed.
INTEGER, INTENT(IN) ::                                                  &
 cycleno
                ! Iteration no

INTEGER, INTENT(IN) ::                                                  &
 ntiles,                                                                &
                                 ! IN No. of land-surface tiles
 land_pts,                                                              &
                                 ! IN No.of land points in whole grid.
 bl_levels,                                                             &
                                 ! IN Max. no. of "boundary" levels
 !time information for current timestep
 val_year, val_day_number, val_hour, val_minute, val_second

! (b) Defining vertical grid of model atmosphere.
REAL, INTENT(IN) ::                                                     &
 p(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,         &
   pdims_s%k_start:bl_levels+1),                                        &
 p_theta_levels(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end, 0:bl_levels+1),              &
 rho_wet_rsq(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,                             &
             pdims_s%k_start:bl_levels+1),                              &
                                       ! IN Density * R**2
 rho_wet(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels+1),                                                  &
                      ! wet density on rho levels (kg/m3)
 rho_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         pdims%k_start:bl_levels+1),                                    &
                      ! dry density on rho levels (kg/m3)
 rho_wet_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                              ! density on TQ (ie. theta) levels;
                              !    used in RHOKM so wet density
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1), &
                                  ! Z_uv(*,K) is height of half
                                  ! level k-1/2.
 z_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),   &
                              ! Z_tq(*,K) is height of full level k.
 u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,               &
     bl_levels),                                                        &
 v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,               &
     bl_levels),                                                        &
 u_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,      &
      bl_levels),                                                       &
 v_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,      &
      bl_levels),                                                       &
 u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),   &
 v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)

! (c) Soil/vegetation/land surface parameters (mostly constant).
INTEGER, INTENT(IN) ::                                                  &
 land_index(land_pts)      ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

REAL, INTENT(IN) ::                                                     &
 canopy(land_pts,ntiles),                                               &
                                 ! IN Surface/canopy water for
                                 !    snow-free land tiles (kg/m2)
 catch(land_pts,ntiles),                                                &
                                 ! IN Surface/canopy water capacity
                                 !    of snow-free land tiles (kg/m2).
 catch_snow(land_pts,ntiles),                                           &
                                 ! IN Snow interception capacity of
                                 !    NLT tile (kg/m2).
 hcon(land_pts),                                                        &
                               ! IN Soil thermal conductivity
!                                     (W/m/K).
   smvccl(land_pts,sm_levels),                                          &
                                 ! IN Critical volumetric SMC (m3/m3
!                                     of soil).
   smvcst(land_pts,sm_levels),                                          &
                                 ! IN Volumetric saturation point
!                                     (m3/m3 of soil).
   smvcwt(land_pts,sm_levels),                                          &
                                 ! IN Volumetric wilting point (m3/m3
!                                     of soil).
   sthf(land_pts,sm_levels),                                            &
                                 ! IN Frozen soil moisture content of
!                                     each layer as a fraction of
!                                     saturation.
   sthu(land_pts,sm_levels),                                            &
                                 ! IN Unfrozen soil moisture content
!                                     of each layer as a fraction of
!                                     saturation.
   sil_orog_land(land_pts),                                             &
                                 ! IN Silhouette area of unresolved
!                                     orography per unit horizontal area
!                                     on land points only.
   ho2r2_orog(land_pts),                                                &
                                 ! IN Standard Deviation of orography.
!                                     equivilent to peak to trough
!                                     height of unresolved orography
!                                     devided by 2SQRT(2) on land
!                                     points only (m)
   sd_orog(land_pts),                                                   &
                                 ! IN Standard Deviation of unresolved
!                                     orography on land points only (m)
   soil_layer_moisture(land_pts,sm_levels),                             &
                                 ! IN soil moisture per layer (kg m-2)
   zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                 ! IN boundary layer height from
!                                     previous timestep
   ddmfx(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! IN Convective downdraught mass-flux
!                                     at cloud base
   zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! IN Height of top of initial
                                 !     parcel ascent
   z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! IN Height of LCL

! (d) Sea/sea-ice data.
REAL, INTENT(IN) ::                                                     &
 ice_fract_cat(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end,nice_use),                     &
                                ! IN Fraction of gridbox covered by
!                                      category sea-ice (decimal fraction).
                       ! If nice_use=1, this is the sum of the categories
   k_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use)
                                   ! IN sea ice surface layer effective
!                                  !    conductivity (W/m2/K)
! (e) Cloud data.
REAL, INTENT(IN) ::                                                     &
 cf_bulk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                                        ! IN Cloud fraction (decimal).
 qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,       &
     tdims_l%k_start:bl_levels),                                        &
                                        ! IN Cloud ice (kg per kg air)
 qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,       &
     tdims_l%k_start:bl_levels),                                        &
                                        ! IN Cloud liquid water
 q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,         &
   tdims_l%k_start:bl_levels),                                          &
                                        ! IN specific humidity
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                        ! IN temperature

! (f) Atmospheric + any other data not covered so far, incl control.
REAL, INTENT(IN) ::                                                     &
 photosynth_act_rad(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end),                         &
                                 ! IN Net downward shortwave radiation
                                 !    in band 1 (w/m2).
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                  ! IN Surface pressure (Pascals).
 rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        2,bl_levels),                                                   &
                                  ! IN (LW,SW) rad heating rate (K/s)
  micro_tends(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              2, bl_levels),                                            &
                         ! Tendencies from microphys within BL levels
                         ! (TL, K/s; QW, kg/kg/s)
  soil_clay (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  soil_sand (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
  dust_mrel1 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel2 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel3 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel4 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel5 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
  dust_mrel6 (pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
 flux_e(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! IN Surf. lat. heat flux   (W/m^2)
 flux_h(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! IN Surf. sens. heat flux  (W/m^2)
 ustar_in(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                 ! IN Surf. friction velocity (m/s)
 tstar_tile(land_pts,ntiles),                                           &
                                 ! IN Surface tile temperatures
 tsurf_elev_surft(land_pts,ntiles),                                     &
                                 ! IN
 z_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! IN    Land height (m).
 z0m_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
 z0h_scm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! IN   Fixed sea-surface roughness
                                 !      lengths for momentum and
                                 !      scalars (m, from SCM namelist)
 t_soil(land_pts,sm_levels),                                            &
                                 ! IN Soil temperatures (K).
 ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
                                 ! IN Sea-ice surface layer
 ti_cat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)

LOGICAL, INTENT(IN) ::                                                  &
 l_spec_z0,                                                             &
                               !    fluxes and roughness lengths
 l_aero_classic,                                                        &
                               ! IN Switch for CLASSIC aerosol scheme
 l_jules_call                  ! switch for whether this call is for jules
                               ! only or the main scheme

! Additional JULES variables
INTEGER, INTENT(IN) ::                                                  &
 asteps_since_triffid
                                 ! IN Number of atmospheric
                                 !    timesteps since last call
                                 !    to TRIFFID.
REAL, INTENT(IN) ::                                                     &
 snow_tile(land_pts,ntiles),                                            &
                                 ! IN Lying snow on tiles (kg/m2)
 z0_tile(land_pts,ntiles),                                              &
                                 ! IN Tile roughness lengths (m).
 z0h_tile_bare(land_pts,ntiles),                                        &
                                 ! IN Tile thermal roughness
                                 ! lengths (m) without snow.
 z0m_soil(land_pts),                                                    &
                                 ! IN bare soil momentum roughness (m)
 lw_down(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! IN Surface downward LW radiation
                                 !    (W/m2).
 sw_tile(land_pts,ntiles),                                              &
                                 ! IN Surface net SW radiation on
                                 !    land tiles (W/m2).
 co2_3d(co2_dim_len,co2_dim_row),                                       &
                                 ! IN 3D CO2 field if required.
 cs(land_pts,dim_cs1),                                                  &
                            ! IN Soil carbon (kg C/m2).
 frac(land_pts,ntype),                                                  &
                                 ! IN Fractions of surface types.
 fland(land_pts),                                                       &
                                 ! IN Land fraction on land tiles.
 flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),   &
                                 ! IN Land fraction on all tiles.
 canht_ft(land_pts,npft),                                               &
                                 ! IN Canopy height (m)
 lai_ft(land_pts,npft),                                                 &
                                 ! IN Leaf area index
 albsoil(land_pts),                                                     &
                                 ! Soil albedo.
 cos_zenith_angle(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! Cosine of the zenith angle

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
 nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
 L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages
!     Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag
!---------------------------------------------------------------------
!  In/outs :-
!---------------------------------------------------------------------
REAL, INTENT(INOUT) ::                                                  &
 Gs(land_pts),                                                          &
                                 ! INOUT "Stomatal" conductance to
                                 !        evaporation (m/s).
 tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                 ! INOUT Surface temperature (K).
 tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                                 ! INOUT  Open sea sfc temperature (K).
 tstar_sice_cat(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,nice_use),                    &
                                 ! INOUT Sea-ice sfc temperature (K).
  z0msea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! INOUT Sea-surface roughness
                                 !       length for momentum (m).
                                 !       NB: same storage is used
                                 !       for Z0V, so the intent is
                                 !       IN for land points.
  w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,0:bl_levels),   &
  etadot(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,           &
         0:bl_levels),                                                  &
  zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! INOUT Height above surface of top of
!                                   boundary layer (metres).
    dzh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! INOUT inversion thickness

! INOUT variables on TKE based turbulence schemes
REAL, INTENT(INOUT) ::                                                  &
 e_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
         bl_levels),                                                    &
                  ! TKE defined on theta levels K-1
 tsq_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                  ! Self covariance of liquid potential temperature
                  ! (thetal'**2) defined on theta levels K-1
 qsq_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                  ! Self covariance of total water
                  ! (qw'**2) defined on theta levels K-1
 cov_trb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                  ! Correlation between thetal and qw
                  ! (thetal'qw') defined on theta levels K-1
 zhpar_shcu(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                  ! Height of mixed layer used to evaluate
                  ! the non-gradient buoyancy flux

LOGICAL, INTENT(INOUT) ::                                               &
  cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                 ! INOUT Logical switch for trade Cu
  l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                         ! INOUT Flag to indicate shallow convection

INTEGER, INTENT(INOUT) ::                                               &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                               ! INOUT Number of model layers in the
                               !    turbulently mixed layer
 ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! INOUT Top level of initial parcel
                               !  ascent. Used in convection scheme.

! Additional JULES variables
REAL, INTENT(INOUT) ::                                                  &
 g_leaf_acc(land_pts,npft),                                             &
                                 ! INOUT Accumulated G_LEAF
 npp_ft_acc(land_pts_trif,npft_trif),                                   &
                                 ! INOUT Accumulated NPP_FT
 resp_w_ft_acc(land_pts_trif,npft_trif),                                &
                                 ! INOUT Accum RESP_W_FT
 resp_s_acc(land_pts_trif,dim_cs1) ! INOUT Accumulated RESP_S

! Variables which are OUT but need to be declared as INOUT because this
! routine is called twice, i.e. they are OUT from its first call, but
! are either needed for the 2nd call or elsewhere in the UM after the
! 2nd call

! ... ones from bdy_expl1
REAL, INTENT(INOUT) ::                                                  &
   dtrdz_charney_grid(pdims%i_start:pdims%i_end,                        &
                      pdims%j_start:pdims%j_end,bl_levels),             &
                                  ! INOUT dt/(rho*r*r*dz) for scalar
                                  !     flux divergence
   rdz_charney_grid(pdims%i_start:pdims%i_end,                          &
                    pdims%j_start:pdims%j_end,bl_levels),               &
                                ! INOUT RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                           ! INOUT dt/(rho*r*r*dz) for
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                           ! INOUT U,V flux divergence
   bq_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! INOUT grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   bt_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! INOUT grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                ! INOUT RDZ (K > 1) on UV-grid.
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels)
                                ! INOUT RDZ (K > 1) on UV-grid.

INTEGER, INTENT(INOUT) ::                                               &
 k_blend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                         ! INOUT Theta level for blending height.
 k_blend_uv(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
                         ! INOUT Rho level for blending height.

! ... ones from surf_couple_explicit
REAL, INTENT(INOUT) ::                                                  &
 uStarGBM(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                               ! GBM surface friction velocity
 fqw(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                                 ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   fqw_tile(land_pts,ntiles),                                           &
                                   ! INOUT Surface FQW for land tiles
   epot_tile(land_pts,ntiles),                                          &
                                   ! INOUT Local EPOT for land tiles.
   ftl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                   ! INOUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   ftl_tile(land_pts,ntiles),                                           &
                                   ! INOUT Surface FTL for land tiles
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                   ! INOUT Exchange coeffs for moisture.
   rhokm(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,   &
         bl_levels),                                                    &
                                   ! INOUT Exchange coefficients for
                                   !     momentum on P-grid
   rib_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! INOUT Mean bulk Richardson number for
!                                     lowest layer.
   vshr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                   ! INOUT Magnitude of surface-to-lowest
!                                     atm level wind shear (m per s).
   vshr_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! INOUT Magnitude of land sfc-to-lowest
!                                     atm level wind shear (m per s).
   vshr_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                   ! INOUT Mag. of mean sea sfc-to-lowest
!                                     atm level wind shear (m per s).
   rho_aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                   ! INOUT RHOSTAR*CD_STD*VSHR
!                                        for CLASSIC aerosol scheme
   aresist(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! INOUT 1/(CD_STD*VSHR)
!                                        for CLASSIC aerosol scheme
   resist_b(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                   ! INOUT (1/CH-1/(CD_STD)/VSHR
!                                         for CLASSIC aerosol scheme
    r_b_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv), &
                                    ! INOUT surface layer resist for dust
    wt_ext(land_pts,sm_levels),                                         &
                                  !INOUT cumulative fraction of transp'n
    tstar_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                   ! INOUT   Land mean sfc temperature (K)
    tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                                   ! INOUT Sea mean sfc temperature (K).
!                                      temperature (K).
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
           pdims_s%j_start:pdims_s%j_end),                              &
  cd10m_n(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end),                               &
         ! variables needed in AP2 for message passing
 gpp(land_pts),                                                         &
                               ! INOUT Gross primary productivity
                               !     (kg C/m2/s).
 npp(land_pts),                                                         &
                               ! INOUT Net primary productivity
                               !     (kg C/m2/s).
 resp_p(land_pts),                                                      &
                               ! INOUT Plant respiration (kg C/m2/s).
 t1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                               ! INOUT Standard deviation of turbulent
!                                   fluctuations of layer 1 temperature;
!                                   for use in initiating convection.
   q1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! INOUT Standard deviation of turbulent
!                                   fluctuations of layer 1 humidity;
!                                   for use in initiating convection.
   fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                ! Surface flux buoyancy over density
                                ! (m^2/s^3)
   z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
   z0h_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                 ! INOUT Effective grid-box roughness
!                                   lengths for momentum and for
!                                   heat, moisture
 alpha1(land_pts,ntiles),                                               &
                                ! INOUT Mean gradient of saturated
                                !     specific humidity with respect
                                !     to temperature between the
                                !     bottom model layer and tile
                                !     surfaces
 ashtf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),   &
                                ! INOUT Coefficient to calculate surface
!                                 heat flux into sea-ice.
   fraca(land_pts,ntiles),                                              &
                                ! INOUT Fraction of surface moisture
                                !     flux with only aerodynamic
                                !     resistance for snow-free land
                                !     tiles.
   rhokh_tile(land_pts,ntiles),                                         &
                                ! INOUT Surface exchange coefficients
                                !     for land tiles
   smc(land_pts),                                                       &
                                ! INOUT Available moisture in the
                                !     soil profile (mm).
   chr1p5m(land_pts,ntiles),                                            &
                                ! INOUT Ratio of coefffs for
                                !     calculation of 1.5m temp for
                                !     land tiles.
   resfs(land_pts,ntiles),                                              &
                                ! INOUT Combined soil, stomatal
                                !     and aerodynamic resistance
                                !     factor for fraction (1-FRACA)
                                !     of snow-free land tiles.
   z0hssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
   z0mssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                ! INOUT Roughness lengths over sea (m).
 radnet_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! INOUT Surface net radiation on
                                 !     open sea (W/m2)
 radnet_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             nice_use),                                                 &
                                 ! INOUT Surface net radiation on
                                 !     sea-ice (W/m2)
 rib_tile(land_pts,ntiles),                                             &
                                 ! INOUT RIB for land tiles.
 rho_aresist_tile(land_pts,ntiles),                                     &
                                 ! INOUT RHOSTAR*CD_STD*VSHR on land
                                 !     tiles for CLASSIC aerosol scheme
 aresist_tile(land_pts,ntiles),                                         &
                                 ! INOUT 1/(CD_STD*VSHR) on land tiles
                                 !     for CLASSIC aerosol scheme
 resist_b_tile(land_pts,ntiles),                                        &
                                 ! INOUT (1/CH-1/CD_STD)/VSHR on land
                                 !     tiles for CLASSIC aerosol scheme
 alpha1_sice(pdims%i_start:pdims%i_end,                                 &
             pdims%j_start:pdims%j_end,nice_use),                       &
                                 ! INOUT ALPHA1 for sea-ice.
 alpha1_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! INOUT ALPHA1 for sea.
 ashtf_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                 ! INOUT Coefficient to calculate surface
                                 !     heat flux into sea.
 ashtf_tile(land_pts,ntiles),                                           &
                                 !INOUT Coefficient to calculate
                                 !     surface heat flux into land
                                 !     tiles.
 fqw_ice(pdims%i_start:pdims%i_end,                                     &
         pdims%j_start:pdims%j_end,nice_use)                 ,          &
                                 ! INOUT Surface FQW for sea-ice
 ftl_ice(pdims%i_start:pdims%i_end,                                     &
         pdims%j_start:pdims%j_end,nice_use),                           &
                                 ! INOUT Surface FTL for sea-ice
 resft(land_pts,ntiles),                                                &
                                 ! INOUT Total resistance factor.
                                 !     FRACA+(1-FRACA)*RESFS for
                                 !     snow-free land, 1 for snow.
 rhokh_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
                                                          nice_use),    &
                                 ! INOUT Surface exchange coefficients
                                 !     for sea-ice
 rhokh_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                 ! INOUT Surface exchange coefficients
                                 !     for sea
 z0h_tile(land_pts,ntiles),                                             &
                                 ! INOUT Tile roughness lengths for heat
                                 !     and moisture (m).
 z0m_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! INOUT Gridbox mean Roughness length
                                 !      for momentum (m).
 z0m_tile(land_pts,ntiles),                                             &
                                 ! INOUT Tile roughness lengths for
                                 !     momentum.
 chr1p5m_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                 ! INOUT CHR1P5M for sea and sea-ice
                                 !     (leads ignored).
 g_leaf(land_pts,npft),                                                 &
                                 ! INOUT Leaf turnover rate (/360days).
 gpp_ft(land_pts,npft),                                                 &
                                 ! INOUT Gross primary productivity
                                 !     on PFTs (kg C/m2/s).
 npp_ft(land_pts,npft),                                                 &
                                 ! INOUT Net primary productivity
                                 !     (kg C/m2/s).
 resp_p_ft(land_pts,npft),                                              &
                                 ! INOUT Plant respiration on PFTs
                                 !     (kg C/m2/s).
 resp_s(land_pts,dim_cs1),                                              &
                             ! INOUT Soil respiration (kg C/m2/s).
 resp_s_tot(dim_cs2),                                                   &
                            ! INOUT Total soil respiration
                            ! (kg C/m2/s).
 resp_w_ft(land_pts,npft),                                              &
                                 ! INOUT Wood maintenance respiration
                                 !     (kg C/m2/s).
 gc(land_pts,ntiles),                                                   &
                                 ! INOUT "Stomatal" conductance to
                                 !      evaporation for land tiles
                                 !      (m/s).
 canhc_tile(land_pts,ntiles),                                           &
                                 ! INOUT Areal heat capacity of canopy
                                 !    for land tiles (J/K/m2).
 wt_ext_tile(land_pts,sm_levels,ntiles),                                &
                                 ! INOUT Fraction of evapotranspiration
                                 !    which is extracted from each
                                 !    soil layer by each tile.
 flake(land_pts,ntiles),                                                &
                                 ! INOUT Lake fraction.
 tile_frac(land_pts,ntiles),                                            &
                                 ! INOUT Tile fractions including
                                 !     snow cover in the ice tile.
 fsmc(land_pts,npft),                                                   &
                                 ! INOUT Moisture availability factor.
 dtstar_tile(land_pts,ntiles),                                          &
                                 ! Change in TSTAR over timestep
                                 ! for land tiles
 dtstar_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! Change is TSTAR over timestep
                                 ! for open sea
 dtstar_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),&
                                 ! Change is TSTAR over timestep
                                 ! for sea-ice
 hcons(land_pts),                                                       &
                                 ! Soil thermal conductivity
                                 ! including water and ice
 emis_tile(land_pts,ntiles),                                            &
                                 ! Emissivity for land tiles
 emis_soil(land_pts)
                                 ! Emissivity of underlying soil

INTEGER, INTENT(INOUT) ::                                               &
 tile_index(land_pts,ntype),                                            &
                                 ! INOUT Index of tile points
 tile_pts(ntype)             ! INOUT Number of tile points

!---------------------------------------------------------------------
!  Outputs :-
!---------------------------------------------------------------------
!  (a) Calculated anyway (use STASH space from higher level) :-
REAL, INTENT(OUT) ::                                                    &
   rib_ssi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                   ! OUT Sea mean bulk Richardson no.
!                                        for lowest layer.
   zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                   ! OUT Max height of turb mixing
   zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                   ! OUT non-local PBL depth
   bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if stable
                                   !     b.l. diagnosed, 0.0 otherwise.
   bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if Sc over
                                   !     stable surface layer diagnosed,
                                   !     0.0 otherwise.
   bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if well
                                   !     mixed b.l. diagnosed,
                                   !     0.0 otherwise.
   bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if
                                   !     decoupled Sc layer (not over
                                   !     cumulus) diagnosed,
                                   !     0.0 otherwise.
   bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if
                                   !     decoupled Sc layer over cumulus
                                   !     diagnosed, 0.0 otherwise.
   bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if a
                                   !     cumulus capped b.l. diagnosed,
                                   !     0.0 otherwise.
   bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                 ! OUT Indicator set to 1.0 if a
                                   !     Shear-dominated unstable b.l.
                                   !     diagnosed, 0.0 otherwise.
    we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    ! OUT rho*entrainment rate implied b
                                    !     placing of subsidence
    zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                    ! OUT (z-z_base)/(z_i-z_base)
    t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                    ! OUT a fraction of the timestep
    we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                    ! OUT rho*entrainment rate implied b
                                    !     placing of subsidence
    zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                    ! OUT (z-z_base)/(z_i-z_base)
    t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),  &
                                    ! OUT a fraction of the timestep
    zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! OUT Top of decoupled layer
    dust_flux(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,ndiv),&
                                    ! OUT dust production flux(kg m-2 s-1)
    dust_emiss_frac(land_pts,ntiles),                                   &
                                    ! OUT dust emiss frac on each tile
    u_s_t_tile(land_pts,ntiles,ndivh),                                  &
                          ! OUT threshold friction vel on tiles for dust
    u_s_t_dry_tile(land_pts,ntiles,ndivh),                              &
                          ! OUT dry threshold friction velocity on tiles
    u_s_std_tile(land_pts,ntiles),                                      &
                          ! OUT Surface friction velocity
    wstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT Convective velocity scale (m/s)
    wthvs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                   ! OUT surface flux of thv (Km/s)
    shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                   ! OUT Shallow Cu diagnostic
                                   !   Indicator set to 1.0 if shallow,
                                   !   0.0 if not shallow or not cumulus
    cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                   ! OUT Indicator for cumulus
                                   !     over steep orography
                                   !   Indicator set to 1.0 if true,
                                   !   0.0 if false. Exclusive.

REAL, INTENT(OUT) :: bl_w_var ( tdims%i_start : tdims%i_end,            &
                                tdims%j_start : tdims%j_end,            &
                                            2 : tdims%k_end+1 )

REAL, INTENT(OUT) ::                                                    &
  tau_fd_x(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end,bl_levels),                    &
  tau_fd_y(pdims_s%i_start:pdims_s%i_end,                               &
           pdims_s%j_start:pdims_s%j_end,bl_levels),                    &
  rhogamu(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
  rhogamv(pdims_s%i_start:pdims_s%i_end,                                &
          pdims_s%j_start:pdims_s%j_end ,bl_levels),                    &
  f_ngstress(pdims_s%i_start:pdims_s%i_end,                             &
             pdims_s%j_start:pdims_s%j_end,2:bl_levels)

INTEGER, INTENT(OUT) ::                                                 &
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                 ! OUT Top level for turb mixing in
!                                           any decoupled Sc layer
   nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT Bottom level of any decoupled
                                   !     turbulently-mixed Sc layer.
    kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                   ! OUT grid-level of SML inversion
    kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                   ! OUT grid-level of DSC inversion

!  (b) Not passed between lower-level routines (not in workspace at this
!      level) :-

!-2 Genuinely output, needed by other atmospheric routines :-
REAL, INTENT(OUT) ::                                                    &
  uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                      !OUT U-component of surface wind stress (P-grid)
  vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                      ! V-component of surface wind stress (P-grid)
REAL, INTENT(OUT) :: taux_p( pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,                 &
                             bl_levels )
REAL, INTENT(OUT) :: tauy_p( pdims%i_start:pdims%i_end,                 &
                             pdims%j_start:pdims%j_end,                 &
                             bl_levels )
         ! Wind stresses on theta-levels, on p-grid (for convection)
REAL, INTENT(OUT) :: rhcpt(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,1:tdims%k_end)

!-----------------------------------------------------------------------
!  Workspace :-
INTEGER ::                                                              &
idiv,                                                                   &
          ! loop counter, mineral dust divisions
m !loop counter

REAL ::                                                                 &
 zblend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! Height of theta level blending height
 zblend_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! Height of u,v level blending height
 qw_blend(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
 tl_blend(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
 u_blend_px(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end),                             &
 v_blend_px(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end),                             &
 bt_blend(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
 bq_blend(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
 zh_jules(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                              ! BL depth to pass to Jules
   dust_flux_tile(land_pts,ntiles,ndiv),                                &
                                            !production flux from tiles
   cd_std_dust(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                   !  Bulk transfer coef. for
                              ! momentum, excluding orographic effects
   clay_land(land_pts),                                                 &
                           ! soil clay fraction on land pts
   sand_land(land_pts),                                                 &
                           ! soil sand fraction on land pts
   pstar_land(land_pts),                                                &
                            ! surface pressure on land pts
   rhostar_land(land_pts),                                              &
                              ! surface air density on land pts
   mrel_land(land_pts,ndivl),                                           &
                                ! soil size fraction on land pts
   rib_land(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                   !     Land mean bulk Richardson no.
!                                        for lowest layer.

REAL ::                                                                 &
rholem,                                                                 &
             ! surface density in LEM
tv1_sd,                                                                 &
             ! virt T standard deviation (approx)
w_m,     &   ! velocity scale
dqsdt_star,& ! derivative of qsat wrt tstar
wthvbar, &   ! flux of theta_v
ch,      &   ! scalar transfer coefficient
theta1       ! level 1 potential temperature

! Variables for computing explicit stresses on p-grid
REAL, ALLOCATABLE :: taux_land(:,:)
REAL, ALLOCATABLE :: tauy_land(:,:)
REAL, ALLOCATABLE :: taux_ssi(:,:)
REAL, ALLOCATABLE :: tauy_ssi(:,:)
LOGICAL :: l_calc_at_p

REAL, ALLOCATABLE :: qs_star(:,:)
             ! qsat of tstar
REAL, ALLOCATABLE :: z1_uv_top(:,:)
             ! Height of top of lowest uv layer above the surface
REAL, ALLOCATABLE :: z1_tq_top(:,:)
             ! Height of top of lowest Tq layer above the surface
INTEGER ::                                                              &
 i,j,k,l,n

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='BDY_LAYR'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_jules_call .OR. .NOT. l_jules_flux) THEN

  IF (flux_bc_opt > interactive_fluxes) THEN
    ! For specified surface fluxes impose uniform TSTAR,
    ! specified or calculated in CONV_DIAG to be consistent with fluxes
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        tstar_land(i,j) = tstar(i,j)
        tstar_sea(i,j)  = tstar(i,j)
        tstar_sice_cat(i,j,:) = tstar(i,j)
        tstar_ssi(i,j)  = tstar(i,j)
      END DO
    END DO
  END IF
  !-----------------------------------------------------------------------
  !     Allocate arrays for conservative discretization of the
  !     surface layer.
  IF (i_modiscopt == on) THEN
    ALLOCATE(z1_uv_top(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))
    ALLOCATE(z1_tq_top(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end))
  ELSE
    ALLOCATE(z1_uv_top(1,1))
    ALLOCATE(z1_tq_top(1,1))
  END IF

  ! OUT arrays from bdy_expl1 which need passing to bdy_expl2
  ALLOCATE(rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, &
                   bl_levels+1))
  ALLOCATE(rho_mix_tq(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,bl_levels))
  ALLOCATE(dzl_charney(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end,bl_levels))
  ALLOCATE(rdz(pdims_s%i_start:pdims_s%i_end,                           &
               pdims_s%j_start:pdims_s%j_end,bl_levels))
  ALLOCATE(qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels))
  ALLOCATE(tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels))
  ALLOCATE(bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels))
  ALLOCATE(bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels))
  ALLOCATE(bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  bl_levels))
  ALLOCATE(bq_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  bl_levels))
  ALLOCATE(a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                bl_levels))
  ALLOCATE(a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, &
                   bl_levels))
  ALLOCATE(dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 bl_levels))

  ! DEPENDS ON: bdy_expl1
  CALL bdy_expl1 (                                                      &
  ! IN values defining vertical grid of model atmosphere :
    bl_levels, p_theta_levels, rho_wet_rsq, rho_wet, rho_dry,           &
  ! IN U, V and W momentum fields.
    u_px, v_px,                                                         &
  ! IN cloud data :
    cf_bulk,q,qcf,qcl,t,                                                &
  ! OUT
    dtrdz_charney_grid,rdz_charney_grid,dtrdz_u,dtrdz_v,                &
    rdz_u,rdz_v,rho_mix,rho_mix_tq,dzl_charney,rdz,                     &
    zblend_tq,z1_tq_top,zblend_uv,z1_uv_top,                            &
    k_blend_tq,k_blend_uv,qw,tl,qw_blend,tl_blend,u_blend_px,v_blend_px,&
    bt,bq,bt_blend,bq_blend,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt&
    )

  ! OUT arrays from surf_couple_explicit which need passing to bdy_expl2
  ALLOCATE(recip_l_mo_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  ALLOCATE(h_blend_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))
  ALLOCATE(rhostar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))

  IF (l_jules_flux) THEN
    ! zh from the previous timestep is all we have available if calling Jules
    ! before convection diagnosis
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        zh_jules(i,j) = zh_prev(i,j)
      END DO
    END DO
  ELSE
    ! use the value from the convection diagnosis
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        zh_jules(i,j) = zh(i,j)
      END DO
    END DO
  END IF

  CALL surf_couple_explicit(                                            &
    !Arguments used by JULES-standalone
    !Misc INTENT(IN) DONE
    bq_blend, bt_blend, zh_jules, photosynth_act_rad,                   &
    val_year, val_day_number, val_hour, val_minute, val_second,         &
    !Forcing INTENT(IN) DONE
    qw_blend, tl_blend, pstar, lw_down,                                 &
    !Fluxes INTENT(IN) DONE
    sw_tile, tstar,                                                     &
    !INOUT Diagnostics, in sf_diags_mod
    sf_diag,                                                            &
    !Fluxes INTENT(OUT) DONE
    fqw,ftl,ftl_tile, fqw_tile, fqw_ice, ftl_ice, fsmc, emis_tile,      &
    !Misc INTENT(OUT)
    !rhokms needed for message passing
    radnet_sea, radnet_sice, rhokm, rhokm_land, rhokm_ssi,              &
    !Out of explicit and into implicit only INTENT(OUT)
    !cdr10m needed for message passing
    cdr10m, cdr10m_n, cd10m_n,                                          &
    alpha1, alpha1_sea, alpha1_sice, ashtf, ashtf_sea, ashtf_tile,      &
    epot_tile,                                                          &
    !rhokh needed in UM
    fraca, resfs, resft, rhokh, rhokh_tile, rhokh_sice, rhokh_sea,      &
    dtstar_tile, dtstar_sea, dtstar_sice,                               &
    z0hssi, z0h_tile, z0mssi, z0m_tile, chr1p5m, chr1p5m_sice,          &
    canhc_tile, wt_ext_tile, flake,                                     &
    !Out of explicit and into extra only INTENT(OUT)
    hcons,                                                              &
    !Out of explicit and into implicit and extra INTENT(OUT)
    tile_frac,                                                          &
    !Additional arguments for the UM-----------------------------------------
    !JULES prognostics module
    !IN
    canopy, snow_tile, k_sice, cs, canht_ft, lai_ft,                    &
    t_soil, tsurf_elev_surft, ti, ti_cat,                               &
    tstar_tile,                                                         &
    !INOUT, OUT, OUT,INOUT
    z0msea, smc, gc, gs,                                                &
    !JULES ancil_info module
    !IN
    land_pts, zblend_uv, zblend_tq, land_index, ntiles, ice_fract_cat,  &
    frac,                                                               &
    !OUT
    tile_index, tile_pts,                                               &
    !JULES coastal module
    !4 IN, vshr_ both OUT
    fland, flandg, tstar_sea, tstar_sice_cat, vshr_land, vshr_ssi,      &
    !JULES aero module
    !co2_3d IN, rest OUT for tracer mixing
    co2_3d, rho_aresist, aresist, resist_b, rho_aresist_tile,           &
    aresist_tile, resist_b_tile, r_b_dust, cd_std_dust, u_s_std_tile,   &
    !JULES trifctl module
    !IN
    asteps_since_triffid,                                               &
    !INOUT
    g_leaf_acc, npp_ft_acc, resp_w_ft_acc, resp_s_acc,                  &
    !OUT
    gpp, npp, resp_p, g_leaf, gpp_ft, npp_ft, resp_p_ft, resp_s,        &
    resp_w_ft,                                                          &
    !JULES p_s_parms module
    !IN
    catch, catch_snow, hcon, smvccl, smvcst, smvcwt, sthf, sthu,z0_tile,&
    !IN
    z0h_tile_bare, z0m_soil, albsoil, cos_zenith_angle, soil_clay,      &
    !JULES orog module
    !2 IN, 2 OUT
    ho2r2_orog, sil_orog_land, h_blend_orog, z0m_eff_gb,                &
    !JULES u_v_grid module
    !IN for message passing
    u_blend_px, v_blend_px, u_0_px, v_0_px,                             &
    !JULES switches module  **squish**
    !IN SCM related
    l_spec_z0,                                                          &
    !JULES c_elevate module
    !IN
    z_land,                                                             &
    !Not in a JULES module
    !IN
    numcycles, cycleno, z1_uv_top, z1_tq_top, ddmfx,                    &
    !3 IN, 1 OUT requiring STASH flag
    l_aero_classic, z0m_scm, z0h_scm,                                   &
    !OUT not requiring STASH flag
    recip_l_mo_sea, rib_gb, rib_tile,                                   &
    !OUT 2 message passing, 1 soil moisture nudging, rest of UM
    flandfac, fseafac, wt_ext, fb_surf, ustargbm, t1_sd, q1_sd, rhostar,&
    !OUT
    z0m_gb, z0h_eff_gb, vshr, resp_s_tot, emis_soil                     &
    )

  DEALLOCATE(z1_uv_top)
  DEALLOCATE(z1_tq_top)

  IF (flux_bc_opt > interactive_fluxes) THEN
    ! Surface fluxes calculated in SFEXPL are substituted with
    ! forcing values.  NOTE: Surface calculation also made
    ! explicit (time weight set to zero).
    IF (flux_bc_opt == specified_fluxes_cd) THEN
      ! Set rhokm using input ustar_in
      DO i = pdims%i_start, pdims%i_end
        DO j = pdims%j_start, pdims%j_end
          uStarGBM(i,j) = ustar_in(i,j)
          rhokm(i,j,1) = rhostar(i,j)*uStarGBM(i,j)*uStarGBM(i,j)/vshr(i,j)
          IF (flandg(i,j) > 0.0) THEN
            rhokm_land(i,j) = rhostar(i,j)*uStarGBM(i,j)*uStarGBM(i,j)/ &
                              vshr_land(i,j)
          END IF
          IF (flandg(i,j) < 1.0) THEN
            rhokm_ssi(i,j)  = rhostar(i,j)*uStarGBM(i,j)*uStarGBM(i,j)/ &
                              vshr_ssi(i,j)
          END IF
        END DO
      END DO
    END IF

    DO i = pdims%i_start, pdims%i_end
      DO j = pdims%j_start, pdims%j_end
        !..Converts Fluxes from W/m^2 to rho*K/s
        rholem = rhostar(i,j)

        !..If comparing against LES with rho ne rhostar then match
        ! w'theta' rather than rho*wtheta
        ! RHOLEM = 1.0

        fqw(i,j,1)   = (rhostar(i,j)*flux_e(i,j))/(lc*rholem)
        ftl(i,j,1)   = (rhostar(i,j)*flux_h(i,j))/(cp*rholem)

        fb_surf(i,j) = g * ( bt_gb(i,j,1)*ftl(i,j,1) +                  &
                             bq_gb(i,j,1)*fqw(i,j,1) ) /rhostar(i,j)
        IF ( fb_surf(i,j)  >   0.0) THEN
          w_m        = ( 0.25*zh(i,j)*fb_surf(i,j) +                    &
                ustargbm(i,j)*ustargbm(i,j)*ustargbm(i,j) ) ** one_third
          t1_sd(i,j) = 1.93 * ftl(i,j,1) / (rhostar(i,j) * w_m)
          q1_sd(i,j) = 1.93 * fqw(i,j,1) / (rhostar(i,j) * w_m)
          tv1_sd     = t(i,j,1) * ( bt_gb(i,j,1)*t1_sd(i,j) +           &
                                    bq_gb(i,j,1)*q1_sd(i,j) )
          t1_sd(i,j) = MAX ( 0.0 , t1_sd(i,j) )
          q1_sd(i,j) = MAX ( 0.0 , q1_sd(i,j) )
          IF (tv1_sd  <=  0.0) THEN
            t1_sd(i,j) = 0.0
            q1_sd(i,j) = 0.0
          END IF
        ELSE
          t1_sd(i,j) = 0.0
          q1_sd(i,j) = 0.0
        END IF
      END DO    ! J
    END DO    ! I
    DO i = 1, land_pts
      DO l = 1, ntiles
        fqw_tile(i,l) = fqw(1,1,1)
        ftl_tile(i,l) = ftl(1,1,1)
      END DO ! L
    END DO ! I

    IF ( l_jules_flux .AND. flux_bc_opt == specified_fluxes_only ) THEN
      ! recalculate the surface temperature here to be consistent with
      ! fixed fluxes, since this isn't done in conv_surf_flux or specified
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          ! start with simple extrapolation from level 1
          tstar(i,j) = t(i,j,1) + grcp * z_tq(i,j,1)
        END DO
      END DO

      ALLOCATE(qs_star(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end))

      IF (l_mr_physics) THEN
        CALL qsat_mix(qs_star,tstar,pstar,pdims%i_len,pdims%j_len)
      ELSE
        CALL qsat(qs_star,tstar,pstar,pdims%i_len,pdims%j_len)
      END IF

      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end

          dqsdt_star = repsilon * lc * qs_star(i,j) /                   &
                       ( r * tstar(i,j) * tstar(i,j) )

          theta1 = t(i,j,1) * (p_zero/p_theta_levels(i,j,1))**kappa

          wthvbar = theta1 *                                            &
                    (1.0+c_virtual*q(i,j,1)-qcl(i,j,1)-qcf(i,j,1)) *    &
                    fb_surf(i,j) / g

          ch = rhokh(i,j,1) / ( vshr(i,j) * rhostar(i,j) )

          ! Now more complicated formula based on fluxes
          tstar(i,j) = ( theta1 + wthvbar/(ch*MAX(0.1,vshr(i,j))) -     &
                       c_virtual * theta1 *                             &
                       (qs_star(i,j)-q(i,j,1)-dqsdt_star*tstar(i,j)) )  &
                       / ( (p_zero/pstar(i,j))**kappa +                 &
                           c_virtual * theta1 * dqsdt_star )

          tstar_land(i,j) = tstar(i,j)
          tstar_sea(i,j)  = tstar(i,j)
          tstar_sice_cat(i,j,:) = tstar(i,j)
          tstar_ssi(i,j)  = tstar(i,j)

        END DO
      END DO

      DEALLOCATE(qs_star)

    END IF ! l_jules_flux

  END IF ! flux_bc_opt

END IF !l_jules_call .OR. .NOT. l_jules_flux

IF (.NOT. l_jules_call  .OR. .NOT. l_jules_flux) THEN

  IF ( i_bl_vn == i_bl_vn_9b .OR. i_bl_vn == i_bl_vn_9c  ) THEN

    CALL bdy_expl2 (                                                    &
    ! IN values defining vertical grid of model atmosphere :
      bl_levels,p_theta_levels,land_pts,land_index,                     &
    ! IN U, V and W momentum fields.
      u_p,v_p, u_0_px, v_0_px,                                          &
    ! IN from other part of explicit boundary layer code
      rho_mix,rho_wet_tq,rho_mix_tq,dzl_charney,rdz,rdz_charney_grid,   &
      z_tq,z_uv,rhostar,bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,   &
      dqsdt,recip_l_mo_sea, flandg,                                     &
      rib_gb, sil_orog_land,z0m_eff_gb,                                 &
    ! IN cloud/moisture data :
      cf_bulk,q,qcf,qcl,t,qw,tl,                                        &
    ! IN everything not covered so far :
      rad_hr,micro_tends,fb_surf,ustargbm,pstar,tstar,h_blend_orog,     &
      zh_prev, zhpar,z_lcl,ho2r2_orog,sd_orog,                          &
    ! SCM Diagnostics (dummy values in full UM) & stash diag
      nSCMDpkgs,L_SCMDiags,BL_diag,                                     &
    ! INOUT variables
      zh,dzh,ntml,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,w,etadot, &
      t1_sd,q1_sd,                                                      &
    ! OUT New variables for message passing
      tau_fd_x, tau_fd_y, f_ngstress,                                   &
    ! OUT Diagnostic not requiring STASH flags :
      zht,zhnl,shallowc,cu_over_orog,                                   &
      bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,      &
      bl_type_7,                                                        &
    ! OUT data for turbulent generation of mixed-phase cloud:
      bl_w_var,                                                         &
    ! OUT data required for tracer mixing :
      kent, we_lim, t_frac, zrzi,                                       &
      kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,                       &
    ! OUT data required elsewhere in UM system :
      zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0,rhcpt                        &
      )

  ELSE IF ( i_bl_vn == i_bl_vn_1a ) THEN

    CALL bdy_expl2_1a (                                                 &
    ! IN values defining vertical grid of model atmosphere :
      bl_levels,p_theta_levels,land_pts,land_index, cycleno,            &
    ! IN U, V and W momentum fields.
      u_p,v_p, u_0_px, v_0_px,                                          &
    ! IN variables for TKE scheme
      pstar,p,                                                          &
    ! IN from other part of explicit boundary layer code
      rho_mix,rho_wet_tq,rdz,rdz_charney_grid,                          &
      z_tq,z_uv,bt,bt_gb,bq_gb, flandg,                                 &
      rib_gb, sil_orog_land,z0m_eff_gb,                                 &
    ! IN cloud/moisture data :
      q,qcf,qcl,t,qw,tl,                                                &
    ! IN everything not covered so far :
      fb_surf,ustargbm,h_blend_orog,                                    &
      zh_prev, ho2r2_orog,sd_orog,                                      &
    ! SCM Diagnostics (dummy values in full UM) & stash diag
      nSCMDpkgs,L_SCMDiags,BL_diag,                                     &
    ! INOUT variables
      zh,ntml,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,              &
    ! INOUT variables on TKE based turbulence schemes
      e_trb, tsq_trb, qsq_trb, cov_trb, zhpar_shcu,                     &
    ! OUT New variables for message passing
      tau_fd_x, tau_fd_y, rhogamu, rhogamv,                             &
    ! OUT Diagnostic not requiring STASH flags :
      zht,shallowc,cu_over_orog,                                        &
      bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,      &
      bl_type_7,                                                        &
    ! OUT data required for tracer mixing :
      kent, we_lim, t_frac, zrzi,                                       &
      kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,                       &
    ! OUT data required elsewhere in UM system :
      zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0                              &
      )

  END IF

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rib_land(i,j)=0.0
      rib_ssi(i,j)=0.0
    END DO
  END DO

  DO n = 1, ntiles
    DO k = 1, tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end
      rib_land(i,j)=rib_land(i,j) +                                     &
        rib_tile(l,n)*tile_frac(l,n)
    END DO
  END DO

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (flandg(i,j) <  1.0)                                           &
        rib_ssi(i,j)=(rib_gb(i,j)-rib_land(i,j)*flandg(i,j))            &
          /(1.0-flandg(i,j))
    END DO
  END DO

  !-----------------------------------------------------------------------
  ! Calculate explicit momentum fluxes on p-grid, on all levels, if requested
  !-----------------------------------------------------------------------
  IF ( l_calc_tau_at_p ) THEN

    ! Allocate temporary space for outputs from bdy_expl3 that we don't need
    ALLOCATE( taux_land( pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end ) )
    ALLOCATE( tauy_land( pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end ) )
    ALLOCATE( taux_ssi( pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end ) )
    ALLOCATE( tauy_ssi( pdims%i_start:pdims%i_end,                      &
                        pdims%j_start:pdims%j_end ) )

    l_calc_at_p = .TRUE.
    CALL bdy_expl3 (                                                    &
    ! IN grid related variables
      bl_levels, l_calc_at_p,                                           &
      pdims_s, pdims_s, pdims, pdims, pdims,                            &
      pdims, pdims,                                                     &
    ! Note: here the inputs all have haloes (pdims_s) while the outputs do not.
    ! IN SCM diags
      nSCMDpkgs,L_SCMDiags,                                             &
    ! IN variables used in flux calculations
      u_p, v_p, u_0_px, v_0_px, rhokm_land,rhokm_land,flandfac,flandfac,&
      rhokm_ssi, rhokm_ssi, fseafac, fseafac, flandg, flandg,           &
      zhnl, rdz(:,:,2:bl_levels), rdz(:,:,2:bl_levels),                 &
      rhokm, rhokm, tau_fd_x, tau_fd_y,                                 &
      rhogamu(:,:,2:bl_levels), rhogamv(:,:,2:bl_levels),               &
      f_ngstress, f_ngstress,                                           &
    ! OUT explicit momentum fluxes
      taux_land, tauy_land, taux_ssi, tauy_ssi, taux_p, tauy_p          &
      )

    DEALLOCATE( tauy_ssi )
    DEALLOCATE( taux_ssi )
    DEALLOCATE( tauy_land )
    DEALLOCATE( taux_land )

  END IF  ! ( l_calc_tau_at_p )


  !-----------------------------------------------------------------------
  ! Mineral dust production
  !-----------------------------------------------------------------------
  IF (l_dust .OR. l_dust_diag) THEN
    !initialisation
    dust_flux(:,:,:)=0.0
    dust_flux_tile(:,:,:) = 0.0
    u_s_t_tile(:,:,:) = 0.0
    u_s_t_dry_tile(:,:,:) = 0.0

    !put fields into land arrays
    DO l = 1, land_pts
      j = (land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end
      pstar_land(l) = pstar(i,j)
      rhostar_land(l) = rhostar(i,j)
      sand_land(l) = soil_sand(i,j)
      clay_land(l) = soil_clay(i,j)
    END DO !LAND_PTS

    DO l = 1, land_pts
      j = (land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end
      mrel_land(l,1) = dust_mrel1(i,j)
      mrel_land(l,2) = dust_mrel2(i,j)
      mrel_land(l,3) = dust_mrel3(i,j)
      mrel_land(l,4) = dust_mrel4(i,j)
      mrel_land(l,5) = dust_mrel5(i,j)
      mrel_land(l,6) = dust_mrel6(i,j)
    END DO !LAND_PTS

    CALL dust_srce(                                                     &
    ! IN arguments
             land_pts,ntiles,tile_pts,tile_index,fland,                 &
             tstar_tile,rhostar_land,soil_layer_moisture,snow_tile,     &
             u_s_std_tile,mrel_land,clay_land,sand_land,ho2r2_orog,     &
    ! OUT arguments
             dust_flux_tile,u_s_t_tile,u_s_t_dry_tile                   &
             )

    ! Get the fraction within each tile which is bare soil, for the purpose
    ! of dust emission:
    CALL dust_calc_emiss_frac(                                          &
  land_pts,ntiles,tile_pts,tile_index,frac,lai_ft,dust_emiss_frac       &
    )

    ! Produce a total dust flux over all tiles, by looping through tiles and
    ! multiplying the flux on that tile by the dust_emiss_frac, and summing.
    DO idiv = 1, ndiv
      DO m = 1, ntiles
        DO n = 1, tile_pts(m)
          l = tile_index(n,m)
          j = (land_index(l)-1)/pdims%i_end + 1
          i = land_index(l) - (j-1)*pdims%i_end
          dust_flux(i,j,idiv) = dust_flux(i,j,idiv) +                   &
           dust_flux_tile(l,m,idiv)*dust_emiss_frac(l,m)
        END DO !TILE_PTS
      END DO !NTILES
    END DO !NDIV

  END IF !L_DUST .OR. L_DUST_DIAG

  ! deallocate variables passed from Jules to bdy_expl2
  DEALLOCATE(rhostar)
  DEALLOCATE(h_blend_orog)
  DEALLOCATE(recip_l_mo_sea)
  ! deallocate variables passed form bdy_expl1 to bdy_expl2
  DEALLOCATE(dqsdt)
  DEALLOCATE(a_dqsdt)
  DEALLOCATE(a_qs)
  DEALLOCATE(bq_cld)
  DEALLOCATE(bt_cld)
  DEALLOCATE(bq)
  DEALLOCATE(bt)
  DEALLOCATE(tl)
  DEALLOCATE(qw)
  DEALLOCATE(rdz)
  DEALLOCATE(dzl_charney)
  DEALLOCATE(rho_mix_tq)
  DEALLOCATE(rho_mix)

END IF !.NOT. l_jules_call .OR. .NOT. l_jules_flux

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bdy_layr
END MODULE bdy_layr_mod
