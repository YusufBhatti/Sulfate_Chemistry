! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE IMP_SOLVER---------------------------------------------

!  Purpose: implicit solver for diffusion equation
!           split from bdy_layr routine

!  Programming standard : UMDP 3

!  Documentation: UMDP 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
!    Arguments :-
SUBROUTINE imp_solver (                                                 &
! IN values defining field dimensions and subset to be processed :
 ntiles, land_pts,                                                      &
! IN values defining vertical grid of model atmosphere :
 bl_levels, alpha_cd,                                                      &
! IN U and V momentum fields and increments
 u, v, du_nt,dv_nt,                                                     &
! IN soil/vegetation/land surface data :
 land_index,tile_frac,canopy,fland,flandg,                              &
! IN sea/sea-ice data :
 ice_fract, di_ncat, ice_fract_ncat, k_sice, u_0, v_0,                  &
! IN cloud data :
 q,qcf,qcl,qcf_latest,qcl_latest, t,                                    &
! IN everything not covered so far :
 rho_wet_theta,pstar,                                                   &
! IN variables from BDY_LAYR
 alpha1,ashtf,bq_gb,bt_gb,dtrdz_charney_grid,rdz_charney_grid,          &
 dtrdz_u,dtrdz_v,rdz_u,rdz_v,fraca,rhokh_tile,smc,chr1p5m,              &
 resfs,z0hssi,z0mssi,cdr10m_u,cdr10m_v,cdr10m_n_u,cdr10m_n_v,           &
 cd10m_n_u,cd10m_n_v,z_theta,zh,rhokm_u,rhokm_v,                        &
 k_blend_tq, k_blend_uv,                                                &
! IN needed for new BL numerical solver
 bl_type_1,bl_type_2,                                                   &
! IN additional variables for JULES
 tile_pts,tile_index,canhc_tile,flake,wt_ext_tile,lw_down,sw_tile,      &
 alpha1_sea,alpha1_sice,ashtf_sea,ashtf_tile,resft,rhokh_sice,rhokh_sea,&
 z0h_tile,z0m_tile,chr1p5m_sice,flandg_u,flandg_v,                      &
 co2_3d, rho1, f3_at_p, uStarGBM,emis_tile,                             &
 t_soil, snow_tile,                                                     &
! INOUT diagnostics :-
 BL_diag, sf_diag,                                                      &
! INOUT data :
 dtstar_tile,dtstar_sea,dtstar_sice,ti,tstar_sice_cat,tstar_ssi,        &
 tstar_tile,tstar_sea,t_latest,q_latest,                                &
! INOUT Diagnostics started in BDY_LAYR not requiring STASH flags :
 fqw_ice,ftl_ice,fqw,fqw_tile,epot_tile,ftl,ftl_tile,rhokh,             &
 taux,tauy,taux_land,taux_ssi,tauy_land,tauy_ssi,                       &
 TScrnDcl_SSI, TScrnDcl_TILE, tStbTrans,                                &
! INOUT additional variables for JULES
 radnet_sice,olr,                                                       &
! OUT Increments to U and V momentum fields and Tl qw
 du,dv,                                                                 &
! OUT Diagnostic not requiring STASH flags :
 rhokh_mix,sea_ice_htf,                                                 &
! OUT additional variables for JULES
 ti_gb,tstar,tstar_land,tstar_sice,e_sea,h_sea,le_tile,radnet_tile,     &
 esoil_tile,surf_ht_flux,surf_ht_flux_land,surf_ht_flux_sice,           &
 surf_htf_tile,ei_tile,ecan_tile,melt_tile,                             &
 e_ssi,ei_sice,ftl_ssi,error_code,                                      &
! OUT data required elsewhere in UM system :
 ecan,ei,es,ext,snowmelt                                                &
 )

USE atm_fields_mod, ONLY: sstfrz
USE atm_fields_bounds_mod, ONLY: udims, vdims, udims_s, vdims_s, pdims, &
                                 pdims_s, tdims, tdims_l, tdims_s
USE bl_option_mod, ONLY: fric_heating, puns, pstb, on, flux_bc_opt,     &
                         interactive_fluxes
USE bl_diags_mod, ONLY: strnewbldiag

USE field_types, ONLY: fld_type_u, fld_type_v
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE jules_sea_seaice_mod, ONLY: nice, nice_use
USE jules_soil_mod, ONLY: sm_levels
USE jules_surface_types_mod, ONLY: ntype
USE level_heights_mod, ONLY:                                            &
  r_theta_levels, r_rho_levels, r_at_u, r_at_v
USE model_domain_mod, ONLY: model_type, mt_global, mt_single_column
USE planet_constants_mod, ONLY: cp, g
USE sf_diags_mod, ONLY: strnewsfdiag
USE surf_couple_implicit_mod, ONLY: surf_couple_implicit
USE swapable_field_mod, ONLY: swapable_field_pointer_type
USE timestep_mod, ONLY: timestep
USE um_parparams, ONLY: PNorth, PSouth
USE um_parvars,   ONLY: at_extremity
USE carbon_options_mod, ONLY: l_co2_interactive
USE u_to_p_mod, ONLY: u_to_p
USE v_to_p_mod, ONLY: v_to_p
#if !defined(LFRIC)
USE swap_bounds_mv_mod, ONLY: swap_bounds_mv
#endif
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

!  Inputs :-
! (a) Defining horizontal grid and subset thereof to be processed.
INTEGER, INTENT(IN) ::                                                  &
 ntiles,                                                                &
                               ! IN number of land tiles
 land_pts
                               ! IN No.of land points in whole grid.

! (b) Defining vertical grid of model atmosphere.
INTEGER, INTENT(IN) ::                                                  &
 bl_levels
                                 ! IN Max. no. of "boundary" levels

REAL, INTENT(IN) ::                                                     &
  bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                ! IN Indicator set to 1.0 if stable
                                !     b.l. diagnosed, 0.0 otherwise.
  bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! IN Indicator set to 1.0 if Sc over
                                !     stable surface layer diagnosed,
                                !     0.0 otherwise.

REAL, INTENT(IN) ::                                                     &
  alpha_cd(bl_levels),                                                  &
  du_nt(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,    &
        bl_levels),                                                     &
                          ! non-turbulent increment to u wind field
  dv_nt(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,    &
        bl_levels),                                                     &
                          ! non-turbulen increment to v wind field
  u(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,        &
      bl_levels),                                                       &
  v(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,        &
      bl_levels)

! (c) Soil/vegetation/land surface parameters (mostly constant).
INTEGER, INTENT(IN) ::                                                  &
 land_index(land_pts)      ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.

REAL, INTENT(IN) ::                                                     &
 tile_frac(land_pts,ntiles),                                            &
                              ! IN fractional coverage for each
                              !    surface tile
 canopy(land_pts,ntiles),                                               &
                              ! IN Surface/canopy water (kg/m2)
 fland(land_pts),                                                       &
                              ! IN Land fraction on land tiles.
 flandg(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! IN Land fraction on all tiles.

! (d) Sea/sea-ice data.
REAL, INTENT(IN) ::                                                     &
 ice_fract(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                 ! IN Fraction of gridbox covered by
                                 !  sea-ice (decimal fraction).
 di_ncat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice),     &
                                 ! IN "Equivalent thickness" of
                                 !   sea-ice on categories(m).
 ice_fract_ncat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,    &
                nice),                                                  &
                                 ! IN Fraction of gridbox
                                 !   covered by sea-ice category.
 k_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice),      &
                                 ! IN sea ice effective conductivity in
                                 !  sfc layer on categories (W/m2/K)
 u_0(udims%i_start:udims%i_end,udims%j_start:udims%j_end),              &
                                 ! IN W'ly component of surface
                                 !  current (m/s).
 v_0(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                 ! IN S'ly component of surface
                                 !  current (m/s).
! (e) Cloud data.
REAL, INTENT(IN) ::                                                     &
 qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,       &
     tdims_l%k_start:bl_levels),                                        &
                                     ! IN Cloud ice (kg per kg air)
 qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,       &
     tdims_l%k_start:bl_levels),                                        &
                                     ! IN Cloud liquid water
 q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,         &
   tdims_l%k_start:bl_levels),                                          &
                                     ! IN specific humidity
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),      &
                                     ! IN temperature
                        ! Latest estimates to time level n+1 values
 qcf_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                              ! IN Cloud ice (kg per kg air)
 qcl_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels) ! IN Cloud liquid water

! (f) Atmospheric + any other data not covered so far, incl control.
REAL, INTENT(IN) ::                                                     &
 pstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                 ! IN Surface pressure (Pascals).
 rho_wet_theta(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,     &
               bl_levels)
                                 ! IN wet density on theta levels

! IN variables from BDY_LAYR (that used to be local arrays)
REAL, INTENT(IN) ::                                                     &
 alpha1(land_pts,ntiles),                                               &
                             ! IN Mean gradient of saturated
!                                 specific humidity with
!                                 respect to temperature between
!                                 the bottom model layer and the
!                                 tile surfaces.
   ashtf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use), &
                                  ! IN Coefficient to calculate surface
!                                 heat flux into sea-ice.
   dtrdz_charney_grid(tdims%i_start:tdims%i_end,                        &
                      tdims%j_start:tdims%j_end,bl_levels),             &
                                ! IN -g.dt/dp for model layers.
   rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,&
                    bl_levels),                                         &
                                ! IN RDZ(,1) is the reciprocal of the
!                                 height of level 1, i.e. of the
!                                 middle of layer 1.  For K > 1,
!                                 RDZ(,K) is the reciprocal
!                                 of the vertical distance
!                                 from level K-1 to level K.
   dtrdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                          ! IN
   dtrdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                           ! IN
!                                 -g.dt/dp for model wind layers.
   bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! IN grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! IN grid-box mean buoyancy parameter
                                !     on p,T,q-levels (full levels).
   rdz_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,           &
          2:bl_levels),                                                 &
                                ! IN 1/(Z_U(K)-Z_U(K-1)) m^{-1}
   rdz_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,           &
          2:bl_levels),                                                 &
                                ! IN 1/(Z_V(K)-Z_V(K-1)) m^{-1}
   fraca(land_pts,ntiles),                                              &
                                ! IN Fraction of surface
                                !            moisture flux with only
                                !            aerodynamic resistance.
   rhokh_tile(land_pts,ntiles),                                         &
                                   ! IN
!                                 Tile surface exchange coefficients
!                                 for heat
   smc(land_pts),                                                       &
                               ! IN Soil moisture content in root depth
!                                  (kg/m2).
   chr1p5m(land_pts,ntiles),                                            &
                                ! IN Ratio of coefficients reqd for
!                                 calculation of 1.5 m T.
   resfs(land_pts,ntiles),                                              &
                              ! IN Combined soil, stomatal
!                                 and aerodynamicresistance
!                                 factor = PSIS/(1+RS/RA) for
!                                 fraction (1-FRACA)
   z0hssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
   z0mssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                ! IN Roughness lengths over sea
   cdr10m_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                ! IN Ratio of CD's reqd for calculation
   cdr10m_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),       &
                                ! IN Ratio of CD's reqd for calculation
   cdr10m_n_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),     &
                                ! IN Ratio of CD's reqd for calculation
   cdr10m_n_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),     &
                                ! IN Ratio of CD's reqd for calculation
   cd10m_n_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
                                ! IN Neutral drag coefficients at 10 m
   cd10m_n_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                                ! IN Neutral drag coefficients at 10 m
   z_theta(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           bl_levels),                                                  &
                                ! IN Height of lowest theta level.
   zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
                                ! IN BL depth (m)
   rhokm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
            bl_levels),                                                 &
                                ! IN Exchange coefficients for u
   rhokm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
            bl_levels)
                                ! IN Exchange coefficients for v

INTEGER, INTENT(IN) ::                                                  &
 k_blend_tq(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                ! IN Theta level for blending height.
 k_blend_uv(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
                                ! IN Rho level for blending height.


! IN Additional variables for screen-level diagnostics
REAL, INTENT(IN)    :: co2_3d(tdims_s%i_start:tdims_s%i_end,            &
                              tdims_s%j_start:tdims_s%j_end)
                              ! 3-D field of CO2
REAL, INTENT(IN)    :: rho1(pdims%i_start:pdims%i_end,                  &
                            pdims%j_start:pdims%j_end)
                              ! Density on lowest level
REAL, INTENT(IN)    :: f3_at_p(pdims%i_start:pdims%i_end,               &
                               pdims%j_start:pdims%j_end)
                              ! Coriolis parameter
REAL, INTENT(IN)    :: uStarGBM(pdims%i_start:pdims%i_end,              &
                                pdims%j_start:pdims%j_end)
                              ! GBM surface friction velocity

! IN additional variables for JULES
INTEGER, INTENT(IN) ::                                                  &
 tile_pts(ntype),                                                       &
                               ! IN Number of tile points.
 tile_index(land_pts,ntype)
                               ! IN Index of tile points.

REAL, INTENT(IN) ::                                                     &
 canhc_tile(land_pts,ntiles),                                           &
                              ! IN Areal heat capacity of canopy
                              !    for land tiles (J/K/m2).
 flake(land_pts,ntiles),                                                &
                              ! IN Lake fraction.
 wt_ext_tile(land_pts,sm_levels,ntiles),                                &
                              ! IN Fraction of evapotranspiration
                              !    which is extracted from each
                              !    soil layer by each tile.
 lw_down(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Surface downward LW radiation
                              !    (W/m2).
 sw_tile(land_pts,ntiles),                                              &
                              ! IN Surface net SW radiation on land
                              !    tiles (W/m2).
 alpha1_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! IN ALPHA1 for sea.
 alpha1_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             nice_use),                                                 &
                              ! IN ALPHA1 for sea-ice.
 ashtf_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                              ! IN Coefficient to calculate surface
                              !    heat flux into sea.
 ashtf_tile(land_pts,ntiles),                                           &
                              ! IN Coefficient to calculate
                              !    surface heat flux into land
                              !    tiles.
 resft(land_pts,ntiles),                                                &
                              ! IN Total resistance factor.
                              !    FRACA+(1-FRACA)*RESFS for
                              !    snow-free land, 1 for snow.
 rhokh_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
                                                           nice_use),   &
                              ! IN Surface exchange coefficients
                              !    for sea and sea-ice
 rhokh_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! IN Surface exchange coefficients
                              !    for sea
 z0h_tile(land_pts,ntiles),                                             &
                              ! IN Tile roughness lengths for heat
                              !    and moisture (m).
 z0m_tile(land_pts,ntiles),                                             &
                              ! IN Tile roughness lengths for
                              !    momentum.
 chr1p5m_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                              ! IN CHR1P5M for sea and sea-ice
                              !    (leads ignored).
 flandg_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end),         &
                              ! IN Land frac (on U-grid, with 1st
                              !    and last rows undefined or, at
                              !    present, set to "missing data")
 flandg_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                              ! IN Land frac (on V-grid, with 1st
                              !    and last rows undefined or, at
                              !    present, set to "missing data")

! Additional variables for JULES
REAL, INTENT(IN) ::                                                     &
 emis_tile(land_pts,ntiles),                                            &
                                 ! IN  Emissivity for land tiles
 snow_tile(land_pts,ntiles),                                            &
                              ! IN Snow on tiles (kg/m2).
 t_soil(land_pts,sm_levels)
                                ! IN Soil temperatures (K).

!     Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag
TYPE (strnewsfdiag), INTENT(INOUT) :: sf_diag

!  In/outs :-
REAL, INTENT(INOUT) ::                                                  &
 q_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                          ! INOUT specific humidity
 t_latest(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                                          ! INOUT temperature
 dtstar_tile(land_pts,ntiles),                                          &
                                 ! INOUT  Change in TSTAR over timestep
                                 !     for land tiles
 dtstar_sea(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! INOUT  Change is TSTAR over timestep
                                 !        for open sea
 dtstar_sice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),&
                                 ! INOUT  Change is TSTAR over timestep
                                 !        for sea-ice
 ti(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice),          &
                                ! INOUT Sea-ice category surface layer
                                !    temperature (K)
                                ! (IN only if l_sice_multilayers=T)
 tstar_sice_cat(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,    &
                nice_use),                                              &
                                ! INOUT Sea-ice sfc temperature (K).
 tstar_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                                ! INOUT Sea mean sfc temperature (K).
 tstar_tile(land_pts,ntiles),                                           &
                                ! INOUT Surface tile temperature
 tstar_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! INOUT Open sea sfc temperature (K).

!  In/Outputs :-
!-1 Diagnostic (or effectively so - includes coupled model requisites):-
!  (a) Calculated anyway (use STASH space from higher level) :-
REAL, INTENT(INOUT) ::                                                  &
fqw_ice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),  &
                              ! INOUT Surface FQW for sea-ice
ftl_ice(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,nice_use),  &
                              ! INOUT Surface FTL for sea-ice
 fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                              ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   fqw_tile(land_pts,ntiles),                                           &
                                ! INOUT surface tile moisture flux
   epot_tile(land_pts,ntiles),                                          &
                                ! INOUT surface tile potential
                                !       evaporation
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                ! INOUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   ftl_tile(land_pts,ntiles),                                           &
                                ! INOUT surface tile heat flux
   ftl_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                                ! OUT sea mean surface heat flux
   rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                ! INOUT Exchange coeffs for moisture.
   taux(udims%i_start:udims%i_end,udims%j_start:udims%j_end,            &
         bl_levels),                                                    &
                                ! INOUT W'ly component of surface wind
!                                     stress (N/sq m).(On UV-grid with
!                                     first and last rows undefined or
!                                     at present, set to missing data
   taux_land(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
                                 ! INOUT W'ly component of land sfc wind
                                 !     stress (N/sq m). (On U-grid
                                 !     with first and last rows
                                 !     undefined or, at present,
                                 !     set to missing data
   taux_ssi(udims%i_start:udims%i_end,udims%j_start:udims%j_end),       &
                                 ! INOUT W'ly compt of mean sea sfc wind
                                 !     stress (N/sq m). (On U-grid
                                 !     with first and last rows
                                 !     undefined or, at present,
                                 !     set to missing data
   tauy(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,            &
         bl_levels),                                                    &
                                ! INOUT S'ly component of surface wind
!                                     stress (N/sq m).  On UV-grid;
!                                     comments as per TAUX.
   tauy_land(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                                 ! INOUT S'ly component of land sfc wind
                                 !     stress (N/sq m).  On V-grid;
                                 !     comments as per TAUX.
   tauy_ssi(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)
                                 ! INOUT S'ly compt of mean sea sfc wind
                                 !     stress (N/sq m).  On V-grid;
                                 !     comments as per TAUX.

REAL, INTENT(INOUT) :: TScrnDcl_SSI(tdims%i_start:tdims%i_end,          &
                                    tdims%j_start:tdims%j_end)
                          !    Decoupled screen-level temperature
                          !    over sea or sea-ice
REAL, INTENT(INOUT) :: TScrnDcl_TILE(land_pts,ntiles)
                          !    Decoupled screen-level temperature
                          !    over land tiles
REAL, INTENT(INOUT) :: tStbTrans(tdims%i_start:tdims%i_end,             &
                                 tdims%j_start:tdims%j_end)
                          !    Time since transition to stable
                          !    conditions

! INOUT additional variables for JULES
REAL, INTENT(INOUT) ::                                                  &
 radnet_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             nice_use),                                                 &
                              ! INOUT Sea-ice surface net radiation.
 olr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                              ! IN    TOA - surface upward LW on
                              !       last radiation timestep
                              ! OUT   Corrected TOA outward LW
!  Outputs :-
REAL, INTENT(OUT) ::                                                    &
 du(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end,        &
      bl_levels),                                                       &
                              ! OUT BL increment to u wind field
 dv(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end,        &
      bl_levels)          ! OUT BL increment to u wind field

!-1 Diagnostic (or effectively so - includes coupled model requisites):-
!  (a) Calculated anyway (use STASH space from higher level) :-
REAL, INTENT(OUT) ::                                                    &
 rhokh_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
                              ! OUT Exchange coeffs for moisture.
! for use in tracer mixing routines
   sea_ice_htf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice)
                                ! OUT Heat flux through sea-ice
!                                     (W/m2, positive downwards).

!-2 Genuinely output, needed by other atmospheric routines :-
REAL, INTENT(OUT) ::                                                    &
 ti_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                ! OUT Sea-ice sfc layer temperature
                                !   (ice mean) (K)
 tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                ! OUT Surface temperature (K).
 tstar_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                                ! OUT   Land mean sfc temperature (K)
 tstar_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                                ! OUT Sea-ice sfc temperature (K)
                                ! (ice mean over categories)
 e_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                              ! OUT Evaporation from sea times
!                                     leads fraction. Zero over land.
!                                     (kg per square metre per sec).
   h_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                ! OUT Surface sensible heat flux over
!                                     sea times leads fraction. (W/m2)
   le_tile(land_pts,ntiles),                                            &
                                ! OUT Surface latent heat flux for
                                !       land tiles (W/m2).
   radnet_tile(land_pts,ntiles),                                        &
                                ! OUT Tile surface net radiation.
   e_ssi(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                ! OUT   Surface FQW for mean sea.
  ei_sice(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,nice_use),&
                                ! OUT   Sea-ice sumblimation
                                !       (sea mean).
   ei(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
                                ! OUT Sublimation from lying snow or
!                                   sea-ice (kg/m2/s).
   ecan(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                                ! OUT Gridbox mean evaporation from
!                                   canopy/surface store (kg/m2/s).
!                                   Zero over sea.
   es(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),             &
                                ! OUT Surface evapotranspiration from
                                !     soil moisture store (kg/m2/s).
   ext(land_pts,sm_levels),                                             &
                                ! OUT Extraction of water from each
                                !     soil layer (kg/m2/s).
   snowmelt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! OUT Snowmelt (kg/m2/s).

! OUT additional variables for JULES
REAL, INTENT(OUT) ::                                                    &
 esoil_tile(land_pts,ntiles),                                           &
                              ! OUT Evaporation from bare soil (kg/m2)
 surf_ht_flux(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                              ! OUT Net downward heat flux at surface
!                                     over land and sea-ice fraction of
!                                     gridbox (W/m2).
   surf_ht_flux_land(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end),                        &
                                ! OUT Net downward heat flux at
                                !     surface over land
                                !     fraction of gridbox (W/m2).
   surf_ht_flux_sice(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end,nice),                   &
                                ! OUT Net category downward heat flux at
                                !     surface over sea-ice
                                !     fraction of gridbox (W/m2).
   surf_htf_tile(land_pts,ntiles),                                      &
                                ! OUT Net downward surface heat flux
                                !     on tiles (W/m2)
   ei_tile(land_pts,ntiles),                                            &
                                ! OUT EI for land tiles
   ecan_tile(land_pts,ntiles),                                          &
                                ! OUT ECAN for land tiles
   melt_tile(land_pts,ntiles)
                                ! OUT Snowmelt on tiles (kg/m2/s).

INTEGER, INTENT(OUT) ::                                                 &
 error_code                   ! OUT 0 - AOK;
                              ! 1 to 7  - bad grid definition detected;
!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

! Derived local parameters.

REAL :: pnonl,p1,p2,i1,e1,e2 ! parameters for new BL solver
REAL :: sqrt2                ! SQRT(2.)

REAL ::                                                                 &
  gamma1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
  gamma2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
REAL :: r_gamma              ! level 1 of alpha_cd, for JULES
!-----------------------------------------------------------------------
!  Workspace :-
INTEGER ::                                                              &
 nblyr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                      ! number of levels in boundary
                                      !  layer, for frictional heating

REAL ::                                                                 &
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),    &
                                      ! LOCAL total water
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),    &
                                      ! LOCAL liquid water temperature
 dqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                                      ! LOCAL BL increment to q field
 dtl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                                      ! LOCAL BL increment to T field
! DU_STAR, DV_STAR: 1st stage Temporary BL incr to u, v wind
! components from new (stable) BL solver.
   du_star(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end, &
           bl_levels),                                                  &
   dv_star(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end, &
           bl_levels),                                                  &
  dqw_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                        ! OUT NT incr to qw
  dtl_nt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                        ! OUT NT incr to TL
  ct_ctq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                        ! LOCAL Coefficient in T and q
                                        !       tri-diagonal implicit
                                        !       matrix
   cq_cm_u(udims%i_start:udims%i_end,udims%j_start:udims%j_end,         &
             bl_levels),                                                &
                                        ! LOCAL Coefficient in U
                                        !       tri-diagonal implicit
                                        !       matrix
   cq_cm_v(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,         &
             bl_levels),                                                &
                                        ! LOCAL Coefficient in V
                                        !       tri-diagonal implicit
                                        !       matrix
   fric_heating_dz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end, &
                   bl_levels),                                          &
                                        !frictional heating increment
                                        !multiplied by layer thickness
   ctctq1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                                        ! LOCAL Coefficient of H*, E*
                                        !       for implicit coupling
                                        !       at level k_blend_tq
   dqw1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                        ! LOCAL Coefficient needed
                                        !       for implicit coupling
                                        !       at level k_blend_tq
   dtl1_1(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                        ! LOCAL Coefficient needed
                                        !       for implicit coupling
                                        !       at level k_blend_tq
   cq_cm_u_1(udims%i_start:udims%i_end,udims%j_start:udims%j_end),      &
                                        ! LOCAL Coefficient of taux*
                                        !       for implicit coupling
                                        !       at level k_blend_uv
   du_1(udims_s%i_start:udims_s%i_end,udims_s%j_start:udims_s%j_end),   &
                                        ! LOCAL Coefficient needed
                                        !       for implicit coupling
                                        !       at level k_blend_uv
   cq_cm_v_1(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),      &
                                        ! LOCAL Coefficient of tauy*
                                        !       for implicit coupling
                                        !       at level k_blend_uv
   dv_1(vdims_s%i_start:vdims_s%i_end,vdims_s%j_start:vdims_s%j_end),   &
                                        ! LOCAL Coefficient needed
                                        !       for implicit coupling
                                        !       at level k_blend_uv
   e_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),         &
                                        ! LOCAL FQW over mean land
   ei_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),        &
                                        ! LOCAL EI over mean land
   ftl_land(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                                        ! LOCAL FTL over mean land
   fric_heating_blyr(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end),                        &
                                        ! LOCAL Frictional heating rate in
                                        ! surface layer
   weight1, weight2, weight3,                                           &
                                        ! LOCAL Weights for f_buoy interp
   ftl_m, fqw_m, f_buoy_m,                                              &
                                        ! LOCAL Fluxes interpd to theta-levs
   dissip_mol,                                                          &
                                        ! LOCAL Molecular dissipation rate
   dissip_0_int,                                                        &
                                        ! LOCAL Integral to rho-level 1
   z_blyr,                                                              &
                                        ! LOCAL Height of surface layer
   fric_heating_inc,                                                    &
                                        ! LOCAL heating rate
   fric_heating_incv(pdims%i_start:pdims%i_end)
                                        ! LOCAL heating rate, 1D array

! Arrays below are needed for frictional dissipation heating source
REAL, ALLOCATABLE, TARGET ::                                            &
  dissip_u(:,:,:),dissip_v(:,:,:),dissip_u_p(:,:,:),dissip_v_p(:,:,:)

! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)
REAL, ALLOCATABLE :: fqw_star(:,:,:), ftl_star(:,:,:)
REAL ::                                                                 &
  taux_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end,        &
               bl_levels),                                              &
  tauy_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end,        &
               bl_levels)

! Automatic arrays for land and sea surface stress diagnostics (used
! as inputs to the ocean when the coupled model is selected)
REAL :: taux_land_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end),&
        taux_ssi_star(udims%i_start:udims%i_end,udims%j_start:udims%j_end)

REAL :: tauy_land_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end),&
        tauy_ssi_star(vdims%i_start:vdims%i_end,vdims%j_start:vdims%j_end)

LOGICAL :: l_correct
INTEGER :: i,j,k,l,n
INTEGER :: i_field ! for swap_bounds_mv

TYPE(swapable_field_pointer_type) :: fields_to_swap(2)

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='IMP_SOLVER'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)
IF (BL_diag%l_fqw) THEN
  ALLOCATE (fqw_star(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end,bl_levels))
ELSE
  ALLOCATE (fqw_star(1,1,1))
END IF
IF (BL_diag%l_ftl) THEN
  ALLOCATE (ftl_star(tdims%i_start:tdims%i_end,                         &
                     tdims%j_start:tdims%j_end,bl_levels))
ELSE
  ALLOCATE (ftl_star(1,1,1))
END IF

!----------------------------------------------------------------------
! Compute 1st stage solution (predictor).
!----------------------------------------------------------------------

! First compute the scheme coefficients for the 1st stage. Make
! coefficients dependent on the BL type for achieving better balance
! between stability-accuracy: stable BL can be strongly nonlinear and
! stiff and thus numerically unstable, so choose a large P value.
! Unstable BL are weakly nonlinear so the solver should be able to cope
! with small P.

!----------------------------------------------------------------------
sqrt2 = SQRT(2.0)

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( pdims, bl_type_1, pstb, puns, bl_type_2, sqrt2, gamma1,  &
!$OMP&         gamma2)  PRIVATE( i, j, p1, p2, pnonl, i1, e1)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    p1=bl_type_1(i,j)*pstb+(1.0-bl_type_1(i,j))*puns
    p2=bl_type_2(i,j)*pstb+(1.0-bl_type_2(i,j))*puns
    pnonl=MAX(p1,p2)
    i1 = (1.0+1.0/sqrt2)*(1.0+pnonl)
    e1 = (1.0+1.0/sqrt2)*( pnonl + (1.0/sqrt2) +                        &
                        SQRT(pnonl*(sqrt2-1.0)+0.5) )
    gamma1(i,j) = i1
    gamma2(i,j) = i1 - e1
  END DO
END DO
!$OMP END PARALLEL DO

l_correct = .FALSE.

! DEPENDS ON: BDY_IMPL3
CALL bdy_impl3 (                                                        &
! IN levels/switches
   bl_levels, l_correct,                                                &
! IN fields
   q,qcl,qcf,q_latest,qcl_latest,qcf_latest,t,t_latest,                 &
   dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                  &
   rhokh, rhokm_u, rhokm_v,                                             &
   rdz_charney_grid,rdz_u,rdz_v,gamma1,gamma2,alpha_cd,                 &
   du_nt,dv_nt,                                                         &
   k_blend_tq,k_blend_uv,                                               &
! INOUT fields
   fqw,ftl,taux,tauy,du,dv,dqw,dtl,                                     &
! OUT fields
   dqw_nt,dtl_nt,qw,tl,ct_ctq,cq_cm_u,cq_cm_v,                          &
   cq_cm_u_1,cq_cm_v_1,du_1,dv_1,                                       &
   dqw1_1,dtl1_1,ctctq1_1                                               &
    )

r_gamma=alpha_cd(1)
IF (flux_bc_opt > interactive_fluxes) THEN
  ! explicit scalar surface fluxes means surface evaporation, snow melt,
  ! etc will all be consistent with specified surface fluxes
  r_gamma=0.0
END IF ! flux_bc_opt > interactive_fluxes

CALL surf_couple_implicit(                                              &
  !Important switch
  l_correct,                                                            &
  !Forcing INTENT(IN)
  pstar, lw_down, qw, tl, u, v, u_0, v_0,                               &
  !Fluxes INTENT(IN)
  sw_tile, emis_tile,                                                   &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  !Things with _u/v get interpolated in UM
  rhokm_u, rhokm_v, gamma1, gamma2, alpha1, alpha1_sea, alpha1_sice,    &
  ashtf, ashtf_sea, ashtf_tile, du_1, dv_1, fraca, resfs, resft, rhokh, &
  rhokh_tile, rhokh_sice, rhokh_sea,                                    &
  z0hssi, z0mssi, z0h_tile, z0m_tile, chr1p5m,                          &
  chr1p5m_sice, canhc_tile, flake, tile_frac, wt_ext_tile,              &
  cdr10m_u, cdr10m_v, cdr10m_n_u,cdr10m_n_v, r_gamma,                   &
  !INOUT diagnostics
  sf_diag,                                                              &
  !Fluxes INTENT(INOUT)
  fqw_ice, ftl_ice, fqw_tile, fqw, ftl, ftl_tile,                       &
  !Misc INTENT(INOUT)
  epot_tile, dtstar_tile, dtstar_sea, dtstar_sice, radnet_sice, olr,    &
  !Fluxes INTENT(OUT)
  tstar, le_tile, radnet_tile, e_sea, h_sea, taux, tauy, ecan_tile, ei, &
  es, ext, snowmelt, melt_tile,                                         &
  ecan, ei_tile, esoil_tile, sea_ice_htf, surf_ht_flux, surf_htf_tile,  &
  !Misc INTENT(OUT)
  error_code,                                                           &
  !UM-only arguments
  !JULES ancil_info module
  !IN
  ntiles, land_pts, land_index, tile_index, tile_pts, ice_fract, sstfrz,&
  ice_fract_ncat, z_theta,                                              &
  !JULES prognostics module
!IN, tstar_tile INOUT
  canopy, smc, k_sice, t_soil, ti, snow_tile, di_ncat, tstar_tile,      &
  !JULES coastal module
  !IN
  fland, flandg,                                                        &
  !INOUT
  tstar_sea, tstar_sice_cat, tstar_ssi, taux_land, tauy_land,           &
  !2 INOUT, 3 OUT
  taux_ssi, tauy_ssi, surf_ht_flux_land, surf_ht_flux_sice, tstar_land, &
  !INOUT
  tstar_sice,                                                           &
  !JULES u_v_grid module
  !IN
  dtrdz_charney_grid,                                                   &
  !JULES switches module
  !IN
  l_co2_interactive, l_mr_physics,                                      &
  !JULES aero module
  !IN
  co2_3d,                                                               &
  !Arguments without a JULES module
!IN
  ctctq1_1,dqw1_1,dtl1_1,du_star,dv_star,cq_cm_u_1,cq_cm_v_1,           &
  flandg_u,flandg_v,                                                    &
!3 IN, 3 INOUT
  rho1, f3_at_p, uStarGBM,TScrnDcl_SSI,TScrnDcl_TILE,tStbTrans,         &
!OUT
  taux_land_star,tauy_land_star,taux_ssi_star,                          &
!OUT
  tauy_ssi_star,ei_sice,rhokh_mix, ti_gb                                &
  )

! DEPENDS ON: BDY_IMPL4
CALL bdy_impl4 (                                                        &
! IN levels, switches
   bl_levels,  l_correct,                                               &
! IN data :
   gamma1, gamma2, rhokm_u, rhokm_v,                                    &
   rdz_charney_grid,dtrdz_charney_grid,rdz_u,rdz_v,                     &
   ct_ctq,cq_cm_u,cq_cm_v,dqw_nt,dtl_nt,                                &
! INOUT data :
   qw,tl,fqw,ftl,taux,tauy,fqw_star,ftl_star,taux_star,tauy_star,       &
   du,dv,du_star,dv_star,dqw,dtl,rhokh,bl_diag,                         &
! OUT data
   t_latest,q_latest,rhokh_mix                                          &
    )
!----------------------------------------------------------------------
! Compute 2nd stage (final) solution (corrector).
!----------------------------------------------------------------------

! First compute the scheme coefficients for the 2nd stage. Make
! coefficients dependent on the BL type for achieving better balance
! between stability-accuracy: stable BL can be strongly nonlinear and
! stiff and thus numerically unstable, so choose a large P value.
! Unstable BL are weakly nonlinear so the solver should be able to cope
! with small P.

!----------------------------------------------------------------------
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( pdims, bl_type_1, pstb, puns, bl_type_2, sqrt2, gamma1,  &
!$OMP&         gamma2)  PRIVATE( i, j, p1, p2, pnonl, i1, e2)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    p1=bl_type_1(i,j)*pstb+(1.0-bl_type_1(i,j))*puns
    p2=bl_type_2(i,j)*pstb+(1.0-bl_type_2(i,j))*puns
    pnonl=MAX(p1,p2)
    i1=(1.0+1.0/sqrt2)*(1.0+pnonl)
    e2=(1.0+1.0/sqrt2)*( pnonl+(1.0/sqrt2) -                            &
                      SQRT(pnonl*(sqrt2-1.0)+0.5))
    gamma1(i,j) = i1
    gamma2(i,j) = i1 - e2
  END DO
END DO
!$OMP END PARALLEL DO

l_correct = .TRUE.
! DEPENDS ON: BDY_IMPL3
CALL bdy_impl3 (                                                        &
! IN levels/switches
   bl_levels,l_correct,                                                 &
! IN fields
   q,qcl,qcf,q_latest,qcl_latest,qcf_latest,t,t_latest,                 &
   dtrdz_charney_grid,dtrdz_u,dtrdz_v,                                  &
   rhokh, rhokm_u, rhokm_v,                                             &
   rdz_charney_grid,rdz_u,rdz_v,gamma1,gamma2,alpha_cd,                 &
   du_nt,dv_nt,                                                         &
   k_blend_tq,k_blend_uv,                                               &
! INOUT fields
   fqw,ftl,taux,tauy,du,dv,dqw,dtl,                                     &
! OUT fields
   dqw_nt,dtl_nt,qw,tl,ct_ctq,cq_cm_u,cq_cm_v,                          &
   cq_cm_u_1,cq_cm_v_1,du_1,dv_1,                                       &
   dqw1_1,dtl1_1,ctctq1_1                                               &
    )

CALL surf_couple_implicit(                                              &
  !Important switch
  l_correct,                                                            &
  !Forcing INTENT(IN)
  pstar, lw_down, qw, tl, u, v, u_0, v_0,                               &
  !Fluxes INTENT(IN)
  sw_tile, emis_tile,                                                   &
  !Misc INTENT(IN) Many of these simply come out of explicit and into here.
  !Things with _u/v get interpolated in UM
  rhokm_u, rhokm_v, gamma1, gamma2, alpha1, alpha1_sea, alpha1_sice,    &
  ashtf, ashtf_sea, ashtf_tile, du_1, dv_1, fraca, resfs, resft, rhokh, &
  rhokh_tile, rhokh_sice, rhokh_sea,                                    &
  z0hssi, z0mssi, z0h_tile, z0m_tile, chr1p5m,                          &
  chr1p5m_sice, canhc_tile, flake, tile_frac, wt_ext_tile,              &
  cdr10m_u, cdr10m_v, cdr10m_n_u, cdr10m_n_v, r_gamma,                  &
  !INOUT diagnostics
  sf_diag,                                                              &
  !Fluxes INTENT(INOUT)
  fqw_ice, ftl_ice, fqw_tile, fqw, ftl, ftl_tile,                       &
  !Misc INTENT(INOUT)
  epot_tile, dtstar_tile, dtstar_sea, dtstar_sice, radnet_sice, olr,    &
  !Fluxes INTENT(OUT)
  tstar, le_tile, radnet_tile, e_sea, h_sea, taux, tauy, ecan_tile, ei, &
  es, ext, snowmelt, melt_tile,                                         &
  ecan, ei_tile, esoil_tile, sea_ice_htf, surf_ht_flux, surf_htf_tile,  &
  !Misc INTENT(OUT)
  error_code,                                                           &
  !UM-only arguments
  !JULES ancil_info module
  !IN
  ntiles, land_pts, land_index, tile_index, tile_pts, ice_fract, sstfrz,&
  ice_fract_ncat, z_theta,                                              &
  !JULES prognostics module
!IN, tstar_tile INOUT
  canopy, smc, k_sice, t_soil, ti, snow_tile, di_ncat, tstar_tile,      &
  !JULES coastal module
  !IN
  fland, flandg,                                                        &
  !INOUT
  tstar_sea, tstar_sice_cat, tstar_ssi, taux_land, tauy_land,           &
  !2 INOUT, 3 OUT
  taux_ssi, tauy_ssi, surf_ht_flux_land, surf_ht_flux_sice, tstar_land, &
  !INOUT
  tstar_sice,                                                           &
  !JULES u_v_grid module
  !IN
  dtrdz_charney_grid,                                                   &
  !JULES switches module
  !IN
  l_co2_interactive, l_mr_physics,                                      &
  !JULES aero module
  !IN
  co2_3d,                                                               &
  !Arguments without a JULES module
!IN
  ctctq1_1,dqw1_1,dtl1_1,du_star,dv_star,cq_cm_u_1,cq_cm_v_1,           &
  flandg_u,flandg_v,                                                    &
!3 IN, 3 INOUT
  rho1, f3_at_p, uStarGBM,TScrnDcl_SSI,TScrnDcl_TILE,tStbTrans,         &
!OUT
  taux_land_star,tauy_land_star,taux_ssi_star,                          &
!OUT
  tauy_ssi_star,ei_sice,rhokh_mix, ti_gb                                &
  )

! DEPENDS ON: BDY_IMPL4
CALL bdy_impl4 (                                                        &
! IN levels, switches
   bl_levels, l_correct,                                                &
! IN data :
   gamma1, gamma2, rhokm_u, rhokm_v,                                    &
   rdz_charney_grid,dtrdz_charney_grid,rdz_u,rdz_v,                     &
   ct_ctq,cq_cm_u,cq_cm_v,dqw_nt,dtl_nt,                                &
! INOUT data :
   qw,tl,fqw,ftl,taux,tauy,fqw_star,ftl_star,taux_star,tauy_star,       &
   du,dv,du_star,dv_star,dqw,dtl,rhokh,bl_diag,                         &
! OUT data, NB these are really tl and qt on exit!
   t_latest,q_latest,rhokh_mix                                          &
    )

! Following arrays: stable & non-oscillatory solver fluxes for 1st
!                   stage (predictor)

DEALLOCATE (fqw_star)
DEALLOCATE (ftl_star)

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( pdims, ftl_land, ftl_ssi, e_land, e_ssi, ei_land,       &
!$OMP&         ntiles, tile_pts, tile_index, land_index, ftl_tile,     &
!$OMP&         tile_frac, fqw_tile, ei_tile, ftl, flandg, fqw )        &
!$OMP& PRIVATE( i, j, n, k, l )
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ftl_land(i,j)=0.0
    ftl_ssi(i,j)=0.0
    e_land(i,j)=0.0
    e_ssi(i,j)=0.0
    ei_land(i,j)=0.0
  END DO
END DO
!$OMP END DO

DO n = 1, ntiles
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tile_pts(n)
    l = tile_index(k,n)
    j=(land_index(l)-1)/pdims%i_end + 1
    i = land_index(l) - (j-1)*pdims%i_end
    ftl_land(i,j)=ftl_land(i,j) +                                       &
      ftl_tile(l,n)*tile_frac(l,n)
    e_land(i,j)=e_land(i,j) +                                           &
      fqw_tile(l,n)*tile_frac(l,n)
    ei_land(i,j)=ei_land(i,j) +                                         &
      ei_tile(l,n)*tile_frac(l,n)
  END DO
!$OMP END DO
END DO

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF (flandg(i,j) <  1.0) THEN
      ftl_ssi(i,j)=(ftl(i,j,1)-ftl_land(i,j)*flandg(i,j))               &
        /(1.0-flandg(i,j))
      e_ssi(i,j)=(fqw(i,j,1)-e_land(i,j)*flandg(i,j))                   &
        /(1.0-flandg(i,j))
    END IF
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

IF (sf_diag%sq1p5 .AND. l_mr_physics) THEN
  !----------------------------------------------------------------------
  ! Convert 1.5m mixing ratio to specific humidity
  !-----------------------------------------------
  ! Assume qcf is the same as level 1, equivalent to what is assumed
  ! in diagnostics_bl where the cloud scheme is used to separate
  ! qT into qv + qcl
  !----------------------------------------------------------------------
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      sf_diag%q1p5m(i,j)=sf_diag%q1p5m(i,j)/                            &
                         (1.0+sf_diag%q1p5m(i,j)+qcf(i,j,1))
    END DO
  END DO

  DO n = 1, ntiles
    DO k = 1, tile_pts(n)
      l = tile_index(k,n)
      j=(land_index(l)-1)/pdims%i_end + 1
      i = land_index(l) - (j-1)*pdims%i_end

      sf_diag%q1p5m_surft(l,n) = sf_diag%q1p5m_surft(l,n)/              &
                  (1.0 + sf_diag%q1p5m_surft(l,n)+qcf(i,j,1) )
    END DO
  END DO
END IF

IF ((sf_diag%l_q10m .OR. sf_diag%l_t10m) .AND. l_mr_physics) THEN
  ! Convert 10m mixing ratio to specific humidity
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      sf_diag%q10m(i,j)=sf_diag%q10m(i,j)/                              &
                         (1.0+sf_diag%q10m(i,j)+qcf(i,j,1))
    END DO
  END DO
END IF

IF ( fric_heating == on ) THEN

  ALLOCATE (dissip_u(udims_s%i_start:udims_s%i_end,                     &
                     udims_s%j_start:udims_s%j_end,bl_levels))
  ALLOCATE (dissip_v(vdims_s%i_start:vdims_s%i_end,                     &
                     vdims_s%j_start:vdims_s%j_end,bl_levels))
  ALLOCATE (dissip_u_p(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end,bl_levels))
  ALLOCATE (dissip_v_p(pdims%i_start:pdims%i_end,                       &
                       pdims%j_start:pdims%j_end,bl_levels))
  !----------------------------------------------------------------------
  ! Add heating increment from turbulence dissipation
  !--------------------------------------------------
  ! Calculate dissipation rate on theta-levels,
  ! recalling that TAU(K) is defined on theta_level(K-1)
  !----------------------------------------------------------------------

  !  Calculate as TAU * DU/DZ
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& SHARED( bl_levels, udims, dissip_u, taux, rdz_u, u, du )         &
!$OMP& PRIVATE( i, j, k )
  DO k = 1, bl_levels-1
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        dissip_u(i,j,k) = taux(i,j,k+1) * rdz_u(i,j,k+1) *              &
                    (u(i,j,k+1)+du(i,j,k+1)-u(i,j,k)-du(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  SELECT CASE (model_type)

  CASE (mt_single_column)
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        ! Also need to account for dissipation between the surface and
        ! uv-level 1 in the theta-level(K=1) heating term.
        ! Total theta-level 1 dissipation is the cell-average between
        ! the surface and rho-level 2.
        ! First calculate integral to rho-level 1:
        dissip_0_int = taux(i,j,1) * ( u(i,j,1)+du(i,j,1)-u_0(i,j) )

        ! Note integration implies multiplication by zrho1-zth0
        ! hence no division by zrho1-zth0
        ! Then take cell average up to rho-level 2:
        dissip_u(i,j,1) = ( dissip_0_int + dissip_u(i,j,1)*             &
                      (r_rho_levels(i,j,2)-r_rho_levels(i,j,1)) )       &
                   / ( r_rho_levels(i,j,2)-r_theta_levels(i,j,0) )

        ! Stress equals zero on top theta-level
        dissip_u(i,j,bl_levels) = 0.0

      END DO
    END DO
  CASE DEFAULT
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& SHARED( udims, taux, u, du, u_0, dissip_u, r_at_u,               &
!$OMP&         r_theta_levels, bl_levels )                              &
!$OMP& PRIVATE( i, j, dissip_0_int )
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        ! Also need to account for dissipation between the surface and
        ! uv-level 1 in the theta-level(K=1) heating term.
        ! Total theta-level 1 dissipation is the cell-average between
        ! the surface and rho-level 2.
        ! First calculate integral to rho-level 1:
        dissip_0_int = taux(i,j,1) * ( u(i,j,1)+du(i,j,1)-u_0(i,j) )

        ! Note integration implies multiplication by zrho1-zth0
        ! hence no division by zrho1-zth0
        ! Then take cell average up to rho-level 2:
        dissip_u(i,j,1) = ( dissip_0_int + dissip_u(i,j,1)*             &
                      (r_at_u(i,j,2)-r_at_u(i,j,1)) )                   &
                   / ( r_at_u(i,j,2)-                                   &
                 0.5*(r_theta_levels(i,j,0)+r_theta_levels(i+1,j,0)) )

        ! Stress equals zero on top theta-level
        dissip_u(i,j,bl_levels) = 0.0
      END DO
    END DO
!$OMP END PARALLEL DO

  END SELECT ! model_type

!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& SHARED( bl_levels, vdims, dissip_v, tauy, rdz_v, v, dv )         &
!$OMP& PRIVATE( i, j, k )
  DO k = 1, bl_levels-1
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        dissip_v(i,j,k) = tauy(i,j,k+1) * rdz_v(i,j,k+1) *              &
                    (v(i,j,k+1)+dv(i,j,k+1)-v(i,j,k)-dv(i,j,k))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  SELECT CASE (model_type)
  CASE DEFAULT
!$OMP  PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP& SHARED( vdims, tauy, v, dv, v_0, dissip_v, r_at_v,               &
!$OMP&         r_theta_levels, bl_levels )                              &
!$OMP& PRIVATE( i, j, dissip_0_int )
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        ! Also need to account for dissipation between the surface and
        ! uv-level 1 in the theta-level(K=1) heating term
        ! Total theta-level 1 dissipation is the cell-average between
        ! the surface and rho-level 2.
        ! First calculate integral to rho-level 1:
        dissip_0_int = tauy(i,j,1) * ( v(i,j,1)+dv(i,j,1)-v_0(i,j) )
        dissip_v(i,j,1) = ( dissip_0_int + dissip_v(i,j,1)*             &
                    (r_at_v(i,j,2)-r_at_v(i,j,1)) )                     &
                 / ( r_at_v(i,j,2)-                                     &
            0.5*(r_theta_levels(i,j,0)+r_theta_levels(i,j+1,0)) )

        ! Stress equals zero on top theta-level
        dissip_v(i,j,bl_levels) = 0.0
      END DO
    END DO
!$OMP END PARALLEL DO

  CASE (mt_single_column)
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        ! Also need to account for dissipation between the surface and
        ! uv-level 1 in the theta-level(K=1) heating term
        ! Total theta-level 1 dissipation is the cell-average between
        ! the surface and rho-level 2.
        ! First calculate integral to rho-level 1:
        dissip_0_int = tauy(i,j,1) * ( v(i,j,1)+dv(i,j,1)-v_0(i,j) )
        dissip_v(i,j,1) = ( dissip_0_int + dissip_v(i,j,1)*             &
                     (r_rho_levels(i,j,2)-r_rho_levels(i,j,1)) )        &
                   / (r_rho_levels(i,j,2)-r_theta_levels(i,j,0))

        ! Stress equals zero on top theta-level
        dissip_v(i,j,bl_levels) = 0.0

      END DO
    END DO

  END SELECT


  SELECT CASE (model_type)
#if !defined(LFRIC)
  CASE DEFAULT
    i_field = 1
    fields_to_swap(i_field) % field       => dissip_u(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_u
    fields_to_swap(i_field) % levels      =  bl_levels
    fields_to_swap(i_field) % rows        =  udims%j_end
    fields_to_swap(i_field) % vector      =  .TRUE.

    i_field = i_field + 1
    fields_to_swap(i_field) % field       => dissip_v(:,:,:)
    fields_to_swap(i_field) % field_type  =  fld_type_v
    fields_to_swap(i_field) % levels      =  bl_levels
    fields_to_swap(i_field) % rows        =  vdims%j_len
    fields_to_swap(i_field) % vector      =  .TRUE.

    CALL swap_bounds_mv(fields_to_swap, i_field, pdims%i_end,           &
         udims_s%halo_i, udims_s%halo_j)

    ! Interpolate u and v components to p grid.
!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( dissip_u, udims_s, pdims, bl_levels, at_extremity,       &
!$OMP         dissip_u_p, dissip_v, vdims_s, dissip_v_p )
    CALL u_to_p(dissip_u,                                               &
                  udims_s%i_start,udims_s%i_end,                        &
                  udims_s%j_start,udims_s%j_end,                        &
                  pdims%i_start,pdims%i_end,                            &
                  pdims%j_start,pdims%j_end,                            &
                  bl_levels, at_extremity,dissip_u_p)

    CALL v_to_p(dissip_v,                                               &
                  vdims_s%i_start,vdims_s%i_end,                        &
                  vdims_s%j_start,vdims_s%j_end,                        &
                  pdims%i_start,pdims%i_end,                            &
                  pdims%j_start,pdims%j_end,                            &
                  bl_levels, at_extremity,dissip_v_p)
!$OMP END PARALLEL

#endif
  CASE (mt_single_column)
    DO k = 1, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          dissip_u_p(i,j,k) = dissip_u(i,j,k)
          dissip_v_p(i,j,k) = dissip_v(i,j,k)
        END DO
      END DO
    END DO

  END SELECT ! model_type

  !-----------------------------------------------------------------------
  ! First, estimate molecular dissipation rate by assuming steady state
  ! subgrid KE, so that     dissip_mol = dissip_rke + f_buoy
  ! where f_buoy is the buoyancy flux interpolated to theta levels
  !-----------------------------------------------------------------------
  ! Then convert dissipation rate to heating rate,
  ! noting that dissipation rates are zero on BL_LEVELS,
  ! and redistribute within the boundary layer to account
  ! for lack of BL mixing of these increments
  !-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED( pdims, nblyr, r_rho_levels, r_theta_levels, ftl, fqw,    &
!$OMP&         g, bt_gb, bq_gb, cp, dissip_u_p, dissip_v_p, timestep,   &
!$OMP&         rho_wet_theta, fric_heating_blyr, bl_levels, zh,z_theta, &
!$OMP&         fric_heating_dz, t_latest, BL_diag )                     &
!$OMP& PRIVATE( i, j, k, weight1, weight2, weight3, ftl_m, fqw_m,       &
!$OMP&          f_buoy_m, dissip_mol, fric_heating_inc,                 &
!$OMP&          fric_heating_incv, z_blyr )
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      nblyr(i,j) = 1
      k = 1

      weight1 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,0)
      weight2 = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
      weight3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
      ftl_m = weight2 * ftl(i,j,k+1) + weight3 * ftl(i,j,k)
      fqw_m = weight2 * fqw(i,j,k+1) + weight3 * fqw(i,j,k)
      f_buoy_m = g*( bt_gb(i,j,k)*(ftl_m/cp) +                          &
                     bq_gb(i,j,k)*fqw_m )/weight1

      dissip_mol = dissip_u_p(i,j,k)+dissip_v_p(i,j,k)                  &
                            + f_buoy_m
      fric_heating_inc = MAX (0.0, timestep * dissip_mol                &
                                     / ( cp*rho_wet_theta(i,j,k) ) )

        ! Save level 1 heating increment for redistribution over
        ! boundary layer
      fric_heating_blyr(i,j) = fric_heating_inc *                       &
                       (r_rho_levels(i,j,2)-r_theta_levels(i,j,0))

    END DO
  END DO
!$OMP END DO NOWAIT

  DO k = 2, bl_levels-1
!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end

!$DIR VECTOR ALWAYS
      DO i = pdims%i_start, pdims%i_end
        weight1 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
        weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
        weight3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
        ftl_m = weight2 * ftl(i,j,k+1) + weight3 * ftl(i,j,k)
        fqw_m = weight2 * fqw(i,j,k+1) + weight3 * fqw(i,j,k)

        f_buoy_m = g*( bt_gb(i,j,k)*(ftl_m/cp) +                        &
                       bq_gb(i,j,k)*fqw_m )/weight1

        dissip_mol = dissip_u_p(i,j,k)+dissip_v_p(i,j,k)                &
                            + f_buoy_m
        fric_heating_incv(i) = MAX (0.0, timestep * dissip_mol          &
                                     / ( cp*rho_wet_theta(i,j,k) ) )
      END DO

!$DIR NOFUSION
      DO i = pdims%i_start, pdims%i_end

        IF ( z_theta(i,j,k) <= zh(i,j) ) THEN
            !------------------------------------------------------
            ! Sum increments over boundary layer to avoid
            ! adding large increments in level 1
            !------------------------------------------------------
          nblyr(i,j) = k
          fric_heating_blyr(i,j) = fric_heating_blyr(i,j) +             &
                       fric_heating_incv(i) *                           &
                       (r_rho_levels(i,j,k+1)-r_rho_levels(i,j,k))
        ELSE
          t_latest(i,j,k) = t_latest(i,j,k) + fric_heating_incv(i)
          IF (BL_diag%l_dtfric) THEN
            BL_diag%dTfric(i,j,k) = fric_heating_incv(i)
          END IF
        END IF

      END DO
    END DO
!$OMP END DO NOWAIT
  END DO

  !-----------------------------------------------------------------------
  ! Redistribute heating within the boundary layer to account
  ! for lack of BL mixing of these increments
  !-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      z_blyr = r_rho_levels(i,j,nblyr(i,j)+1)                           &
                                - r_theta_levels(i,j,0)
      fric_heating_blyr(i,j) = fric_heating_blyr(i,j) / z_blyr

      DO k = 1, nblyr(i,j)

          ! Linearly decrease heating rate across surface layer
        fric_heating_inc = 2.0 * fric_heating_blyr(i,j) *               &
                              (1.0-z_theta(i,j,k)/z_blyr)

        t_latest(i,j,k) = t_latest(i,j,k) + fric_heating_inc

        IF (BL_diag%l_dtfric) THEN
          BL_diag%dTfric(i,j,k) = fric_heating_inc
        END IF

      END DO
    END DO
  END DO
!$OMP END DO
!$OMP END PARALLEL

  DEALLOCATE (dissip_v_p)
  DEALLOCATE (dissip_u_p)
  DEALLOCATE (dissip_v)
  DEALLOCATE (dissip_u)

END IF

! Calculate pseudostresses
IF (sf_diag%suv10m_n) THEN
  sf_diag%mu10m_n(:,:) = (taux(:,:,1)) / cd10m_n_u(:,:)
  sf_diag%mv10m_n(:,:) = (tauy(:,:,1)) / cd10m_n_v(:,:)
END IF

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE imp_solver
