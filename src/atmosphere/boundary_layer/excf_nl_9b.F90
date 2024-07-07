! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  SUBROUTINE EXCF_NL ----------------------------------------------

!  Purpose: To calculate non-local exchange coefficients for
!           boundary layer subroutine KMKH.

!  Programming standard: UMDP 3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
SUBROUTINE excf_nl_9b (                                                 &
! IN levels
 bl_levels,                                                             &
! IN fields
 rdz,z_uv,z_tq,rho_mix,rho_wet_tq,rhostar_gb,v_s,fb_surf,db_top,        &
 bflux_surf,bflux_surf_sat,zeta_s,bt_top,btt_top,                       &
 df_top_over_cp,zeta_r,btc_top,db_top_cld,chi_s_top,                    &
 bt_dsct,btt_dsct,df_dsct_over_cp,zeta_r_dsc,                           &
 db_dsct_cld,chi_s_dsct,d_siems,d_siems_dsc,                            &
 db_ksurf_dry,db_ktop_dry,db_ksurf_cld,db_ktop_cld,db_dsct,             &
! INOUT fields
 coupled,cumulus,dsc,ntml,ntdsc,zh,zc,zhsc,dscdepth,zc_dsc,BL_diag,     &
! OUT fields
 rhokm,rhokh,rhokm_top,rhokh_top,f_ngstress,zdsc_base,                  &
 rhokh_top_ent,rhokh_dsct_ent,rhokh_surf_ent,nbdsc                      &
)

USE atm_fields_bounds_mod, ONLY: pdims, tdims, pdims_s
USE bl_diags_mod, ONLY: strnewbldiag
USE bl_option_mod, ONLY: ng_stress, BrownGrant97,                       &
                         BrownGrant97_limited, on, off,                 &
                         one_third, two_thirds
USE planet_constants_mod, ONLY: vkman
USE stochastic_physics_run_mod, ONLY: l_rp2, i_rp_scheme, i_rp2b,       &
                                      a_ent_1_rp, a_ent_1_rp_max,       &
                                      a_ent_1_rp_min, g1_rp
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
IMPLICIT NONE

! IN fields
INTEGER, INTENT(IN) ::                                                  &
 bl_levels
                 ! IN maximum number of boundary layer levels

REAL, INTENT(IN) ::                                                     &
 rdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                              ! IN Reciprocal of distance between
                              !    T,q-levels (m^-1). 1/RDZ(,K) is
                              !    the vertical distance from level
                              !    K-1 to level K, except that for
                              !    K=1 it is just the height of the
                              !    lowest atmospheric level.
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1), &
                              ! IN For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, Z_UV(*,K) is the height
                              !    of the k-th u,v-level (half level
                              !    k-1/2) above the surface;
                              !    for an unstaggered grid the
                              !    heights of the half-levels
                              !    0.5 to BL_LEVELS-0.5 should be
                              !    input to elements 1 to BL_LEVELS.
                              !    (1st value not used in either
                              !     case.)
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! IN For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, Z_TQ(*,K) is the height
                              !    of the k-th T,q-level (full level
                              !    k) above the surface;
                              !    code no longer works for an
                              !    unstaggered grid as Z_TQ is used
                              !    to calculate K_SURF
 rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
                              ! IN For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, RHO_MIX(*,K) is the
                              !    density at the k-th u,v-level
                              !    above the surface;
                              !    for an unstaggered grid the
                              !    densities at the layer interfaces
                              !    (half-levels) 0.5 to BL_LEVELS-0.5
                              !    should be input to elements 1 to
                              !    BL_LEVELS.
                              !    (1st value not used in either
                              !    case.)
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                              ! IN For a vertically staggered grid
                              !    with a u,v-level first above the
                              !    surface, RHO_WET_TQ(*,K) is the
                              !    density of the k-th T,q-level
                              !    above the surface;
                              !    for an unstaggered grid the
                              !    densities at the layer interfaces
                              !    (half-levels) 1.5 to BL_LEVELS+0.5
                              !    should be input to elements 1 to
                              !    BL_LEVELS.
                              !    (value for BL_LEVELS not used
                              !    in either case.)
 rhostar_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! IN Surface density (kg/m3)
 v_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                              ! IN Surface friction velocity (m/s).
 fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Buoyancy flux at the surface over
                              !    density (m^2/s^3).
 bflux_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! IN Surface buoyancy flux (kg/m/s^3).
 bflux_surf_sat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                              ! IN Saturated-air surface buoyancy
                              !    flux.
 db_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! IN Buoyancy jump across top of b.l
                              !    (m/s^2)
 df_top_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                              ! IN Radiative flux change at cloud top
                              !    divided by c_P (K.kg/m^2/s).
 bt_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! IN Buoyancy parameter at the top of
                              !    the b.l. (m/s^2/K).
 btt_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN In-cloud buoyancy parameter at
                              !    the top of the b.l. (m/s^2/K).
 btc_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Cloud fraction weighted buoyancy
                              !    parameter at the top of the b.l.
 db_top_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                ! IN In-cloud buoyancy jump at the
                              !    top of the b.l. (m/s^2).
 chi_s_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! IN Mixing fraction of just saturated
                              !    mixture at top of the b.l.
 zeta_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! IN Non-cloudy fraction of mixing
                              !    layer for surface forced
                              !    entrainment term.
 zeta_r(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! IN Non-cloudy fraction of mixing
                              !    layer for cloud top radiative
                              !    cooling entrainment term.

REAL, INTENT(IN) ::                                                     &
 db_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! IN Buoyancy parameter at the top of
                              !    the DSC layer (m/s^2/K).
 df_dsct_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                              ! IN Radiative flux change at DSC top
                              !    divided by c_P (K.kg/m^2/s).
 bt_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                               ! IN Buoyancy parameter at the top of
                               !    the DSC  (m/s^2/K).
 btt_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                               ! IN In-cloud buoyancy parameter at
                               !    the top of the DSC (m/s^2/K).
 db_dsct_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                              !     IN In-cloud buoyancy jump at the
                              !    top of the DSC (m/s^2).
 chi_s_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! IN Mixing fraction of just saturated
                              !    mixture at top of the DSC
 zeta_r_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                              ! IN Non-cloudy fraction of DSC
                              !    for cloud top radiative
                              !    cooling entrainment term.
 d_siems(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                               ! IN Siems (1990) et al. cloud-top
                               !    entr.t instab. parm
 d_siems_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                               !IN Siems (1990) et al. cloud-top
                               !   entr.t instab. parm for DSC layer
 db_ksurf_dry(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              2:bl_levels),                                             &
                                               ! IN Dry buoyancy jump
                               ! flux integral calculation (m/s2)
 db_ktop_dry(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels),                                              &
                                              ! IN Sat. buoyancy jump
                               ! flux integral calculation (m/s2)
 db_ksurf_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              2:bl_levels),                                             &
                                               ! IN Dry buoyancy jump
                               ! flux integral calculation (m/s2)
 db_ktop_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels) ! IN Sat. buoyancy jump
                               ! flux integral calculation (m/s2)
! INOUT fields
LOGICAL, INTENT(INOUT) ::                                               &
 coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! INOUT  Flag to indicate Sc layer weakly
                              !     coupled to surface (ie weakly
                              !     decoupled)
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! INOUT Flag for cumulus
 dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                 ! INOUT Flag set if decoupled stratocu layer found.

!     Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag

INTEGER, INTENT(INOUT) ::                                               &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              ! INOUT  Number of turbulently mixed
                              !     layers.
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
        ! INOUT  Top level of any decoupled turbulently mixed Sc layer.

REAL, INTENT(INOUT) ::                                                  &
 zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                              ! INOUT Boundary layer height (m).
 zc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                              ! INOUT Cloud depth (not cloud fraction
                              !    weighted) (m).
 zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                              ! INOUT Cloud-layer height (m).
 dscdepth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                              ! INOUT Decoupled cloud-layer depth (m).
 zc_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! INOUT Cloud depth (not cloud fraction
                               !    weighted) for DSC (m).
! OUT fields
INTEGER, INTENT(OUT) ::                                                 &
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
             ! OUT Bottom level of any decou turbulently mixed Sc layer.
REAL, INTENT(OUT) ::                                                    &
 rhokm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
       2:bl_levels),                                                    &
                               ! OUT Layer k-1 - to - layer k
                               !     turbulent mixing coefficient
                               !     for momentum (kg/m/s).
 rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
       2:bl_levels),                                                    &
                                       ! OUT Layer k-1 - to - layer k
                              !     turbulent mixing coefficient
                              !     for heat and moisture (kg/m/s).
 rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                              ! OUT exchange coefficient for
                              !     momentum due to top-down mixing
 rhokh_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           2:bl_levels),                                                &
                              ! OUT exchange coefficient for
                              !     heat and moisture due to top-down
                              !     mixing
 f_ngstress(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,&
            2:bl_levels),                                               &
                              ! OUT dimensionless function for
                              !     non-gradient stresses
 zdsc_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                               ! OUT Height of base of K_top in DSC
 rhokh_surf_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                    ! OUT SML surf-driven entr. KH
 rhokh_top_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                    ! OUT SML top-driven entr. KH
 rhokh_dsct_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! OUT DSC top-driven entr. KH
!  ---------------------------------------------------------------------
!    Local and other symbolic constants :-
REAL :: a_ent_1,a_ent_2,c_t,a_ent_shr,dec_thres_clear,dec_thres_cloud
REAL :: g1
INTEGER :: n_steps
PARAMETER (                                                             &
 a_ent_2=0.056,                                                         &
                              ! Entrainment parameter.
 c_t=1.0,                                                               &
                              ! Parameter in Zilitinkevich term.
 dec_thres_clear=1.0,                                                   &
                              ! Decoupling thresholds for clear and
 dec_thres_cloud=0.05,                                                  &
                              ! cloudy layers (larger makes
                              ! decoupling less likely)
 n_steps=3                                                              &
                              ! Number of steps through the mixed
                              ! layer per sweep
)

REAL :: s_m,a_ngs
PARAMETER (                                                             &
 s_m   = 1.0,                                                           &
                              ! empirical parameters in
 a_ngs = 2.7                                                            &
                              ! non-gradient stresses
)

!  Define local storage.

!  (a) Workspace.

LOGICAL ::                                                              &
 scbase(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! Flag to signal base of CML reached
 ksurf_iterate(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                    ! Flag to perform iteration to
                                    ! find top of Ksurf
 ktop_iterate(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                    ! Flag to perform iteration to find base of Ktop

REAL ::                                                                 &
 w_m_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! Turbulent velocity scale for momentum
                              ! evaluated at the top of the b.l.
 w_h_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! Turbulent velocity scale for scalars
                              ! evaluated at the top of the b.l.
 prandtl_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                              ! Turbulent Prandtl number
                              ! evaluated at the top of the b.l.
 kh_top_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                              ! Factor to ensure K_H profile is
                              ! continuous at z_h.
 km_top_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                              ! Factor to ensure K_M profile is
                              ! continuous at z_h.
 v_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                              ! velocity scale for top-down convection
 v_top_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! velocity scale for top-down convection
 kh_sct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                              ! Factor to ensure K_H profile is
                              ! continuous at z_h.
 km_sct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                              ! Factor to ensure K_M profile is
                              ! continuous at z_h.
 kh_dsct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                              ! Factor to ensure K_H profile is
                              ! continuous at z_h.
 km_dsct_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                              ! Factor to ensure K_M profile is
                              ! continuous at z_h.
 zsml_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                              ! Height of top of surf-driven K in S
 zsml_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! Height of base of top-driven K in S
 v_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                              ! Velocity scale for surface-up conve
 kh_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         2:bl_levels),                                                  &
                              ! Shape factor for non-local
                              ! turbulent mixing coefficient
 scdepth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              ! Depth of top-driven mixing in
!   wbmix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),&
                                        ! WB*DZ if were diag as mixed
!   wbend(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels)
                                        ! WB*DZ after dec diag

INTEGER ::                                                              &
 ksurf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
         ! First Theta-level above surface layer well-mixed SC layer

!  (b) Scalars.

REAL ::                                                                 &
 Prandtl,                                                               &
                  ! Turbulent Prandtl number.
 zk_uv,                                                                 &
                  ! Height above surface of u,v-level.
 zk_tq,                                                                 &
                  ! Height above surface of T,q-level.
 w_s_cubed_uv,                                                          &
                  ! Cube of free-convective velocity scale (u,v-level)
 w_s_cubed_tq,                                                          &
                  ! Cube of free-convective velocity scale (T,q-level)
 c_tke,                                                                 &
                  ! Empirical constant in tke diagnostic
 w_m_uv,                                                                &
                  ! Turbulent velocity scale for momentum (u,v-level).
 w_m_tq,                                                                &
                  ! Turbulent velocity scale for momentum (T,q-level).
 w_h_uv,                                                                &
                  ! Turbulent velocity scale for scalars (u,v-level).
 sf_term,                                                               &
                  ! Surface flux term for entrainment parametrization.
 sf_shear_term,                                                         &
                  ! Surface shear term for entrainment paramn.
 ir_term,                                                               &
                  ! Indirect radiative term for entrainment paramn.
 dr_term,                                                               &
                  ! Direct radiative term for entrainment paramn.
 evap_term,                                                             &
                  ! Evaporative term in entrainment parametrization.
 zil_corr,                                                              &
                  ! Zilitinkevich correction term in entrn. paramn.
 zeta_s_fac,                                                            &
                  ! Factor involving ZETA_S.
 zeta_r_sq,                                                             &
                  ! ZETA_R squared.
 zr,                                                                    &
                  ! Ratio ZC/ZH.
 z_pr,                                                                  &
                  ! Height above surface layer
 zh_pr,                                                                 &
                  ! Height of layer top above surface
 rhokh_ent,                                                             &
                  ! entrainment eddy viscosity
 frac_top,                                                              &
                  ! Fraction of turbulent mixing driven from the top
 factor,                                                                &
                  ! Temporary scalar
 ent_factor,                                                            &
                  ! Factor to weight entrainment by CF
 v_sum,                                                                 &
                  ! generalised turbulent velocity scale (m/s)
 alpha_t,                                                               &
                  ! Parametrized fraction of cloud-top
                  ! radiative cooling within the inversion
 dz_inv,                                                                &
                  ! Parametrizzed inversion thickness (m)
 l_rad,                                                                 &
                  ! Estimate of e-folding radiative flux
                  ! decay depth (assumed >= 25m)
 wb_cld,                                                                &
                   ! Cloud layer buoyancy flux
 wb_scld,                                                               &
                   ! Sub-cloud layer buoyancy flux
 cld_frac,                                                              &
                   ! Vertical fraction of layer containing cloud
 zb_ktop,                                                               &
                   ! height of base of K_top profile
 db_ratio,                                                              &
                   ! Temporary in ZWB0 calculation
 gamma_wbs     ! Surface layer wb gradient

INTEGER ::                                                              &
 i,j,                                                                   &
                  ! Loop counter (horizontal field index).
 k,                                                                     &
                  ! Loop counter (vertical level index).
 n_sweep,                                                               &
                  ! sweep counter
 ns           ! step counter

! 2D arrays for optimisation

REAL :: z_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                   ! inversion height (top of K profiles)
        wb_surf_int(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),&
                   ! Estimate of wb intqegrated over surface layer
        v_ktop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                   ! velocity scale for K_top profile
        z_cbase(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                   ! cloud base height
        wb_ratio(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                   ! WBN_INT/WBP_INT
        wbp_int(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                   ! Positive part of buoyancy flux integral
        wbn_int(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                   ! Negative part of buoyancy flux integral
        dec_thres(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
                   ! Local decoupling threshold
        zinv_pr(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                   ! Height of layer top above surface
        khtop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                   ! temporary KH_top in wb integration
        khsurf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                   ! temporary KH_surf in wb integration
        zwb0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                   ! height at which wb assumed to go to zero
        z_top_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
                   ! upper height limit on K profile
        z_bot_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end), &
                   ! lower height limit on K profile
        z_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! Step size (m)

INTEGER :: ntop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                   ! top level of surf-driven K profile
           kwb0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                   ! level at which wb assumed to go to zero

INTEGER :: up(pdims%i_end*pdims%j_end)    ! indicator of upward/downward sweep

LOGICAL :: status_ntop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)

INTEGER :: ksurf_min, ntop_max, ntdsc_max

! Array introduced to calculate kwb0

LOGICAL :: kstatus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
INTEGER                             :: max_ntml

! Variables for vector compression

INTEGER :: ij_len
INTEGER :: ic
INTEGER :: c_len
LOGICAL :: to_do(pdims%i_end*pdims%j_end)
INTEGER :: ind_todo(pdims%i_end*pdims%j_end)
INTEGER :: c_len_i
LOGICAL :: todo_inner(pdims%i_end*pdims%j_end)
INTEGER :: ind_todo_i(pdims%i_end*pdims%j_end)

INTEGER :: i1, j1, l, jj
INTEGER, PARAMETER :: j_block = 4
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EXCF_NL_9B'

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!     Set values of A_ENT_1, G1 and A_ENT_SHR

IF (l_rp2) THEN
  a_ent_1 = a_ent_1_rp       ! Entrainment parameter
  g1 = g1_rp                 ! Velocity scale parameter
ELSE
  a_ent_1 = 0.23             ! Entrainment parameter
  g1 = 0.85                  ! Velocity scale parameter
END IF

IF ( l_rp2 .AND. i_rp_scheme == i_rp2b ) THEN
  !
  ! Alternative calculation of a_ent_shr to extend the range
  ! over which a_ent_shr varies with a_ent_1_rp.
  !
  IF ( a_ent_1_rp_min < 0.23 .AND. a_ent_1_rp_max > 0.23 ) THEN
    a_ent_shr =                                                         &
     ( ( a_ent_1 - a_ent_1_rp_min ) / ( a_ent_1_rp_max-0.23 ) ) *       &
     ( ( 8.7*( a_ent_1 - 0.23 ) / ( a_ent_1_rp_max - a_ent_1_rp_min ) ) &
     - ( 5*( a_ent_1 - a_ent_1_rp_max ) / ( 0.23 - a_ent_1_rp_min ) ) )
  ELSE
    ! Case a_ent_1_rp_min = a_ent_1 = a_ent_1_rp_max
    a_ent_shr = 5.0 * a_ent_1 / 0.23
  END IF
ELSE
  a_ent_shr = 5.0 * a_ent_1 / 0.23   ! Entrainment parameter.
END IF

!-----------------------------------------------------------------------
! 0.  Calculate top-of-b.l. velocity scales and Prandtl number.
!-----------------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(w_s_cubed_uv, zeta_s_fac,       &
!$OMP& sf_term, sf_shear_term, zeta_r_sq, ir_term, zr, evap_term,       &
!$OMP& ent_factor, v_sum, dz_inv, l_rad, alpha_t, dr_term, zil_corr,    &
!$OMP& rhokh_ent, frac_top, i, j, k, jj, db_ratio, wb_scld,             &
!$OMP& wb_cld, zb_ktop, z_pr, cld_frac, i1, j1, l, zh_pr, zk_tq,        &
!$OMP& w_s_cubed_tq, gamma_wbs, factor, w_m_uv, w_m_tq, w_h_uv, ic,     &
!$OMP& zk_uv, Prandtl, c_tke, n_sweep, ns)

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    rhokh_surf_ent(i,j) = 0.0
    rhokh_top_ent(i,j) = 0.0
    rhokh_dsct_ent(i,j) = 0.0
    v_top(i,j) = 0.0
    v_surf(i,j) = 0.0
    v_top_dsc(i,j) = 0.0

    IF (fb_surf(i,j)  >=  0.0) THEN

      ! By definition the top of the b.l. is in the 'outer layer' so
      !  the free-convective velocity scale cubed is

      IF (coupled(i,j)) THEN
        w_s_cubed_uv = 0.25 * zhsc(i,j) * fb_surf(i,j)
      ELSE
        w_s_cubed_uv = 0.25 * zh(i,j) * fb_surf(i,j)
      END IF

      !         Turbulent velocity scale for momentum

      w_m_top(i,j) = (v_s(i,j)*v_s(i,j)*v_s(i,j) +                      &
                      w_s_cubed_uv)**one_third

      !         Turbulent Prandtl number and velocity scale for scalars

      prandtl_top(i,j) = 0.75 *                                         &
                      ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +           &
                          (4.0/25.0)*w_s_cubed_uv*w_m_top(i,j) ) /      &
                           ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +      &
                          (8.0/25.0)*w_s_cubed_uv*w_m_top(i,j) )
      w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
    ELSE
      w_m_top(i,j) = v_s(i,j)
      prandtl_top(i,j) = 0.75
      w_h_top(i,j) = w_m_top(i,j) / prandtl_top(i,j)
    END IF
  END DO
END DO
!$OMP END DO
!-----------------------------------------------------------------------
! 1.  Loop round levels; calculate the top-of-b.l. entrainment
!     mixing coefficients.
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, j_block
  DO k = 2, bl_levels
    !cdir collapse
    DO j = jj, MIN(jj+j_block-1, pdims%j_end)
      DO i = pdims%i_start, pdims%i_end

        !----------------------------------------------------------------
        ! 1.2 Calculate top-of-b.l. entrainment mixing coefficients
        !     and store b.l. top quantities for later use.
        !----------------------------------------------------------------
        !      FIRST the top of the SML (if not coupled)
        !-----------------------------------------------
        !..Initialise RHOKs: entrainment now added later for KH, in KMKHZ
        rhokh(i,j,k) = 0.0
        rhokm(i,j,k) = 0.0
        rhokh_top(i,j,k) = 0.0
        rhokm_top(i,j,k) = 0.0
        f_ngstress(i,j,k) = 0.0
        kh_surf(i,j,k) = 0.0
        IF ( k  ==  ntml(i,j)+1 .AND. .NOT. coupled(i,j) .AND.          &
             ( fb_surf(i,j)  >=  0.0 ) ) THEN
          !-----------------------------------------------------------
          ! Calculate the surface buoyancy flux term
          !-----------------------------------------------------------
          zeta_s_fac = (1.0 - zeta_s(i,j)) * (1.0 - zeta_s(i,j))
          sf_term = a_ent_1 * MAX ( 0.0 ,                               &
               ( (1.0 - zeta_s_fac) * bflux_surf(i,j)                   &
               + zeta_s_fac * bflux_surf_sat(i,j) ) )
          !-----------------------------------------------------------
          ! Calculate the surface shear term
          !-----------------------------------------------------------
          IF (fb_surf(i,j)  >=  0.0) THEN
            sf_shear_term =  a_ent_shr * v_s(i,j) * v_s(i,j) *          &
               v_s(i,j)  * rho_mix(i,j,k)  / zh(i,j)
          ELSE
            sf_shear_term = 0.0
          END IF
          !-----------------------------------------------------------
          ! Calculate the indirect radiative term
          !-----------------------------------------------------------
          zeta_r_sq = zeta_r(i,j)*zeta_r(i,j)
          ir_term = ( bt_top(i,j)*zeta_r_sq +                           &
               btt_top(i,j)*(1.0-zeta_r_sq) )                           &
               * a_ent_1 * df_top_over_cp(i,j)
          !-----------------------------------------------------------
          ! Calculate the evaporative term
          !-----------------------------------------------------------
          IF ( db_top(i,j)  >   0.0) THEN
            zr = SQRT( zc(i,j) / zh(i,j) )
            evap_term = a_ent_2 * rho_mix(i,j,k)                        &
                 * chi_s_top(i,j) * chi_s_top(i,j)                      &
                 * zr * zr * zr * db_top_cld(i,j)                       &
                 * SQRT( zh(i,j) * db_top(i,j) )
          ELSE
            evap_term = 0.0
          END IF
          !            IF (CF(I,j,K-1) >= 0.9) THEN
          ent_factor = 1.0
          !            ELSE
          !              ENT_FACTOR = EXP(-((0.90-CF(I,j,K-1))**3.0)/0.075)
          !            END IF
          !-----------------------------------------------------------
          ! Combine forcing terms to calculate the representative
          ! velocity scales
          !-----------------------------------------------------------
          v_sum   = ( (sf_term + sf_shear_term +                        &
               ent_factor*(ir_term + evap_term))                        &
               * zh(i,j) /(a_ent_1*rho_mix(i,j,k)) )**one_third
          v_top(i,j) = ( ent_factor*(ir_term+evap_term) * zh(i,j)       &
               / (a_ent_1*rho_mix(i,j,k)) )**one_third
          v_surf(i,j) = ( (sf_term) * zh(i,j)                           &
               / (a_ent_1*rho_mix(i,j,k)) )**one_third
          !-----------------------------------------------------------
          ! Calculate the direct radiative term
          !  can only calculate for DB_TOP > 0
          !-----------------------------------------------------------
          IF ( db_top(i,j)  >   0.0) THEN
            dz_inv  = MIN( v_sum * v_sum / db_top(i,j) ,100.0 )
            l_rad   = 15.0 * MAX( 1.0 , 200.0/(zc(i,j)+1.0e-14) )
            alpha_t = 1.0 - EXP(-0.5*dz_inv/l_rad)
            IF (d_siems(i,j)  >   0.0) alpha_t =                        &
                 MIN( 1.0, alpha_t + 10.0*d_siems(i,j)*(1.0-alpha_t) )
            dr_term = btc_top(i,j) * alpha_t * df_top_over_cp(i,j)
            !-----------------------------------------------------------
            ! Combine terms to calculate the entrainment
            ! mixing coefficients
            !-----------------------------------------------------------
            zil_corr = c_t * ( (sf_term + sf_shear_term + ent_factor *  &
                 (ir_term + evap_term)) /                               &
                 (rho_mix(i,j,k) * SQRT(zh(i,j))) )**two_thirds

            rhokh_ent = (sf_term + sf_shear_term                        &
                 + ent_factor*(ir_term + evap_term + dr_term))          &
                 / ((db_top(i,j) + zil_corr) * rdz(i,j,k) )

            frac_top = v_top(i,j) / (v_top(i,j)+w_h_top(i,j)+1.0e-14)

            rhokh_surf_ent(i,j) = rhokh_ent * ( 1.0 - frac_top )
            rhokh_top_ent(i,j) = rhokh_ent * frac_top

            ! APL change for C-P grid:
            !  RHOKM(I,j,K) = PRANDTL_TOP(I,j) * RHOKH_SURF_ENT(I,j)
            !               * RHO_WET_TQ(I,j,K-1) / RHO_MIX(I,j,K)
            rhokm(i,j,k) = prandtl_top(i,j) * rhokh_surf_ent(i,j)       &
                 * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))             &
                 * rho_wet_tq(i,j,k-1) / rho_mix(i,j,k)
            !  RHOKM_TOP(I,j,K) = 0.75 * RHOKH_TOP_ENT(I,j)
            !                   * RHO_WET_TQ(I,j,K-1) / RHO_MIX(I,j,K)
            rhokm_top(i,j,k) = 0.75 * rhokh_top_ent(i,j)                &
                 * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))             &
                 * rho_wet_tq(i,j,k-1) / rho_mix(i,j,k)
          END IF    ! test on DB_TOP GT 0
        END IF
        !----------------------------------------------------------------
        !      THEN the top of the DSC (if coupled use ZHSC length-scale)
        !----------------------------------------------------------------
        IF ( (k  ==  ntdsc(i,j)+1) ) THEN
          IF (coupled(i,j)) THEN
            !--------------------------------------------------------
            ! Calculate the surface buoyancy flux term
            !--------------------------------------------------------
            zeta_s_fac = (1.0 - zeta_s(i,j)) * (1.0 - zeta_s(i,j))
            sf_term = a_ent_1 * MAX ( 0.0 ,                             &
                 ( (1.0 - zeta_s_fac) * bflux_surf(i,j)                 &
                 + zeta_s_fac * bflux_surf_sat(i,j) ) )
            !--------------------------------------------------------
            ! Calculate the surface shear term
            !--------------------------------------------------------
            IF (fb_surf(i,j)  >=  0.0) THEN
              sf_shear_term = a_ent_shr * v_s(i,j)*v_s(i,j)*v_s(i,j)    &
                   * rho_mix(i,j,k)  / zhsc(i,j)
            ELSE
              sf_shear_term = 0.0
            END IF
            v_surf(i,j) = ( (sf_term) * zhsc(i,j)                       &
                 / (a_ent_1*rho_mix(i,j,k)) )**one_third
          ELSE
            sf_term = 0.0
            sf_shear_term = 0.0
          END IF
          !-----------------------------------------------------------
          ! Calculate the indirect radiative term
          !-----------------------------------------------------------
          zeta_r_sq = zeta_r_dsc(i,j)*zeta_r_dsc(i,j)
          ir_term = ( bt_dsct(i,j)*zeta_r_sq +                          &
               btt_dsct(i,j)*(1.0-zeta_r_sq) )                          &
               * a_ent_1 * df_dsct_over_cp(i,j)
          !-----------------------------------------------------------
          ! Calculate the evaporative term
          !-----------------------------------------------------------
          IF (db_dsct(i,j)  >   0.0) THEN
            zr = SQRT( zc_dsc(i,j) / dscdepth(i,j) )
            evap_term = a_ent_2 * rho_mix(i,j,k)                        &
                 * chi_s_dsct(i,j) * chi_s_dsct(i,j)                    &
                 * zr * zr * zr * db_dsct_cld(i,j)                      &
                 * SQRT( dscdepth(i,j) * db_dsct(i,j) )
          ELSE
            evap_term = 0.0
          END IF
          !             IF (CF(I,j,K-1) >= 0.9) THEN
          ent_factor = 1.0
          !             ELSE
          !              ENT_FACTOR = EXP(-((0.90-CF(I,j,K-1))**3.0)/0.075)
          !             END IF
          !-----------------------------------------------------------
          ! Combine forcing terms to calculate the representative
          ! velocity scales
          !-----------------------------------------------------------
          v_sum   = ( (sf_term + sf_shear_term +                        &
               ent_factor*(ir_term + evap_term))                        &
               * dscdepth(i,j) / (a_ent_1*rho_mix(i,j,k)) )**one_third
          v_top_dsc(i,j) =( ent_factor * (ir_term + evap_term) *        &
               dscdepth(i,j) /                                          &
               (a_ent_1*rho_mix(i,j,k)) )**one_third
          !-----------------------------------------------------------
          ! Calculate the direct radiative term
          !-----------------------------------------------------------
          IF (db_dsct(i,j)  >   0.0) THEN
            dz_inv  = MIN( v_sum*v_sum / db_dsct(i,j) ,100.0 )
            l_rad   = 15.0 * MAX( 1.0 , 200.0/(zc_dsc(i,j)+1.0) )
            alpha_t = 1.0 - EXP(-0.5*dz_inv/l_rad)
            IF (d_siems_dsc(i,j)  >   0.0) alpha_t =                    &
                 MIN( 1.0, alpha_t + 10.0*d_siems_dsc(i,j)*             &
            (1.0-alpha_t) )
            dr_term = btc_top(i,j) * alpha_t * df_dsct_over_cp(i,j)
            !----------------------------------------------------------
            ! Finally combine terms to calculate the entrainment
            ! mixing coefficients
             !----------------------------------------------------------
            zil_corr = c_t * ( (sf_term + sf_shear_term +               &
                 ent_factor*(ir_term + evap_term)) /                    &
                 (rho_mix(i,j,k) * SQRT(dscdepth(i,j))) )**two_thirds
            rhokh_dsct_ent(i,j) = ( sf_term + sf_shear_term             &
                 + ent_factor * (ir_term + evap_term + dr_term) )       &
                 / ((db_dsct(i,j) + zil_corr) * rdz(i,j,k) )

            !             RHOKM_TOP(I,j,K) = 0.75 * RHOKH_DSCT_ENT(I,j)
            !     &                   * RHO_WET_TQ(I,j,K-1) / RHO_MIX(I,j,K)
            rhokm_top(i,j,k) = 0.75 * rhokh_dsct_ent(i,j)               &
                 * rdz(i,j,k) * (z_uv(i,j,k)-z_uv(i,j,k-1))             &
                 * rho_wet_tq(i,j,k-1) / rho_mix(i,j,k)
          END IF   ! test on DB_DSCT gt 0
        END IF
      END DO
    END DO
  END DO
END DO
!$OMP END DO
!  If there is no turbulence generation in DSC layer, ignore it.

!cdir collapse
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( v_top_dsc(i,j)  <=  0.0 ) THEN
      dsc(i,j) = .FALSE.
      ntdsc(i,j) = 0
      zhsc(i,j) = 0.0
      zc_dsc(i,j) = 0.0
      dscdepth(i,j) = 0.0
      coupled(i,j) = .FALSE.
    END IF
  END DO
END DO
!$OMP END DO

! ----------------------------------------------------------------------
! 2.0 Estimate the depths of top-down and surface-up mixing.
!     These amount to diagnoses of recoupling and decoupling.
!     The K_H profiles are applied over layers such that the ratio
!        WBN_INT/WBP_INT = DEC_THRES (decoupling parameter),
!     where WBN_INT and WBP_INT are the magnitudes of the integrals of
!     the negative and positive parts, respectively, of the resulting
!     buoyancy flux profile (given by - KH * DB_FOR_FLUX).
! ----------------------------------------------------------------------
! 2.1 First test for well-mixed boundary layer
!     (ie. both KH profiles extending from cloud-top to the surface).
!     If the parcel ascent diagnosed:
!        DSC    - test for well-mixed up to ZHSC = recoupling
!        no DSC - test for well-mixed up to ZH   = decoupling
! -----------------------------------------------------------
!cdir collapse
! IAB comment this out for now as unused
!  DO j = pdims%j_start,pdims%j_end
!    DO i = pdims%i_start, pdims%i_end
!      DO k = 1, bl_levels
            ! potentially useful diagnostics
!        wbmix(i,j,k)=0.0  ! WB if were diag as well-mixed
!        wbend(i,j,k)=0.0  ! WB after dec diag
!      END DO
!    END DO
!  END DO

! Default settings

!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    zsml_top(i,j)  = zh(i,j)
    zsml_base(i,j) = 0.1 * zh(i,j)
    zdsc_base(i,j) = 0.1 * zhsc(i,j)
    z_inv(i,j) = 0.0    ! inversion height (top of K profiles)
    zwb0(i,j)  = 0.0    ! height at which WB goes to zero
    kstatus(i,j)= .TRUE.
    kwb0(i,j)  = 2
    ntop(i,j)  = -1
    ksurf(i,j) = 1
  END DO
END DO
!$OMP END DO

!This variable was taken outside the main initialisation loop to correct a
!segmentation fault seen with the Intel compiler. It is not clear why this
!helps.
!$OMP MASTER
dec_thres(:,:) = dec_thres_cloud ! use cloudy by default
!$OMP END MASTER
!$OMP BARRIER


! Find KSURF, the first theta-level above the surface layer

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, j_block
  DO k = 2, bl_levels
    !cdir collapse
    DO j = jj, MIN(jj+j_block-1, pdims%j_end)
      DO i = pdims%i_start, pdims%i_end
        IF ( z_tq(i,j,k-1)  <   0.1*zh(i,j) ) ksurf(i,j) = k
      END DO
    END DO
  END DO
END DO
!$OMP END DO


! Find kwb0, level with lowest positive cloud-free buoyancy gradient

!$OMP BARRIER
!$OMP SINGLE
max_ntml=MAXVAL(ntml)
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, j_block
  DO k = 2, max_ntml
    DO j = jj, MIN(jj+j_block-1, pdims%j_end)
      DO i = pdims%i_start, pdims%i_end
        IF (kstatus(i,j)) THEN
          IF ( (db_ksurf_dry(i,j,k) <=  0.0) .OR.                       &
               (k >= ntml(i,j)) ) THEN
            kstatus(i,j)=.FALSE.
            kwb0(i,j)=k
          END IF
        END IF
      END DO
    END DO
  END DO
END DO
!$OMP END DO

! Set flags for iteratiing wb integral to calculate depth of mixing,
! one each for KSURF and K_TOP.  Note these will be updated depending on
! what happpens on testing for a well-mixed layer in section 2.2.

!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ksurf_iterate(i,j)= .FALSE.
    ktop_iterate(i,j) = .FALSE.
    IF ( bflux_surf(i,j)  >   0.0 .AND.                                 &
           ! Otherwise: surface mixing generated by stable scheme
         .NOT. cumulus(i,j) .AND.                                       &
           ! Rule out CUMULUS layers from iteration of ZH
         ntdsc(i,j)  >   2                                              &
           ! Otherwise: layer too shallow to resolve decoupling
       ) THEN
      ksurf_iterate(i,j)= .TRUE.
    END IF
    IF ( ntdsc(i,j)  >   2 ) THEN
      ktop_iterate(i,j) = .TRUE.
    END IF
  END DO ! I
END DO ! J
!$OMP END DO

!------------------------------------------------------------
! First test buoyancy flux integral for well-mixed layer
! -----------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( bflux_surf(i,j)  >   0.0) THEN
        ! can only be coupled to an unstable SML
      IF ( ntdsc(i,j)  >   2 ) THEN
         ! well-resolved DSC layer
          ! ...test for recoupling with SML
        z_inv(i,j)  = zhsc(i,j)
        z_cbase(i,j)= z_inv(i,j) - zc_dsc(i,j)
        v_ktop(i,j) = v_top_dsc(i,j)
        ntop(i,j)   = ntdsc(i,j)
      ELSE IF ( .NOT. dsc(i,j) .AND. .NOT. cumulus(i,j) .AND.           &
                ntml(i,j)  >   2) THEN
         ! well-resolved SML
          ! ...test for decoupling
          ! Note: code can only deal with one DSC layer at a time so
          ! can't decouple SML if a DSC layer already exists.
        IF ( zc(i,j)  ==  0.0) THEN
            ! If the BL is cloud-free then use a less restrictive
            ! threshold - ideally, the parcel ascent would have
            ! found the correct BL top in this case but this test is
            ! kept to keep negative buoyancy fluxes under control
            ! (ie. DEC_THRES_CLEAR=1 ensures wbn_int < |wbp_int|)
          dec_thres(i,j) = dec_thres_clear
        END IF
        z_inv(i,j)  = zh(i,j)
        z_cbase(i,j)= z_inv(i,j) - zc(i,j)
        v_ktop(i,j) = v_top(i,j)
        ntop(i,j)   = ntml(i,j)
      END IF
        !----------------------------------------------------
        ! estimate wb integral over surface layer
        ! (and up to next theta-level, namely Z_TQ(KSURF) )
        ! assuming linear profile going to zero at ZWB0
        !----------------------------------------------------
      IF ( kwb0(i,j)  ==  ntml(i,j) ) THEN
        zwb0(i,j) = zh(i,j)
      ELSE IF ( kwb0(i,j)  ==  2 ) THEN
        zwb0(i,j) = z_uv(i,j,2)
      ELSE
        k=kwb0(i,j)
          ! now DB_KSURF_DRY(K) LE 0 and DB_KSURF_DRY(K-1) GT 0
          ! so interpolate:
        db_ratio = db_ksurf_dry(i,j,k-1)                                &
                    /(db_ksurf_dry(i,j,k-1)                             &
                    -db_ksurf_dry(i,j,k))
        db_ratio = MAX( 0.0, db_ratio )  ! trap for rounding error
        zwb0(i,j)=z_uv(i,j,k-1)+db_ratio*                               &
                    (z_uv(i,j,k)-z_uv(i,j,k-1))
      END IF
      wb_surf_int(i,j) = bflux_surf(i,j) * z_tq(i,j,ksurf(i,j)) *       &
                    ( 1.0 - z_tq(i,j,ksurf(i,j))/(2.0*zwb0(i,j)))
      wb_surf_int(i,j) = MAX( 1.0e-14, wb_surf_int(i,j) )
    ELSE
        ! only include surface layer contribution for unstable mixing
      wb_surf_int(i,j) = 1.0e-14
    END IF

    wbp_int(i,j) = wb_surf_int(i,j)  ! must be > 0
    wbn_int(i,j) = 0.0

  END DO ! I
END DO ! J
!$OMP END DO

!$OMP SINGLE
ksurf_min=MINVAL(ksurf)
ntop_max=MAXVAL(ntop)
!$OMP END SINGLE

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, j_block
  DO k = ksurf_min+1, ntop_max+1
    DO j = jj, MIN(jj+j_block-1, pdims%j_end)
      DO i = pdims%i_start, pdims%i_end

        IF ( z_inv(i,j)  >   0.0 ) THEN
          ! ----------------------------------------------
          ! worth testing layer as well-mixed to cloud-top
          ! ----------------------------------------------
          !          wbmix(i,j,ksurf(i,j)) = wb_surf_int(i,j)
          !          wbend(i,j,ksurf(i,j)) = wb_surf_int(i,j)

          zb_ktop = 0.1*z_inv(i,j)
          zinv_pr(i,j) = z_inv(i,j) - zb_ktop
          ! DB(K)is the K to K-1 difference and already
          ! integrated up to K_SURF, so start this loop at KSURF+1
          IF ( (k >= ksurf(i,j)+1) .AND.  (k <= ntop(i,j)+1) ) THEN

            khtop(i,j) = 0.0
            khsurf(i,j) = 0.0
            z_pr = z_uv(i,j,k) - zb_ktop
            IF (z_pr  >   0.0 .AND. z_pr  <   zinv_pr(i,j)) THEN
              khtop(i,j)=vkman*rho_mix(i,j,k)*g1*v_ktop(i,j)*           &
                   (( 1.0 - z_pr/zinv_pr(i,j) )**0.8)                   &
                   * z_pr * z_pr / zinv_pr(i,j)
            END IF

            z_pr = z_uv(i,j,k)
            IF ( z_pr  <   z_inv(i,j)) THEN
              !--------------------------------
              ! include surface-driven profile
              !--------------------------------
              khsurf(i,j) = vkman * rho_mix(i,j,k) *                    &
                   w_h_top(i,j)*z_pr*( 1.0 - z_pr/z_inv(i,j) )          &
                   *( 1.0 - z_pr/z_inv(i,j) )
            END IF

            IF (z_cbase(i,j)  >   z_tq(i,j,k)) THEN
              ! cloud-base above this integration range so use dry WB
              wb_scld = ( khsurf(i,j) * db_ksurf_dry(i,j,k) +           &
                   khtop(i,j) * db_ktop_dry(i,j,k) )
              wb_cld  = 0.0
            ELSE IF (z_cbase(i,j)  <   z_tq(i,j,k-1)) THEN
              ! cloud-base below this integration range so use cloudy WB
              wb_cld = ( khsurf(i,j) * db_ksurf_cld(i,j,k) +            &
                   khtop(i,j)  * db_ktop_cld(i,j,k) )
              wb_scld = 0.0
            ELSE
              ! cloud-base within this integration range
              ! so treat cloud and sub-cloud layer wb separately
              cld_frac = (z_tq(i,j,k)-z_cbase(i,j))                     &
                   /(z_tq(i,j,k)-z_tq(i,j,k-1))
              wb_cld = cld_frac                                         &
                   * ( khsurf(i,j) * db_ksurf_cld(i,j,k) +              &
                   khtop(i,j) * db_ktop_cld(i,j,k) )
              wb_scld = (1.0-cld_frac) *                                &
                   ( khsurf(i,j) * db_ksurf_dry(i,j,k) +                &
                   khtop(i,j) * db_ktop_dry(i,j,k) )
            END IF

            !            wbmix(i,j,k) = wb_cld+wb_scld
            !            wbend(i,j,k) = wb_cld+wb_scld

            IF (wb_cld  >=  0.0) THEN
              wbp_int(i,j) = wbp_int(i,j) + wb_cld
            ELSE
              wbn_int(i,j) = wbn_int(i,j) - wb_cld
            END IF
            IF (wb_scld  >=  0.0) THEN
              wbp_int(i,j) = wbp_int(i,j) + wb_scld
            ELSE
              wbn_int(i,j) = wbn_int(i,j) - wb_scld
            END IF

          END IF ! K

          wb_ratio(i,j) = wbn_int(i,j)/wbp_int(i,j)

        END IF ! Z_INV GT 0
      END DO ! I
    END DO ! J
  END DO ! K
END DO
!$OMP END DO

! Test WB_Ratio to see if layer should be well-mixed

!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( z_inv(i,j)  >   0.0 ) THEN
      IF ( wb_ratio(i,j)  <=  dec_thres(i,j) ) THEN
          ! No need to test depth of mixing any further as
          ! well-mixed layer buoyancy flux integral criteria.
          ! SML will simply stay well-mixed (and so use defaults)
        ksurf_iterate(i,j)= .FALSE.
        ktop_iterate(i,j) = .FALSE.
        IF ( dsc(i,j) ) THEN
            ! Recouple DSC with SML:
            ! move surface driven entrainment
            ! RHOKH(z_i) = rho * w_e * DZL and w_e ~ 1/DB_TOP, so:
          IF ( db_top(i,j) >  0.0 .AND. db_dsct(i,j) >  0.01 ) THEN
                                        ! can't calc Zil. term
            rhokh_surf_ent(i,j) = rhokh_surf_ent(i,j) *                 &
                   ( rho_mix(i,j,ntdsc(i,j)+1) * db_top(i,j) *          &
                     rdz(i,j,ntml(i,j)+1) ) /                           &
                   ( rho_mix(i,j,ntml(i,j)+1) * db_dsct(i,j) *          &
                                        rdz(i,j,ntdsc(i,j)+1) )
            rhokm(i,j,ntdsc(i,j)+1) = rhokm(i,j,ntml(i,j)+1) *          &
                   ( rho_wet_tq(i,j,ntdsc(i,j)) * db_top(i,j) *         &
                     rdz(i,j,ntml(i,j)+1) ) /                           &
                   ( rho_wet_tq(i,j,ntml(i,j)) * db_dsct(i,j) *         &
                                      rdz(i,j,ntdsc(i,j)+1) )
          END IF
            ! redesignate top-driven entrainment at ZHSC
            ! (ignore that calculated at ZH)
          rhokh_top_ent(i,j) = rhokh_dsct_ent(i,j)
          zh(i,j) = zhsc(i,j)
          ntml(i,j) = ntdsc(i,j)
          v_top(i,j) = v_top_dsc(i,j)
          zsml_base(i,j) = 0.1 * zh(i,j)
          zc(i,j) = zc_dsc(i,j)
          zhsc(i,j) = 0.0
          ntdsc(i,j) = 0
          v_top_dsc(i,j) = 0.0
          zdsc_base(i,j) = 0.0
          zc_dsc(i,j) = 0.0
          dsc(i,j) = .FALSE.
          cumulus(i,j) = .FALSE.
          coupled(i,j) = .FALSE.
        END IF  ! recoupled DSC layer
      ELSE   ! buoyancy flux threshold violated
          !---------------------------------
          ! Extent of mixing must be reduced
          !---------------------------------
        IF ( .NOT. cumulus(i,j) ) ksurf_iterate(i,j) = .TRUE.
        ktop_iterate(i,j)  = .TRUE.
        IF (.NOT. dsc(i,j)) THEN
            ! Set up a `COUPLED' decoupled layer
            ! Note a new ZH (and thence NTML) will be calculated by
            ! wb integral iteration.
          IF (cumulus(i,j)) zk_uv=SQRT(zh(i,j)-1000000.0)
                            ! APLTEST: shouldn't ever happen!
          dsc(i,j) = .TRUE.
          coupled(i,j) = .TRUE.  ! as full entrainment applied at ZH
          ntdsc(i,j) = ntml(i,j)
          zhsc(i,j) = zh(i,j)
          zc_dsc(i,j) = zc(i,j)
          v_top_dsc(i,j) = v_top(i,j)
            ! put all entrainment into RHOKH_TOP
          rhokh_dsct_ent(i,j) = rhokh_top_ent(i,j)                      &
                              + rhokh_surf_ent(i,j)
          rhokh_top_ent(i,j) = 0.0
          rhokh_surf_ent(i,j) = 0.0
          rhokm_top(i,j,ntml(i,j)+1) = rhokm_top(i,j,ntml(i,j)+1)       &
                                     + rhokm(i,j,ntml(i,j)+1)
          rhokm(i,j,ntml(i,j)+1) = 0.0
        END IF
      END IF   ! test on WB_RATIO LE DEC_THRES
    END IF   ! testing for well-mixed layer (Z_INV GT 0)

  END DO ! I
END DO ! J
!$OMP END DO

! ----------------------------------------------------------------------
! 2.2 Start iteration to find top of surface-driven mixing, ZSML_TOP,
!     within predetermined maximum and minimum height limits.
!     The solution is the height that gives WB_RATIO = DEC_THRES.
!     Procedure used makes 3 sweeps (up, down and up again), using
!     progressively smaller increments (Z_INC), each time stopping when
!     the buoyancy flux threshold or the height limits are reached.
!--------------------------------------------------------------------
!     If boundary layer is stable then ignore surface driven mixing.
!--------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    IF ( ksurf_iterate(i,j) ) THEN
        !-----------------------------------------------------------
        ! Mixing must extend just above surface layer
        ! (not clear precisely how to define this here: for now use
        !  K_SURF calculated from ZH)
        !-----------------------------------------------------------
      z_bot_lim(i,j)=z_uv(i,j,ksurf(i,j)+1)                             &
             + 0.1 * (z_uv(i,j,ksurf(i,j)+2)-z_uv(i,j,ksurf(i,j)+1))
        ! limit K-surf to below cloud-base
      z_top_lim(i,j)=MAX(z_bot_lim(i,j),zhsc(i,j) )

      z_cbase(i,j) = zhsc(i,j) - zc_dsc(i,j)
        !-----------------------------------------------------
        ! Initial increment to ZSML_TOP found by dividing
        ! up depth of layer within which it is allowed:
        ! Start with ZSML_TOP at lower limit and work upwards
        !-----------------------------------------------------
      z_inc(i,j)=(z_top_lim(i,j)-z_bot_lim(i,j))                        &
                   / REAL(n_steps)
      zsml_top(i,j) = z_bot_lim(i,j)

      n_sweep = 1
      wb_ratio(i,j) = dec_thres(i,j) - 1.0 ! to be < DEC_THRES
    END IF ! KSURF_ITERATE
  END DO
END DO
!$OMP END DO

!$OMP BARRIER

!$OMP MASTER
ij_len=pdims%i_end*pdims%j_end
DO i = 1, ij_len
  to_do(i)    = .FALSE.
  ind_todo(i) = i
  up(i)       = 1
END DO

c_len=ij_len

l=0
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    l=l+1
    IF (ksurf_iterate(i,j)) THEN
      to_do(l)=.TRUE.
    END IF
  END DO
END DO

ksurf_min=MINVAL(ksurf)
ntdsc_max=MAXVAL(ntdsc)

!$OMP END MASTER
!$OMP BARRIER

DO n_sweep = 1, 3
!$OMP BARRIER

!$OMP MASTER
          ! Compress to_do and ind_todo (will have new length c_len)
  ! DEPENDS ON: excfnl_cci
  CALL excfnl_cci(c_len, to_do, ind_todo)

      ! Restart inner interation with the points of outer
  c_len_i = c_len
  todo_inner(1:c_len_i) = to_do(1:c_len_i)
  ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)
!$OMP END MASTER
!$OMP BARRIER

  DO ns = 1, n_steps
!$OMP BARRIER

!$OMP MASTER
              ! Calculate active elements and compress
    ! DEPENDS ON: excfnl_compin
    CALL excfnl_compin(up, wb_ratio, dec_thres, 1,                      &
                       c_len_i, ind_todo_i, todo_inner)

!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    !cdir nodep
    DO ic = 1, c_len_i
      j1=(ind_todo_i(ic)-1)/pdims%i_end+1
      i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

      zsml_top(i1,j1)=zsml_top(i1,j1)+z_inc(i1,j1)

          ! assume wb goes to zero at ZSML_TOP
      wb_surf_int(i1,j1) =                                              &
           bflux_surf(i1,j1) * z_tq(i1,j1,ksurf(i1,j1)) *               &
           ( 1.0 - z_tq(i1,j1,ksurf(i1,j1))/                            &
                                      (2.0*zsml_top(i1,j1)) )
      wb_surf_int(i1,j1) = MAX(1.0e-14,wb_surf_int(i1,j1))
      !        wbend(i1,j1,ksurf(i1,j1)) = wb_surf_int(i1,j1)
      wbp_int(i1,j1) = wb_surf_int(i1,j1)  ! must be > 0
      wbn_int(i1,j1) = 0.0

      z_inv(i1,j1) = zsml_top(i1,j1)

    END DO ! ic c_len_i
!$OMP END DO

    !..Integrate buoyancy flux profile given this ZSML_TOP
!$OMP DO SCHEDULE(STATIC)
    DO jj = 1, c_len_i, j_block
      DO k = ksurf_min+1, ntdsc_max+1
        !cdir nodep
        DO ic = jj, MIN(jj+j_block-1,c_len_i)
          j1=(ind_todo_i(ic)-1)/pdims%i_end+1
          i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

          IF ( k  >=  ksurf(i1,j1)+1 .AND.                              &
               k  <=  ntdsc(i1,j1)+1 ) THEN

            z_pr = z_uv(i1,j1,k)
            IF (z_pr  <   z_inv(i1,j1)) THEN
              kh_surf(i1,j1,k) = w_h_top(i1,j1)                         &
                   *z_pr*(1.0-z_pr/z_inv(i1,j1) )                       &
                   *(1.0-z_pr/z_inv(i1,j1) )
            ELSE
              kh_surf(i1,j1,k) = 0.0
            END IF

            IF (z_cbase(i1,j1)  >   z_tq(i1,j1,k)) THEN
              ! cloud-base above this range so use dry WB
              wb_scld= kh_surf(i1,j1,k) * db_ksurf_dry(i1,j1,k)
              wb_cld = 0.0
            ELSE IF (z_cbase(i1,j1)  <   z_tq(i1,j1,k-1)) THEN
              ! cloud-base below this range so use cloudy WB
              wb_cld = kh_surf(i1,j1,k) * db_ksurf_cld(i1,j1,k)
              wb_scld=0.0
            ELSE
              ! cloud-base within this integration range
              ! so treat cloud and sub-cloud layer wb separately
              cld_frac = (z_tq(i1,j1,k)-z_cbase(i1,j1))                 &
                   /(z_tq(i1,j1,k)-z_tq(i1,j1,k-1))
              wb_cld  = cld_frac                                        &
                   * kh_surf(i1,j1,k)*db_ksurf_cld(i1,j1,k)
              wb_scld = (1.0-cld_frac)                                  &
                   * kh_surf(i1,j1,k)*db_ksurf_dry(i1,j1,k)
            END IF

            !            wbend(i1,j1,k) = wb_cld+wb_scld

            IF (wb_cld  >=  0.0) THEN
              wbp_int(i1,j1)= wbp_int(i1,j1) + wb_cld
            ELSE
              wbn_int(i1,j1)= wbn_int(i1,j1) - wb_cld
            END IF
            IF (wb_scld  >=  0.0) THEN
              wbp_int(i1,j1)= wbp_int(i1,j1) + wb_scld
            ELSE
              wbn_int(i1,j1) = wbn_int(i1,j1)- wb_scld
            END IF

          END IF ! K

        END DO ! ic c_len_i
      END DO ! K
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO ic = 1, c_len_i
      j1=(ind_todo_i(ic)-1)/pdims%i_end+1
      i1=ind_todo_i(ic)-(j1-1)*pdims%i_end
      wb_ratio(i1,j1) = wbn_int(i1,j1)/wbp_int(i1,j1)
    END DO ! ic c_len_i
!$OMP END DO
  END DO  ! loop stepping up through ML (N_steps)

!$OMP DO SCHEDULE(STATIC)
  !cdir nodep
  DO ic = 1, c_len
    l=ind_todo(ic)
    j1=(l-1)/pdims%i_end+1
    i1=l-(j1-1)*pdims%i_end

    !..sub-divide current Z_INC into one more part than there will be steps
    !..as there is no need to calculate WB for ZSML_TOP at a current Z_INC
    z_inc(i1,j1)= z_inc(i1,j1)/REAL(n_steps+1)

    IF ((up(l) == 1 .AND. wb_ratio(i1,j1) >= dec_thres(i1,j1)) .OR.     &
               ! hit thres while working up
        (up(l) == 0 .AND. wb_ratio(i1,j1) <= dec_thres(i1,j1))) THEN
               ! hit thres while working down
      up(l) = 1-up(l)   ! change direction of sweep
      z_inc(i1,j1)= - z_inc(i1,j1)
    ELSE IF (zsml_top(i1,j1) >= z_top_lim(i1,j1)-1.0) THEN
               ! hit upper height limit (give-or-take 1m) without
               ! reaching threshold
      to_do(ic)=.FALSE.
      zsml_top(i1,j1) = z_top_lim(i1,j1)
    ELSE IF (zsml_top(i1,j1)  <=  z_bot_lim(i1,j1)+ 1.0) THEN
               ! hit lower height limit (give-or-take 1m) without
               ! reaching threshold
      to_do(ic)=.FALSE.
      zsml_top(i1,j1) = z_bot_lim(i1,j1)
    END IF

    !..Note that if the threshold has not been passed then the next sweep
    !..continues in the same direction (but with reduced increment).

  END DO ! c_len
!$OMP END DO

END DO ! n_sweep

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ntop(i,j)=2
    status_ntop(i,j)=.TRUE.
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, j_block
  DO k = 2, bl_levels-2
    !cdir collapse
    DO j = jj, MIN(jj+j_block-1, pdims%j_end)
      DO i = pdims%i_start, pdims%i_end
        IF ( ksurf_iterate(i,j) .AND. status_ntop(i,j) ) THEN
          ! -------------
          ! find new NTML
          ! -------------
          IF (z_uv(i,j,k+1)  <   zsml_top(i,j)) THEN
            ntop(i,j) = k+1
          ELSE
            status_ntop(i,j)=.FALSE.
          END IF
          ! --------------------------------------------------------
          ! Rounding error previously found to give
          !      ZSML_TOP > Z_TOP_LIM = ZHSC
          ! Test on ZSML_TOP hitting thresholds consequently changed
          ! but also include the following failsafe tests here.
          ! --------------------------------------------------------
          ntml(i,j) = MIN( ntdsc(i,j), ntop(i,j)-1 )
          zh(i,j)   = MIN(  zhsc(i,j), zsml_top(i,j) )

        END IF  ! KSURF_ITERATE true

      END DO
    END DO
  END DO
END DO
!$OMP END DO

! ----------------------------------------------------------------------
! 2.3 Now repeat the above procedure to find the base of the
!     top-driven K profile, ZDSC_BASE.
! ----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    IF ( ktop_iterate(i,j) ) THEN

      z_bot_lim(i,j) = 0.1 * zsml_top(i,j)
          ! Ensure mixing extends slightly below base of level NTDSC
      z_top_lim(i,j) = z_uv(i,j,ntdsc(i,j)) - 0.1 *                     &
                  ( z_uv(i,j,ntdsc(i,j))-z_uv(i,j,ntdsc(i,j)-1) )
          ! Limit base of top-driven mixing to above ZH if cumulus
      IF ( cumulus(i,j) ) THEN
        z_bot_lim(i,j) = z_uv(i,j,ntml(i,j)+1)
        IF (z_top_lim(i,j) <  z_bot_lim(i,j) )                          &
             z_bot_lim(i,j) = z_top_lim(i,j)
      END IF

      z_cbase(i,j) = zhsc(i,j) - zc_dsc(i,j)

      !..Divide up depth of layer within which ZDSC_BASE is allowed
      z_inc(i,j)=(z_top_lim(i,j)-z_bot_lim(i,j))                        &
                   /REAL(n_steps)
      zdsc_base(i,j) = z_bot_lim(i,j)
                           ! will start at Z_BOT_LIM+Z_INC


      wb_ratio(i,j) = dec_thres(i,j) + 1.0 ! to be > DEC_THRES

    END IF ! KTOP_ITERATE

  END DO
END DO
!$OMP END DO

!$OMP MASTER

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    IF ( ktop_iterate(i,j)) n_sweep = 1

  END DO
END DO

ij_len=pdims%i_end*pdims%j_end
c_len=ij_len

!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
DO i = 1, ij_len
  to_do(i)    = .FALSE.
  ind_todo(i) = i
  up(i)     = 1
END DO
!$OMP END DO

!$OMP SINGLE
l=0
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    l=l+1
    IF (ktop_iterate(i,j)) THEN
      to_do(l)=.TRUE.
    END IF
  END DO
END DO
!$OMP END SINGLE

DO n_sweep = 1, 3
!$OMP BARRIER

!$OMP MASTER
      ! Compress to_do and ind_todo (will have new length c_len)
  ! DEPENDS ON: excfnl_cci
  CALL excfnl_cci(c_len, to_do, ind_todo)

      ! Restart inner interation with the points of outer
  c_len_i = c_len
  todo_inner(1:c_len_i) = to_do(1:c_len_i)
  ind_todo_i(1:c_len_i) = ind_todo(1:c_len_i)
!$OMP END MASTER
!$OMP BARRIER

  DO ns = 1, n_steps
!$OMP BARRIER

!$OMP MASTER

              ! Calculate active elements and compress
    ! DEPENDS ON: excfnl_compin
    CALL excfnl_compin(up, wb_ratio, dec_thres, 2,                      &
                       c_len_i, ind_todo_i, todo_inner)

!$OMP END MASTER
!$OMP BARRIER

!$OMP DO SCHEDULE(STATIC)
    !cdir nodep
    DO ic = 1, c_len_i
      j1=(ind_todo_i(ic)-1)/pdims%i_end+1
      i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

      zdsc_base(i1,j1) = zdsc_base(i1,j1)+z_inc(i1,j1)
      wbn_int(i1,j1) = 0.0
      wbp_int(i1,j1) = 1.0e-14
      !        wbend(i1,j1,ksurf(i1,j1)) = 0.0
      IF ( ksurf_iterate(i1,j1) .AND.                                   &
           zdsc_base(i1,j1)  <   zsml_top(i1,j1) ) THEN
            ! only include surface flux if K_SURF is included
            ! in the wb calculation and K profiles overlap
        wbp_int(i1,j1) = wb_surf_int(i1,j1)
        wbn_int(i1,j1) = 0.0
        !          wbend(i1,j1,ksurf(i1,j1)) = wb_surf_int(i1,j1)
      END IF

      zinv_pr(i1,j1) = zhsc(i1,j1)-zdsc_base(i1,j1)

    END DO ! ic c_len_i
!$OMP END DO

    !..Integrate buoyancy flux profile given this ZDSC_BASE
!$OMP DO SCHEDULE(STATIC)
    DO jj = 1, c_len_i, j_block
      DO k = ksurf_min+1, ntdsc_max+1
        !cdir nodep
        DO ic = jj, MIN(jj+j_block-1,c_len_i)
          j1=(ind_todo_i(ic)-1)/pdims%i_end+1
          i1=ind_todo_i(ic)-(j1-1)*pdims%i_end

          IF ((k >= ksurf(i1,j1)+1) .AND. (k <= ntdsc(i1,j1)+1)) THEN
            z_pr = z_uv(i1,j1,k) - zdsc_base(i1,j1)
            khtop(i1,j1) = 0.0
            IF (z_pr >   0.0 .AND. z_pr <  zinv_pr(i1,j1)) THEN
              khtop(i1,j1) = g1 * v_top_dsc(i1,j1) *                    &
                   (( 1.0 - z_pr/zinv_pr(i1,j1) )**0.8)                 &
                   * z_pr * z_pr / zinv_pr(i1,j1)
            END IF
            khsurf(i1,j1) = 0.0
            IF ( zdsc_base(i1,j1)  <   zsml_top(i1,j1) ) THEN
              ! only include K_surf if profiles overlap
              ! otherwise layers are independent
              khsurf(i1,j1) = kh_surf(i1,j1,k)
            END IF

            IF (z_cbase(i1,j1)  >   z_tq(i1,j1,k)) THEN
              ! cloud-base above this range so use dry WB
              wb_scld=( khsurf(i1,j1)*db_ksurf_dry(i1,j1,k)+            &
                   khtop(i1,j1) *db_ktop_dry(i1,j1,k) )
              wb_cld = 0.0
            ELSE IF (z_cbase(i1,j1)  <   z_tq(i1,j1,k-1)) THEN
              ! cloud-base below this range so use cloudy WB
              wb_cld = ( khsurf(i1,j1)*db_ksurf_cld(i1,j1,k) +          &
                       khtop(i1,j1)*db_ktop_cld(i1,j1,k) )
              wb_scld = 0.0
            ELSE
              ! cloud-base within this integration range
              ! so treat cloud and sub-cloud layer wb separately
              cld_frac= (z_tq(i1,j1,k)-z_cbase(i1,j1))                  &
                   /(z_tq(i1,j1,k)-z_tq(i1,j1,k-1))
              wb_cld  = cld_frac *                                      &
                   ( khsurf(i1,j1) * db_ksurf_cld(i1,j1,k) +            &
                   khtop(i1,j1)  * db_ktop_cld(i1,j1,k) )
              wb_scld = (1.0-cld_frac) *                                &
                   ( khsurf(i1,j1) * db_ksurf_dry(i1,j1,k) +            &
                 khtop(i1,j1) * db_ktop_dry(i1,j1,k) )
            END IF

            !            wbend(i1,j1,k) = wb_cld + wb_scld

            IF (wb_cld  >=  0.0) THEN
              wbp_int(i1,j1) = wbp_int(i1,j1)+wb_cld
            ELSE
              wbn_int(i1,j1) = wbn_int(i1,j1)-wb_cld
            END IF
            IF (wb_scld  >=  0.0) THEN
              wbp_int(i1,j1) = wbp_int(i1,j1)+wb_scld
            ELSE
              wbn_int(i1,j1) = wbn_int(i1,j1)-wb_scld
            END IF

          END IF ! K
        END DO ! ic c_len_i
      END DO ! K
    END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
    DO ic = 1, c_len_i
      j1=(ind_todo_i(ic)-1)/pdims%i_end+1
      i1=ind_todo_i(ic)-(j1-1)*pdims%i_end
      wb_ratio(i1,j1)=wbn_int(i1,j1)/wbp_int(i1,j1)
    END DO ! ic c_len_i
!$OMP END DO

  END DO  ! loop stepping up through ML

!$OMP DO SCHEDULE(STATIC)
  !cdir nodep
  DO ic = 1, c_len
    l=ind_todo(ic)
    j1=(l-1)/pdims%i_end+1
    i1=l-(j1-1)*pdims%i_end

    !..sub-divide current Z_INC into one more part than there will be steps
    !..as there is no need to recalculate WB at a current Z_INC
    z_inc(i1,j1)= z_inc(i1,j1)/REAL(n_steps+1)

    IF (                                                                &
       (up(l) == 1 .AND. wb_ratio(i1,j1) <= dec_thres(i1,j1)) .OR.      &
           ! hit thres while working up
       (up(l) == 0 .AND. wb_ratio(i1,j1) >= dec_thres(i1,j1))) THEN
           ! hit thres while working down
      up(l) = 1-up(l)   ! change direction of sweep
      z_inc(i1,j1)=- z_inc(i1,j1)
    ELSE IF ( zdsc_base(i1,j1) >= z_top_lim(i1,j1)-1.0 .OR.             &
              zdsc_base(i1,j1) <=  z_bot_lim(i1,j1)+1.0 ) THEN
          ! hit height limits (give-or-take 1m) without
          ! reaching threshold
      to_do(ic)=.FALSE.
    END IF

    !..Note that if the threshold has not been passed then the next sweep
    !..continues in the same direction (but with reduced increment).

  END DO ! c_len
!$OMP END DO

END DO  ! loop over sweeps

! ----------------------------------------------------------------------
! 2.4 Set depth of cloud-top driven mixing in SML when there is a DSC
!     layer above (eg. fog under Sc) to be the SML layer depth
! ----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF (cumulus(i,j) .OR. coupled(i,j) ) THEN
        ! ignore SML `cloud-top' driven mixing
      zsml_base(i,j) = zh(i,j)
      v_top(i,j)     = 0.0
    ELSE
      zsml_base(i,j) = 0.1*zh(i,j)
    END IF
  END DO  ! loop over j
END DO  ! loop over I
!$OMP END DO

!-----------------------------------------------------------------------
!     Calculate factors required to ensure that the non-local turbulent
!     mixing coefficient profiles are continuous as the entrainment
!     level is approached.
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    k=ntml(i,j)+1
    kh_top_factor(i,j) = MAX( 0.7 , 1.0 - SQRT(                         &
             rhokh_surf_ent(i,j) /                                      &
                   ( rho_mix(i,j,k)*w_h_top(i,j)*vkman*zh(i,j) ) ) )
    km_top_factor(i,j) = MAX( 0.7 , 1.0 - SQRT( rhokm(i,j,k) /          &
               ( rho_wet_tq(i,j,k-1)*w_m_top(i,j)*vkman*zh(i,j) ) ) )
    scdepth(i,j) = zh(i,j) - zsml_base(i,j)
    factor = g1 * rho_mix(i,j,k) * v_top(i,j) *vkman *scdepth(i,j)
    IF ( factor  >   0.0) THEN
      kh_sct_factor(i,j) = 1.0 -                                        &
                           ( rhokh_top_ent(i,j) / factor )**1.25
                                                        ! 1.25=1/0.8
    ELSE
      kh_sct_factor(i,j) = 1.0
    END IF
    factor = g1 * rho_wet_tq(i,j,k-1) * v_top(i,j) *                    &
                    vkman * scdepth(i,j) * 0.75
    IF ( factor  >   0.0) THEN
      km_sct_factor(i,j) = 1.0 -                                        &
                           ( rhokm_top(i,j,k) / factor )**1.25
                                                        ! 1.25=1/0.8
    ELSE
      km_sct_factor(i,j) = 1.0
    END IF

    IF (ntdsc(i,j)  >   0) THEN
      !-------------------------------------------------------------
      ! Set up factors to ensure K profile continuity at ZHSC;
      ! no need to limit size of factor as precise shape of top-down
      ! mixing profile not important.
      ! Only calculate _DSCT_FACTORs when a decoupled stratocumulus
      ! layer exists, i.e. NTDSC > 0.
      !-------------------------------------------------------------
      k=ntdsc(i,j)+1
      dscdepth(i,j) = zhsc(i,j) - zdsc_base(i,j)
      factor = g1*rho_mix(i,j,k)*v_top_dsc(i,j)*vkman*dscdepth(i,j)
      IF ( factor  >   0.0) THEN
        kh_dsct_factor(i,j) = 1.0 -                                     &
                            ( rhokh_dsct_ent(i,j) / factor )**1.25
                                                        ! 1.25=1/0.8
      ELSE
        kh_dsct_factor(i,j) = 1.0
      END IF

      factor = 0.75 * g1 * rho_wet_tq(i,j,k-1) * v_top_dsc(i,j) *       &
                           vkman * dscdepth(i,j)
      IF ( factor  >   0.0) THEN
        km_dsct_factor(i,j) = 1.0 -                                     &
                           ( rhokm_top(i,j,k) / factor )**1.25
                                                        ! 1.25=1/0.8
      ELSE
        km_dsct_factor(i,j) = 1.0
      END IF
    END IF
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 2.  Loop around levels again calculating height dependent turbulent
!     transport coefficients within the mixing layer.
!-----------------------------------------------------------------------
c_tke = 1.33/(vkman*g1)

! Reset identifiers of base of decoupled layer mixing

!$OMP DO SCHEDULE(STATIC)
!cdir collapse
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    scbase(i,j) = .FALSE.
    nbdsc(i,j)  = 0
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO jj = pdims%j_start, pdims%j_end, j_block
  DO k = 2, bl_levels
    !cdir collapse
    DO j = jj, MIN(jj+j_block-1, pdims%j_end)
      DO i = pdims%i_start, pdims%i_end

        !           Calculate the height of u,v-level above the surface
        ! *APL: z0m removed from z in K(z)
        zk_uv = z_uv(i,j,k)

        !           Calculate the height of T,q-level above the surface

        zk_tq = z_tq(i,j,k-1)

        !-------------------------------------------------------------
        ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing profiles
        ! for the surface mixed layer.
        ! This is a variation on an up-side-down version of the cubic
        ! surface-forced profiles below.  Implement between at least
        ! the top of the `surface layer' (at Z=0.1*ZH) and ZH.
        ! Note this may well include NTML+1: entrainment fluxes will
        ! be dealt with in KMKHZ.
        !-------------------------------------------------------------

        IF ( zk_uv  <   zh(i,j) .AND.                                   &
             zk_uv  >   zsml_base(i,j) ) THEN
          z_pr  = zk_uv - zsml_base(i,j)
          zh_pr = zh(i,j) - zsml_base(i,j)
          rhokh_top(i,j,k) = rho_mix(i,j,k) * v_top(i,j) * g1 *         &
               vkman * (( 1.0 - kh_sct_factor(i,j)*z_pr/zh_pr )**0.8)   &
               * z_pr * z_pr / zh_pr
        END IF
        IF ( zk_tq  <   zh(i,j) .AND.                                   &
             zk_tq  >   zsml_base(i,j) ) THEN
          z_pr = zk_tq - zsml_base(i,j)
          zh_pr = zh(i,j) - zsml_base(i,j)
          rhokm_top(i,j,k) = 0.75 * rho_wet_tq(i,j,k-1) * v_top(i,j) *  &
               g1 * vkman *                                             &
               ( ( 1.0 - km_sct_factor(i,j)*z_pr/zh_pr )**0.8 )         &
               * z_pr * z_pr / zh_pr
          ! PRANDTL=0.75
        END IF
        !-------------------------------------------------------------
        ! Add contribution to top-down mixing coefficient
        ! profiles for decoupled stratocumulus layers when
        ! one exists
        !-------------------------------------------------------------
        IF ( zk_uv  <   zhsc(i,j) .AND.                                 &
             zk_uv  >   zdsc_base(i,j) ) THEN
          IF (.NOT. scbase(i,j) ) THEN
            scbase(i,j) = .TRUE.
            ! identifies lowest layer below which there is mixing
            nbdsc(i,j) = k
          END IF
          !-----------------------------------------------------------
          ! Calculate RHOK(H/M)_TOP, top-down turbulent mixing
          ! profiles and add to any generated in the surface mixing
          ! layer.
          ! This is a variation on an up-side-down version of the
          ! cubic surface-forced profiles above.  Implement between
          ! at least the top of the `surface layer' (at Z=0.1*ZH) and
          ! ZHSC.
          !-----------------------------------------------------------
          z_pr = zk_uv - zdsc_base(i,j)
          zh_pr = zhsc(i,j) - zdsc_base(i,j)
          rhokh_top(i,j,k) = rhokh_top(i,j,k) +                         &
               rho_mix(i,j,k)*v_top_dsc(i,j)*g1*vkman*                  &
               ( ( 1.0 - kh_dsct_factor(i,j)*z_pr/zh_pr )**0.8 )        &
               * z_pr * z_pr / zh_pr
        END IF
        !-------------------------------------------------------------
        ! Now momentum
        !-------------------------------------------------------------
        IF ( zk_tq  <   zhsc(i,j) .AND.                                 &
             zk_tq  >   zdsc_base(i,j) ) THEN
          z_pr = zk_tq - zdsc_base(i,j)
          zh_pr = zhsc(i,j) - zdsc_base(i,j)
          rhokm_top(i,j,k) = rhokm_top(i,j,k) +                         &
               0.75*rho_wet_tq(i,j,k-1)*v_top_dsc(i,j)*g1*vkman*        &
               ( ( 1.0 - km_dsct_factor(i,j)*z_pr/zh_pr )**0.8 )        &
               * z_pr * z_pr / zh_pr
        END IF

        IF (BL_diag%l_tke) THEN
          ! save 1/timescale for TKE diag, completed in bdy_expl2
          IF ( zk_tq  <   zh(i,j) .AND.                                 &
               zk_tq  >   zsml_base(i,j) ) THEN
            BL_diag%tke(i,j,k) = c_tke*v_top(i,j)/zh(i,j)
          END IF
          IF ( zk_tq  <   zhsc(i,j) .AND.                               &
               zk_tq  >   zdsc_base(i,j) ) THEN
            ! save 1/timescale for TKE diag, completed in bdy_expl2
            BL_diag%tke(i,j,k) = MAX( BL_diag%tke(i,j,k),               &
                            c_tke*v_top_dsc(i,j)/dscdepth(i,j) )
          END IF
        END IF

        IF (fb_surf(i,j)  >=  0.0) THEN

          !           Calculate the free-convective scaling velocity at z(k)

          IF (zk_uv  <=  0.1*zh(i,j)) THEN

            !             Surface layer calculation

            w_s_cubed_uv = 2.5 * zk_uv * fb_surf(i,j)
          ELSE

            !             Outer layer calculation

            IF (coupled(i,j)) THEN  !  coupled and cloudy
              w_s_cubed_uv = 0.25 * zhsc(i,j) * fb_surf(i,j)
            ELSE
              w_s_cubed_uv = 0.25 * zh(i,j) * fb_surf(i,j)
            END IF
          END IF

          IF (zk_tq  <=  0.1*zh(i,j)) THEN

            !             Surface layer calculation

            w_s_cubed_tq = 2.5 * zk_tq * fb_surf(i,j)
          ELSE

            !             Outer layer calculation

            IF (coupled(i,j)) THEN  !  coupled and cloudy
              w_s_cubed_tq = 0.25 * zhsc(i,j) * fb_surf(i,j)
            ELSE
              w_s_cubed_tq = 0.25 * zh(i,j) * fb_surf(i,j)
            END IF
          END IF

          !           Turbulent velocity scale for momentum

          w_m_uv = (v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_uv)          &
               **one_third

          w_m_tq = (v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_tq)          &
               **one_third

          !           Turbulent Prandtl number and velocity scale for scalars

          Prandtl = 0.75 * ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +      &
               (4.0/25.0)*w_s_cubed_uv*w_m_uv ) /                       &
               ( v_s(i,j)*v_s(i,j)*v_s(i,j)*v_s(i,j) +                  &
               (8.0/25.0)*w_s_cubed_uv*w_m_uv )
          w_h_uv = w_m_uv / Prandtl

          IF ( zk_uv  <   zh(i,j) ) THEN
            !---------------------------------------------------------
            ! Calculate RHOKH(w_h,z/z_h)
            !---------------------------------------------------------

            rhokh(i,j,k) = rho_mix(i,j,k) * w_h_uv * vkman * zk_uv *    &
                 ( 1.0 - kh_top_factor(i,j) * ( zk_uv / zh(i,j) ) ) *   &
                 ( 1.0 - kh_top_factor(i,j) * ( zk_uv / zh(i,j) ) )

          END IF
          IF ( zk_tq  <   zh(i,j) ) THEN
            !---------------------------------------------------------
            ! Calculate RHOKM(w_m,z/z_h)
            !---------------------------------------------------------

            rhokm(i,j,k) = rho_wet_tq(i,j,k-1)*w_m_tq*vkman * zk_tq *   &
                 ( 1.0 - km_top_factor(i,j) * ( zk_tq / zh(i,j) ) ) *   &
                 ( 1.0 - km_top_factor(i,j) * ( zk_tq / zh(i,j) ) )

          END IF

          IF (BL_diag%l_tke) THEN
            ! save 1/timescale for TKE diag, completed in bdy_expl2
            IF ( zk_tq  <  zsml_top(i,j) ) THEN
              BL_diag%tke(i,j,k) = MAX( BL_diag%tke(i,j,k),             &
                                        c_tke*w_m_tq/zh(i,j) )
            END IF
          END IF

        END IF
      END DO
    END DO
  END DO
END DO
!$OMP END DO

IF ( ng_stress  ==  BrownGrant97 .OR.                                   &
     ng_stress  ==  BrownGrant97_limited ) THEN

!$OMP DO SCHEDULE(STATIC)
  DO jj = pdims%j_start, pdims%j_end, j_block
    DO k = 2, bl_levels
      !cdir collapse
      DO j = jj, MIN(jj+j_block-1, pdims%j_end)
        DO i = pdims%i_start, pdims%i_end
          zk_tq = z_tq(i,j,k-1)   ! stresses are calc on theta-levs
          IF ( fb_surf(i,j)  >   0.0 .AND. zk_tq  <   zh(i,j) ) THEN
            !---------------------------------------------------------
            ! Calculate non-gradient stress function
            ! (Brown and Grant 1997)
            ! Shape function chosen such that non-gradient stress
            ! goes to zero at 0.1*ZH and ZH
            !---------------------------------------------------------
            IF ( zk_tq  >   0.1*zh(i,j) ) THEN
              z_pr = zk_tq - 0.1*zh(i,j)
              zh_pr = 0.9*zh(i,j)

              ! Outer layer calculation

              IF (coupled(i,j)) THEN  !  coupled and cloudy
                w_s_cubed_tq = 0.25 * zhsc(i,j) * fb_surf(i,j)
              ELSE
                w_s_cubed_tq = 0.25 * zh(i,j) * fb_surf(i,j)
              END IF

              ! 4*W_S_CUBED_TQ = the convective boundary layer
              ! velocity scale cubed
              ! V_S = the friction velocity

              f_ngstress(i,j,k) =(rho_wet_tq(i,j,k-1)/rhostar_gb(i,j))  &
                 * s_m * ( a_ngs * 4.0*w_s_cubed_tq/                    &
                 (v_s(i,j)*v_s(i,j)*v_s(i,j) + w_s_cubed_tq*4.0*0.6 ) ) &
                   * ( z_pr / zh_pr ) * ( 1.0 -  ( z_pr / zh_pr ) ) *   &
                   ( 1.0 -  ( z_pr / zh_pr ) )

            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
!$OMP END DO
END IF

!$OMP END PARALLEL

! IAB commented out for now as unused
!  DO k = 1, bl_levels
!    DO j = pdims%j_start,pdims%j_end
!      DO i = pdims%i_start, pdims%i_end
          ! convert to m2/s-3
!        IF (k >  ksurf(i,j)) THEN
!          wbmix(i,j,k)=wbmix(i,j,k)*rdz(i,j,k)
!          wbend(i,j,k)=wbend(i,j,k)*rdz(i,j,k)
!        ELSE
!          gamma_wbs = ( (wbmix(i,j,ksurf(i,j))/z_tq(i,j,ksurf(i,j)))    &
!                        - bflux_surf(i,j)  )*2.0/z_tq(i,j,ksurf(i,j))
!          wbmix(i,j,k) = bflux_surf(i,j) + gamma_wbs*z_uv(i,j,k)

!          gamma_wbs = ( (wbend(i,j,ksurf(i,j))/z_tq(i,j,ksurf(i,j)))    &
!                        - bflux_surf(i,j)  )*2.0/z_tq(i,j,ksurf(i,j))
!          wbend(i,j,k) =  bflux_surf(i,j) + gamma_wbs*z_uv(i,j,k)
!        END IF
!      END DO
!    END DO
!  END DO

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE excf_nl_9b
