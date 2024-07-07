! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: Calculate the explicit turbulent fluxes of heat, moisture
!           and momentum between atmospheric levels
!           within the boundary layer, and/or the effects of these
!           fluxes on the primary model variables.

!  Programming standard : UMDP 3

!  Documentation: UMDP 24.

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
MODULE bdy_expl2_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'BDY_EXPL2_MOD'
CONTAINS

SUBROUTINE bdy_expl2 (                                                  &
! IN values defining vertical grid of model atmosphere :
 bl_levels,p_theta_levels,land_pts,land_index,                          &
! IN U, V and W momentum fields.
 u_p,v_p,u_0_px,v_0_px,                                                 &
! IN from other part of explicit boundary layer code
 rho_mix,rho_wet_tq,rho_mix_tq,dzl_charney,rdz,rdz_charney_grid,        &
 z_tq,z_uv,rhostar,bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt,  &
 recip_l_mo_sea,flandg,rib_gb, sil_orog_land, z0m_eff_gb,               &
! IN cloud/moisture data :
 cf_bulk,q,qcf,qcl,t,qw,tl,                                             &
! IN everything not covered so far :
 rad_hr,micro_tends,fb_surf,u_s,pstar,tstar,h_blend_orog,               &
 zh_prev,zhpar,z_lcl,ho2r2_orog,sd_orog,                                &
! SCM Diagnostics (dummy values in full UM) & stash diagnostics
 nSCMDpkgs,L_SCMDiags,BL_diag,                                          &
! INOUT variables
 zh,dzh,ntml,ntpar,l_shallow,cumulus,fqw,ftl,rhokh,rhokm,w,etadot,      &
 t1_sd,q1_sd,                                                           &
! OUT new variables for message passing
 tau_fd_x, tau_fd_y, f_ngstress,                                        &
! OUT Diagnostic not requiring STASH flags :
 zht,zhnl,shallowc,cu_over_orog,                                        &
 bl_type_1,bl_type_2,bl_type_3,bl_type_4,bl_type_5,bl_type_6,bl_type_7, &
! Out data for turbulent generation of mixed-phase cloud:
 bl_w_var,                                                              &
! OUT data required for tracer mixing :
 kent, we_lim, t_frac, zrzi, kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc,&
! OUT data required elsewhere in UM system :
 zhsc,ntdsc,nbdsc,wstar,wthvs,uw0,vw0,rhcpt                             &
 )

USE atm_fields_bounds_mod, ONLY: pdims, tdims, wdims, tdims_l,          &
                                 pdims_s,ScmRowLen,ScmRow, tdims_s
USE bl_diags_mod, ONLY: strnewbldiag
USE bl_option_mod, ONLY:                                                &
    i_bl_vn, i_bl_vn_9b, i_bl_vn_9c,                                    &
    off, max_t_grad, a_grad_adj, sg_orog_mixing, l_use_surf_in_ri,      &
    h_scale, t_drain, idyndiag, DynDiag_ZL,                             &
    DynDiag_ZL_corrn, DynDiag_ZL_CuOnly, DynDiag_Ribased,               &
    RiCrit_sharp, zhloc_depth_fac, non_local_bl, on, l_full_lambdas,    &
    nl_bl_levels, local_fa, free_trop_layers, to_sharp_across_1km,      &
    sbl_op, equilibrium_sbl, one_third, blending_option,                &
    blend_except_cu, sg_shear, sg_shear_enh_lambda, max_tke
USE cloud_inputs_mod, ONLY: i_rhcpt, forced_cu
USE cderived_mod, ONLY: delta_lambda, delta_phi
USE cv_run_mod, ONLY: l_param_conv, l_wvar_for_conv
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE jules_surface_mod, ONLY: formdrag, explicit_stress
USE level_heights_mod, ONLY: r_theta_levels, r_rho_levels
USE model_domain_mod, ONLY: model_type, mt_single_column
USE mphys_inputs_mod, ONLY: l_subgrid_qcl_mp
USE pc2_constants_mod, ONLY: rhcpt_tke_based
USE planet_constants_mod, ONLY: vkman, grcp, pref, kappa, g
USE s_scmop_mod,   ONLY: default_streams,                               &
                         t_inst, t_avg, d_bl, d_sl, d_point, scmdiag_bl
USE scmoutput_mod, ONLY: scmoutput
USE science_fixes_mod, ONLY: l_fix_dyndiag, l_fix_zh
USE stochastic_physics_run_mod, ONLY: l_rp2, par_mezcla
USE trignometric_mod, ONLY: cos_theta_latitude
USE turb_diff_ctl_mod, ONLY:                                            &
    visc_m, visc_h, rneutml_sq, max_diff, delta_smag
USE turb_diff_mod, ONLY:                                                &
    l_subfilter_vert, l_subfilter_horiz, mix_factor,                    &
    l_blend_isotropic, turb_startlev_vert, turb_endlev_vert

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    qsat_wat_new     => qsat_wat,                       &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_bl !Currently defaults to FALSE

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

!$ USE omp_lib

USE btq_int_mod, ONLY: btq_int
USE ex_coef_mod, ONLY: ex_coef
USE ex_flux_tq_mod, ONLY: ex_flux_tq
USE kmkh_mod, ONLY: kmkh
USE kmkhz_9b_mod, ONLY: kmkhz_9b
USE kmkhz_9c_mod, ONLY: kmkhz_9c
IMPLICIT NONE

!  Inputs :-
INTEGER, INTENT(IN) ::                                                  &
 land_pts,                                                              &
                             ! No.of land points in whole grid.
 bl_levels
                             ! IN Max. no. of "boundary" levels

!     Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag

REAL, INTENT(IN) ::                                                     &
  p_theta_levels(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 0:bl_levels+1),                                        &
 rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels+1),                                                  &
                                 ! IN density on UV (ie. rho) levels;
                                 !    used in RHOKH so dry density if
                                 !    L_mr_physics is true
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                                 ! IN density on TQ (ie. theta) levels;
                                 !    used in RHOKM so wet density
 rho_mix_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                                 ! IN density on TQ (ie. theta) levels;
                                 !    used in non-turb flux integration
                                 !    so dry density if L_mr_physics is true
 dzl_charney(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             bl_levels),                                                &
                                 ! IN DZL(,K) is depth in m of theta
                                 !    level K, i.e. distance from
                                 !    boundary K-1/2 to boundary K+1/2
 rdz( pdims_s%i_start:pdims_s%i_end,                                    &
      pdims_s%j_start:pdims_s%j_end, bl_levels ),                       &
                                 ! IN RDZ(,1) is the reciprocal of
                                 !    the height of level 1, i.e. of
                                 !    the middle of layer 1.  For
                                 !    K > 1, RDZ(,K) is the
                                 !    reciprocal of the vertical
                                 !    distance from level K-1 to
                                 !    level K.
 rdz_charney_grid(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  bl_levels),                                           &
                                 ! IN RDZ(,1) is the reciprocal of
                                 !       the height of level 1,
                                 !       i.e. of the middle of layer 1
                                 !       For K > 1, RDZ(,K) is the
                                 !       reciprocal of the vertical
                                 !       distance from level K-1 to
                                 !       level K.
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                 ! IN Z_tq(*,K) is height of full
                                 !    level k.
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1), &
                                  ! OUT Z_uv(*,K) is height of half
                                  ! level k-1/2.
 rhostar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! IN Surface air density
 u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                                 ! IN U on P-grid.
 v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                                 ! IN V on P-grid.
 bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),     &
                                 ! IN A buoyancy parameter for clear
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),     &
                                 ! IN A buoyancy parameter for clear
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bt_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        bl_levels),                                                     &
                                 ! IN A buoyancy parameter for cloudy
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bq_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        bl_levels),                                                     &
                                 ! IN A buoyancy parameter for cloudy
                                 !    air on p,T,q-levels
                                 !    (full levels).
 bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                 ! IN A grid-box mean buoyancy param
                                 ! on p,T,q-levels (full levels).
 bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                 ! IN A grid-box mean buoyancy param
                                 ! on p,T,q-levels (full levels).
 a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                                 ! IN Saturated lapse rate factor
                                 !    on p,T,q-levels (full levels).
 a_dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                                 ! IN Saturated lapse rate factor
                                 !    on p,T,q-levels (full levels).
 dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                                 ! IN Derivative of q_SAT w.r.t. T

REAL, INTENT(IN) ::                                                     &
 recip_l_mo_sea(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),   &
                                 ! IN Reciprocal of the surface
                                 !    Obukhov length over sea (m^-1).
 flandg(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end)
                                 ! IN Land fraction on all tiles

! (f) Atmospheric + any other data not covered so far, incl control.
REAL, INTENT(IN) ::                                                     &
 rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        2,bl_levels),                                                   &
                                  ! IN (LW,SW) rad heating rate (K/s)
  micro_tends(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,      &
              2, bl_levels),                                            &
                         ! Tendencies from microphys within BL levels
                         ! (TL, K/s; QW, kg/kg/s)
 fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                  ! IN Surface flux buoyancy over
                                  ! density (m^2/s^3)

 u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                                  ! IN Surface friction velocity
                                  !    (m/s)
 pstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                  ! IN Surface pressure (Pascals).
 tstar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                                  ! IN Surface temperature (K).
 h_blend_orog(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                                  ! IN Blending height used as part
                                  ! of effective roughness scheme
 zh_prev(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                  ! IN boundary layer height from
                                  !    previous timestep
 rib_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),           &
                               ! IN  Bulk Richardson number for lowest
                               ! layer
 sil_orog_land(land_pts),                                               &
                               ! IN Silhouette area of unresolved
                               ! orography per unit horizontal area
 zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                               ! IN Height of top of initial
                               !     parcel ascent
 z_lcl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! IN Height of LCL

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
 nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
 L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

REAL, INTENT(IN) ::                                                     &
 u_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end),   &
                                 ! IN W'ly component of surface
!                                       current (m/s). P grid
   v_0_px(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end), &
                                   ! IN S'ly component of surface
!                                       current (m/s). P grid
   ho2r2_orog(land_pts),                                                &
                                   ! IN peak to trough height of
!                                       unresolved orography
!                                       on land points only (m)
   sd_orog(land_pts),                                                   &
                                   ! IN Standard Deviation of unresolved
!                                       orography on land points only (m)
   z0m_eff_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                ! IN Effective grid-box roughness
!                                 length for momentum

INTEGER, INTENT(IN) ::                                                  &
 land_index(land_pts)        ! IN LAND_INDEX(I)=J => the Jth
!                                     point in P_FIELD is the Ith
!                                     land point.
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
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),      &
                                   ! IN temperature
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),    &
                                 ! IN Total water content
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)
                                 ! IN Ice/liquid water temperature

! INOUT variables
REAL, INTENT(INOUT) ::                                                  &
 zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                                 ! INOUT Height above surface of top
                                 !       of boundary layer (metres).
 dzh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                                 ! INOUT inversion thickness (m)
 fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                                 ! INOUT Moisture flux between layers
!                                     (kg per square metre per sec).
!                                     FQW(,1) is total water flux
!                                     from surface, 'E'.
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                                   ! INOUT FTL(,K) contains net turbulent
!                                     sensible heat flux into layer K
!                                     from below; so FTL(,1) is the
!                                     surface sensible heat, H. (W/m2)
   rhokh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                   ! INOUT Exchange coeffs for moisture.
    w(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,0:bl_levels), &
    etadot(wdims%i_start:wdims%i_end,wdims%j_start:wdims%j_end,         &
           0:bl_levels),                                                &
   t1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),          &
                                    ! INOUT Standard deviation of
                                    ! turbulent fluctuations of layer 1
                                    ! temperature; for use in
                                    ! initiating convection.
   q1_sd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                    ! INOUT Standard deviation of turbulent
                                    !    fluctuations of layer 1
                                    !    humidity; for use in initiating
                                    !    convection.

REAL, INTENT(INOUT) ::                                                  &
 rhokm(pdims_s%i_start:pdims_s%i_end,                                   &
       pdims_s%j_start:pdims_s%j_end ,bl_levels)
!            Exchange coefficients for momentum on P-grid

LOGICAL, INTENT(INOUT) ::                                               &
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! INOUT Logical switch for trade Cu
 l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! INOUT Flag to indicate shallow
                                 !     convection

INTEGER, INTENT(INOUT) ::                                               &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                               ! INOUT Number of model layers in the
                               !    turbulently mixed layer
 ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! INOUT Top level of initial parcel
                               !  ascent. Used in convection scheme.

!  Outputs :-
!  (a) Calculated anyway (use STASH space from higher level) :-
REAL, INTENT(OUT) ::                                                    &
 f_ngstress(pdims_s%i_start:pdims_s%i_end,                              &
            pdims_s%j_start:pdims_s%j_end,2:bl_levels),                 &
 tau_fd_x(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,  &
          bl_levels),                                                   &
 tau_fd_y(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,  &
          bl_levels),                                                   &
  bl_type_1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                               ! OUT Indicator set to 1.0 if stable
                                 !     b.l. diagnosed, 0.0 otherwise.
  bl_type_2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                               ! OUT Indicator set to 1.0 if Sc over
                                 !     stable surface layer diagnosed,
                                 !     0.0 otherwise.
  bl_type_3(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                               ! OUT Indicator set to 1.0 if well
                                 !     mixed b.l. diagnosed,
                                 !     0.0 otherwise.
  bl_type_4(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                               ! OUT Indicator set to 1.0 if
                                 !     decoupled Sc layer (not over
                                 !     cumulus) diagnosed,
                                 !     0.0 otherwise.
  bl_type_5(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                               ! OUT Indicator set to 1.0 if
                                 !     decoupled Sc layer over cumulus
                                 !     diagnosed, 0.0 otherwise.
  bl_type_6(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                               ! OUT Indicator set to 1.0 if a
                                 !     cumulus capped b.l. diagnosed,
                                 !     0.0 otherwise.
  bl_type_7(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! OUT Indicator set to 1.0 if a
                                 !     Shear-dominated unstable b.l.
                                 !     diagnosed, 0.0 otherwise.

REAL, INTENT(OUT) :: bl_w_var(   tdims%i_start : tdims%i_end,           &
                                 tdims%j_start : tdims%j_end,           &
                                             2 : tdims%k_end+1 )

REAL, INTENT(OUT) ::                                                    &
  zht(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                                 ! OUT Max height of turb mixing
  zhnl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                 ! OUT non-local PBL depth
  wstar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! OUT Convective velocity scale (m/s)
  wthvs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! OUT surface flux of thv (Km/s)
  shallowc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                 ! OUT Shallow Cu diagnostic
                                 !   Indicator set to 1.0 if shallow,
                                 !   0.0 if not shallow or not cumulus
  cu_over_orog(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                 ! OUT Indicator for cumulus
                                 !     over steep orography
                                 !   Indicator set to 1.0 if true,
                                 !   0.0 if false. Exclusive.
  we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                  ! OUT rho*entrainment rate implied b
                                  !     placing of subsidence
  zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),          &
                                  ! OUT (z-z_base)/(z_i-z_base)
  t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),        &
                                  ! OUT a fraction of the timestep
  we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                  ! OUT rho*entrainment rate implied b
                                  !     placing of subsidence
  zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),      &
                                  ! OUT (z-z_base)/(z_i-z_base)
  t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),    &
                                  ! OUT a fraction of the timestep
  zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! OUT Top of decoupled layer

INTEGER, INTENT(OUT) ::                                                 &
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                 ! OUT Top level for turb mixing in
                                 !     any decoupled Sc layer
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                 ! OUT Bottom level of any decoupled
                                 !     turbulently-mixed Sc layer.
 kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                                 ! OUT grid-level of SML inversion
 kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! OUT grid-level of DSC inversion

!-2 Genuinely output, needed by other atmospheric routines :-
REAL, INTENT(OUT) ::                                                    &
  uw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                           ! OUT U-component of surface wind stress
                           !     on P-grid
  vw0(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                           ! OUT V-component of surface wind stress
                           !     on P-grid
REAL, INTENT(OUT) :: rhcpt(tdims%i_start:tdims%i_end,                   &
                           tdims%j_start:tdims%j_end,1:tdims%k_end)
!-----------------------------------------------------------------------
!   Symbolic constants (parameters) reqd in top-level routine :-

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'BDY_EXPL2'
REAL :: TmpScm2d(ScmRowLen,ScmRow)      ! Temporary for SCM output
REAL :: sl(ScmRowLen,ScmRow,bl_levels)  ! Static energy

! Parameters also passed to EX_COEF
! Layer interface K_LOG_LAYR-1/2 is the highest which requires log
! profile correction factors to the vertical finite differences.
! The value should be reassessed if the vertical resolution is changed.
! We could set K_LOG_LAYR = BL_LEVELS and thus apply the correction
! factors for all the interfaces treated by the boundary layer scheme;
! this would be desirable theoretically but expensive computationally
! because of the use of the log function.
INTEGER ::    k_log_layr
PARAMETER (k_log_layr=2)
!-----------------------------------------------------------------------
!  Workspace :-
REAL ::                                                                 &
 recip_time_sbl, recip_time_cbl ! inverse timescales for TKE diagnostic
REAL ::                                                                 &
 a_dqsdtm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                              ! Saturated lapse rate factor
                              ! on intermediate levels (half levels).
 a_qsm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                              ! Saturated lapse rate factor
                              ! on intermediate levels (half levels).
 bqm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                              ! A buoyancy parameter for clear air
                              ! on intermediate levels (half levels).
 bqm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                              ! A buoyancy parameter for cloudy air
                              ! on intermediate levels (half levels).
 btm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                              ! A buoyancy parameter for clear air
                              ! on intermediate levels (half levels).
 btm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                              ! A buoyancy parameter for cloudy air
                              ! on intermediate levels (half levels).
 cfm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                              ! Estimate of cloud fraction
                              ! on intermediate levels (half levels).
 dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
      2:bl_levels),                                                     &
                              ! Buoyancy gradient across layer
                              !  interface.
 dbdz_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
                              ! Buoyancy gradient across layer
                              !  interface, inc gradient adjustment
 dvdzm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
       2:bl_levels),                                                    &
                              ! Modulus of wind shear.
 rmlmax2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! Square of asymptotic mixing length 
                              ! for Smagorinsky scheme
 ri(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),   &
                              ! Local Richardson number.
 ri_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
       2:bl_levels),                                                    &
                              ! Local Richardson number, inc grad adj
 grad_q_adj(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),       &
                              ! Humidity gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
   grad_t_adj(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),     &
                                ! Temperature gradient adjustment
!                                 for non-local mixing in unstable
!                                 turbulent boundary layer.
   rhokhz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,          &
          2:bl_levels),                                                 &
                                ! Non-local turbulent mixing
!                                 coefficient for heat and moisture.
   rhokh_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             2:bl_levels),                                              &
                                ! Non-local turbulent mixing coefficient
                                ! for top-down mixing of heat and
                                ! moisture.
   rhokh_th(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                                ! local scheme rhokh on th-levels,
                                ! index k held on th-level(k-1),
                                ! same as rhokm
   rhokmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          2:bl_levels),                                                 &
                                ! Non-local turbulent mixing
!                                 coefficient for momentum.
   rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,       &
             2:bl_levels),                                              &
                                ! Non-local turbulent mixing coefficient
                                ! for top-down mixing of momentum.
   weight_1dbl(pdims%i_start:pdims%i_end,                               &
               pdims%j_start:pdims%j_end ,bl_levels),                   &
                                ! Weighting applied to 1D BL scheme
                                ! to blend with Smagorinsky scheme,
                                ! index k held on theta level (k-1)
   weight_1dbl_rho(pdims%i_start:pdims%i_end,                           &
                   pdims%j_start:pdims%j_end,bl_levels),                &
                                ! weight_1dbl interpolated to rho levels
   elm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),&
                                ! Mixing length for momentum
   elh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,2:bl_levels),&
   elh_rho(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
           2:bl_levels),                                                &
                                ! Mixing length for heat (m),
                                ! held on theta and rho levels, resp.
   fm_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! stability function for momentum transport
                                ! level 1 value is dummy
   fh_3d(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),&
                                ! stability function for heat and moisture.
                                ! level 1 value is dummy
   sigma_h(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! Standard deviation of subgrid
                                ! orography for sg mixing options (m)

REAL, ALLOCATABLE :: visc_h_rho (:,:,:)                 ! visc_h on rho levels

    ! Terms for non-gradient flux parametrization
    !  (=0 unless using 9C code with FLUX_GRAD=LockWhelan2006)
REAL ::                                                                 &
  ft_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels+1),                                                   &
                              ! Non-turbulent heat and moisture flux
  fq_nt(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels+1)          !  (on rho levels, surface flux(K=1)=0)
REAL ::                                                                 &
  rhof2(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        2:bl_levels),                                                   &
                              ! f2 and fsc term shape profiles
  rhofsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         2:bl_levels)

REAL ::                                                                 &
  tothf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! Total heat fluxes at inversions
  tothf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &

  totqf_zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                              ! Total moisture fluxes at inversions
  totqf_zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &

  ft_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                              ! Non-turbulent heat and moisture flux
  fq_nt_dscb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                              !    at the base of the DSC layer.

REAL ::                                                                 &
zh_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                              ! Height above surface of top of
                              !  boundary layer (metres) as
                              !  determined from the local
                              !  Richardson number profile.
riout(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                              ! Gradient Richardson number
                              ! for SCM output
                              ! RIOUT(K) is on theta-level K
dsldz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
      bl_levels),                                                       &
                              ! TL+gz/cp gradient between
                              ! levels K and K-1
dsldz_ga(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                              ! TL+gz/cp gradient between
                              ! levels K and K-1, inc gradient adjust
dqwdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels)
                              ! QW gradient between

INTEGER ::                                                              &
 ntml_local(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! Number of model layers in the
                                 ! turbulently mixed layer as
                                 ! determined from the local
                                 ! Richardson number profile.
 ntml_nl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                 ! Number of model layers in the
                                 ! turbulently mixed layer as
                                 ! determined from the parcel ascent.
 ntml_save(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                 ! saved copy of ntml on entry
 sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                 ! Flags for whether discontinuous
 dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                 ! inversions are diagnosed

LOGICAL ::                                                              &
 unstable(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                               ! Logical switch for unstable
                               !    surface layer.
 dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                               ! Flag set if decoupled
                               ! stratocumulus layer found
 coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                               ! Flag to indicate Sc layer weakly
                               ! coupled to surface (ie weakly
                               ! decoupled)
 dynamic_bl_diag(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                               ! Flag to indicate the dynamic
                               ! diagnosis (iDynDiag) has
                               ! determined the BL type
 topbl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                               ! Flag for having reached
                               ! the top of the turbulently mixed
                               ! layer.

!  Local scalars :-
REAL ::                                                                 &
 dzu,                                                                   &
            ! Westerly wind shear between levels K+1 and K.
 dzv,                                                                   &
            ! Southerly wind shear between levels K+1 and K.
 lambda_min,                                                            &
            ! Min value of length scale LAMBDA.
 lambdah,                                                               &
            ! Asymptotic mixing length for turbulent transport
            ! of heat/moisture.
 vkz,                                                                   &
            ! Temporary in calculation of ELH.
 f_log,                                                                 &
            ! Temporary in calculation of logarithmic correction
 zmaxb_for_dsc,                                                         &
 zmaxt_for_dsc
            ! Max heights to look for DSC cloud base and top

REAL ::                                                                 &
  weight1,                                                              &
  weight2,                                                              &
  weight3,                                                              &
  z_scale,                                                              &
             ! scaling with height
  zpr,                                                                  &
             ! z/sigma_h
  slope,                                                                &
             ! subgrid orographic slope
  dsldzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
             ! TL+gz/cp gradient interpolated to Z_TQ
  dsldzm_ga,                                                            &
             ! dsldzm with gradient adjustment terms
  dqwdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         2:bl_levels),                                                  &
             ! QW gradient interpolated to Z_TQ
  qssurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
             ! qsat of surface

! variables for rhcrit parametrization
REAL :: b2, sh, sgm, qsw, exner, root6, delta_x,                        &
       qsw_arr(tdims%i_start:tdims%i_end),                              &
       max_rhcpt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),  &
       min_rhcpt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)

INTEGER ::                                                              &
 i,j,iScm,jScm,                                                         &
                ! LOCAL Loop counter (horizontal field index).
 k,ient,                                                                &
                ! LOCAL Loop counter (vertical level index).
 kmax(1),                                                               &
                ! level of max rhokm,
 l,                                                                     &
                ! LOCAL Loop counter for land points
 ntop       ! NTPAR restricted below BL_LEVELS

INTEGER ::                                                              &
 omp_block,                                                             &
                 ! for open mp blocking
 jj
                 ! for indexing over open mp block

REAL, PARAMETER :: max_abs_obkhov = 1.0e6
                 ! Maximum permitted magnitude of the Obukhov
                 ! length (m).

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1) Set various diagnostics and switches
!-----------------------------------------------------------------------
! checking of nl_bl_levels has been moved to readsize/scm_shell
! so that it is only executed at initialisation

!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( pdims, dynamic_bl_diag, ntml_save, ntml, bl_levels,     &
!$OMP&         weight_1dbl, weight_1dbl_rho ) PRIVATE( i, j, k )
!$OMP  DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    dynamic_bl_diag(i,j) = .FALSE.
    ntml_save(i,j) = ntml(i,j)
  END DO
END DO
!$OMP END DO NOWAIT

!------------------------------------------------------------------
!  Initialize weighting applied to 1d BL scheme
!  (used to blend with 3D Smagorinsky scheme)
!------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO k = 1, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      weight_1dbl(i,j,k) = 1.0
      weight_1dbl_rho(i,j,k) = 1.0
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

!-----------------------------------------------------------------------
! Set surface scaling diagnostics
!-----------------------------------------------------------------------
 ! Obukhov length
IF (BL_diag%l_oblen) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,BL_diag,u_s,fb_surf)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      !       Limit the magnitude of the Obukhov length to avoid
      !       problems with packing.
      BL_diag%oblen(i,j)= u_s(i,j)*u_s(i,j)*u_s(i,j)
      IF ( BL_diag%oblen(i,j) <                                         &
           max_abs_obkhov*ABS(vkman*fb_surf(i,j)) ) THEN
        BL_diag%oblen(i,j)=-BL_diag%oblen(i,j)/(vkman*fb_surf(i,j))
      ELSE
        BL_diag%oblen(i,j)=-SIGN(max_abs_obkhov, fb_surf(i,j))
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Ustar
IF (BL_diag%l_ustar) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,BL_diag,u_s)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      BL_diag%ustar(i,j)=u_s(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Surface buoyancy flux
IF (BL_diag%l_wbsurf) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,BL_diag,fb_surf)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      BL_diag%wbsurf(i,j)=fb_surf(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO
END IF
!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_bl) .AND.                                       &
     (model_type == mt_single_column) ) THEN

  DO k=1, bl_levels
    DO j=pdims%j_start, pdims%j_end
      jScm =  j - pdims%j_start + 1
      DO i=pdims%i_start, pdims%i_end
        iScm = i - pdims%i_start + 1
        sl(iScm,jScm,k) = tl(i,j,k) + grcp*z_tq(i,j,k)
      END DO ! I
    END DO ! j
  END DO ! K

  !       Output SL

  CALL scmoutput(sl,'SL',                                               &
       'Liquid/frozen water static energy (IN)','K',                    &
       t_avg,d_bl,default_streams,'',routinename)

  !       Output QW

  CALL scmoutput(qw,'qw',                                               &
       'Total water content (IN)','kg/kg',                              &
       t_avg,d_bl,default_streams,'',routinename)

END IF ! scmdiag_bl / model_type

!-----------------------------------------------------------------------
! 2.  Interpolate BT and BQ to half levels and calculate Ri
!-----------------------------------------------------------------------
CALL btq_int (                                                          &
! IN levels
   bl_levels,                                                           &
! IN fields
   z_tq,z_uv,bq,bt,bq_cld,bt_cld,a_qs,a_dqsdt,cf_bulk,                  &
! OUT fields
   bqm,btm,bqm_cld,btm_cld,a_qsm,a_dqsdtm,cfm                           &
    )
!-----------------------------------------------------------------------
! Calculate lapse rates
!-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, weight1,  weight2,     &
!$OMP& weight3, zpr, dzv, dzu, l, slope, dsldzm_ga)

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    grad_t_adj(i,j) = MIN( max_t_grad,                                  &
                           a_grad_adj * t1_sd(i,j) / zh_prev(i,j) )
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      dsldz(i,j,k)    = ( tl(i,j,k) - tl(i,j,k-1) )                     &
                              *rdz_charney_grid(i,j,k) + grcp
      dsldz_ga(i,j,k) = dsldz(i,j,k)
      IF ( z_tq(i,j,k) <= zh_prev(i,j) ) THEN
        dsldz_ga(i,j,k) = dsldz_ga(i,j,k) - grad_t_adj(i,j)
      END IF
      dqwdz(i,j,k)    = ( qw(i,j,k) - qw(i,j,k-1) )                     &
                             * rdz_charney_grid(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO

! We need a grid-box mean surface-to-level-1 buoyancy gradient in order
! to construct Ri on level 1.  Either the gradient from level 1 to 2
! can be extrapolated, or surface properties can be used.  Over a
! heterogeneous land surface this is poorly defined and we can't use Rib
! from the surface scheme as vertically averaging Ri is numerically
! unstable.  So, over land, only the average temperature gradient is used

!$OMP MASTER
IF ( l_new_qsat_bl ) THEN
  IF ( l_mr_physics ) THEN
    CALL qsat_mix_new(qssurf,tstar,pstar,pdims%i_len,pdims%j_len) ! No halos
  ELSE
    CALL qsat_new(qssurf,tstar,pstar,pdims%i_len,pdims%j_len) ! No halos
  END IF
ELSE
  ! DEPENDS ON: qsat_mix
  CALL qsat_mix(qssurf,tstar,pstar,pdims%i_end*pdims%j_end,l_mr_physics)
END IF
!$OMP END MASTER
!$OMP BARRIER

k=1
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    dsldz(i,j,k)    = ( tl(i,j,k) - tstar(i,j) )                        &
                              *rdz_charney_grid(i,j,k) + grcp
    dsldz_ga(i,j,k) = dsldz(i,j,k) ! no GA below level 1
    IF ( flandg(i,j) < 0.2) THEN
      dqwdz(i,j,k)  = ( qw(i,j,k) - qssurf(i,j) )                       &
                           * rdz_charney_grid(i,j,k)
    ELSE
      dqwdz(i,j,k)  = dqwdz(i,j,2) ! extrapolate qw if mainly land
    END IF
  END DO
END DO
!$OMP END DO

!       Output gradient adjusted SL
!$OMP MASTER
IF ( l_scmdiags(scmdiag_bl) .AND. idyndiag == DynDiag_Ribased .AND.     &
     model_type == mt_single_column ) THEN

  DO k=1, bl_levels
    DO j=pdims%j_start, pdims%j_end
      jScm =  j - pdims%j_start + 1
      DO i=pdims%i_start, pdims%i_end
        iScm = i - pdims%i_start + 1

        sl(iScm,jScm,k) = tl(i,j,k) + grcp*z_tq(i,j,k)
        IF ( z_tq(i,j,k) <= zh_prev(i,j) ) THEN
          sl(iScm,jScm,k) = sl(iScm,jScm,k) - grad_t_adj(i,j)*z_tq(i,j,k)
        END IF
      END DO ! i
    END DO ! j
  END DO ! k

  CALL scmoutput(sl,'SL_GA',                                            &
       'Gradient-adjusted SL','K',                                      &
       t_avg,d_bl,default_streams,'',routinename)
END IF ! scmdiag_bl / model_type / idyndiag

!$OMP END MASTER
!$OMP BARRIER

! Local Ri-based calculation of RHOKM and RHOKH:
! Calculate `buoyancy' gradient, DBDZ, on theta-levels
! NOTE: DBDZ(K) is on theta-level K-1

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      weight1 = r_rho_levels(i,j,k) -                                   &
                r_rho_levels(i,j,k-1)
      weight2 = r_theta_levels(i,j,k-1)-                                &
                r_rho_levels(i,j,k-1)
      weight3 = r_rho_levels(i,j,k) -                                   &
                r_theta_levels(i,j,k-1)
      dsldzm(i,j,k) = weight2 * dsldz(i,j,k)                            &
             + weight3 * dsldz(i,j,k-1)
      dqwdzm(i,j,k) = weight2 * dqwdz(i,j,k)                            &
             + weight3 * dqwdz(i,j,k-1)
      dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dsldzm(i,j,k) +                  &
                        bq_gb(i,j,k-1)*dqwdzm(i,j,k) )/weight1
      !       ! Now with gradient adjustment
      dsldzm_ga = weight2 * dsldz_ga(i,j,k)                             &
             + weight3 * dsldz_ga(i,j,k-1)
      dbdz_ga(i,j,k) = g*( bt_gb(i,j,k-1)*dsldzm_ga +                   &
                           bq_gb(i,j,k-1)*dqwdzm(i,j,k) )/weight1
    END DO
  END DO
END DO
!$OMP END DO

! Overwrite element 2 with level 1 to 2 gradients
! (ie instead of using centred value that uses surface parameters)
IF ( .NOT. l_use_surf_in_ri ) THEN
  k = 2
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      dbdz(i,j,k) = g*( bt_gb(i,j,k-1)*dsldz(i,j,k) +                   &
                        bq_gb(i,j,k-1)*dqwdz(i,j,k) )
      dbdz_ga(i,j,k) = g*( bt_gb(i,j,k-1)*dsldz_ga(i,j,k) +             &
                           bq_gb(i,j,k-1)*dqwdz(i,j,k) )
    END DO
  END DO
!$OMP END DO
END IF
!--------------------------------------------------
! Calculate modulus of shear on theta-levels
! dvdzm(k) is on theta-level(k-1)
!--------------------------------------------------
IF (.NOT. l_subfilter_vert) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dzu = u_p(i,j,k) - u_p(i,j,k-1)
        dzv = v_p(i,j,k) - v_p(i,j,k-1)
        dvdzm(i,j,k) = MAX ( 1.0e-12 ,                                  &
                 SQRT(dzu*dzu + dzv*dzv) * rdz(i,j,k)  )
      END DO
    END DO
  END DO
!$OMP END DO

ELSE

  ! On entry, visc_m is 3D shear(k) on theta-level(k)

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        dvdzm(i,j,k) = MAX( 1.0e-12 , visc_m(i,j,k-1) )
      END DO
    END DO
  END DO
!$OMP END DO

END IF

IF (l_subfilter_horiz .OR. l_subfilter_vert) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      rmlmax2(i,j) = ( mix_factor * delta_smag(i,j) )**2
   END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rneutml_sq(i,j,k) = 1.0 / (                                    &
                 1.0/( vkman*(z_tq(i,j,k) + z0m_eff_gb(i,j)) )**2      &
               + 1.0/rmlmax2(i,j) )
      END DO
    END DO
  END DO
!$OMP END DO

END IF
!-----------------------------------------------------------------------
! 2.1 Orographic enhancement of subgrid mixing
!-----------------------------------------------------------------------
!  Set-up 2D array for standard deviation of subgrid orography.
!-----------------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    sigma_h(i,j) = 0.0
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO l = 1, land_pts
  j=(land_index(l)-1)/pdims%i_end + 1
  i=land_index(l) - (j-1)*pdims%i_end
  sigma_h(i,j) = MIN( sd_orog(l), 300.0 )
END DO
!$OMP END DO
!-----------------------------------------------------------------------
!  Enhance resolved shear through unresolved subgrid drainage flows.
!-----------------------------------------------------------------------
IF (sg_orog_mixing == sg_shear .OR.                                     &
    sg_orog_mixing == sg_shear_enh_lambda) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        IF (sigma_h(i,j) > 1.0 ) THEN
          zpr = z_tq(i,j,k-1)/sigma_h(i,j)
          ! Height dependence, to reduce effect to zero with height
          !   gives z_scale~[1,0.95,0.5,0] at zpr=[0,0.6,1,1.7]
          weight1 = 0.5*( 1.0 - TANH(4.0*(zpr-1.0) ) )

          ! Take slope ~ sd/h_scale for small sd;
          !            tends to 0.2 for large sd
          slope = 1.0 / SQRT( 25.0 + (h_scale/sigma_h(i,j))**2 )

          dvdzm(i,j,k) = MAX ( dvdzm(i,j,k),                            &
                               weight1*slope*t_drain*dbdz(i,j,k) )

          IF (k==2 .AND. BL_diag%l_dvdzm)                               &
            BL_diag%dvdzm(i,j,1)=weight1*slope*t_drain*dbdz(i,j,k)

        END IF
      END DO
    END DO
  END DO
!$OMP END DO
END IF     ! sg_orog_mixing

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      ri(i,j,k)    = dbdz(i,j,k)    / ( dvdzm(i,j,k)*dvdzm(i,j,k) )
      ri_ga(i,j,k) = dbdz_ga(i,j,k) / ( dvdzm(i,j,k)*dvdzm(i,j,k) )
    END DO
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

IF (BL_diag%l_gradrich) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,BL_diag,pdims,ri)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%gradrich(i,j,k)=ri(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF (BL_diag%l_dbdz) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,pdims,BL_diag,dbdz)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%dbdz(i,j,k)=dbdz(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF (BL_diag%l_dvdzm) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,pdims,BL_diag,dvdzm)
  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%dvdzm(i,j,k)=dvdzm(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF
!-----------------------------------------------------------------------
! 3.  Orographic formdrag - distributed drag option
!-----------------------------------------------------------------------
IF (formdrag ==  explicit_stress) THEN
  !------------------------------------------------------------------
  !      Set stresses to zero
  !------------------------------------------------------------------
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,tdims,tau_fd_x,tau_fd_y)
  DO k = 1, bl_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        tau_fd_x(i,j,k) = 0.0
        tau_fd_y(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
  !------------------------------------------------------------------
  !      Calculate stress profiles
  !------------------------------------------------------------------
  ! DEPENDS ON: fm_drag
  CALL fm_drag (                                                        &
  ! IN levels
        land_pts, land_index, bl_levels,                                &
  ! IN fields
        u_p, v_p, rho_wet_tq, z_uv, z_tq, z0m_eff_gb, zh_prev, rib_gb,  &
        sil_orog_land,                                                  &
  ! OUT fields
        tau_fd_x(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 1:bl_levels),                                          &
        tau_fd_y(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 1:bl_levels)                                           &
        )
  !------------------------------------------------------------------
  !      Orographic stress diagnostics
  !------------------------------------------------------------------
  IF (BL_diag%l_ostressx) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,tdims,BL_diag,tau_fd_x)
    DO k = 1, bl_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%ostressx(i,j,k)=tau_fd_x(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  IF (BL_diag%l_ostressy) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,tdims,BL_diag,tau_fd_y)
    DO k = 1, bl_levels
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          BL_diag%ostressy(i,j,k)=tau_fd_y(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

END IF
!-----------------------------------------------------------------------
! 4. Apply dynamic diagnosis of shear-driven layers.
!-----------------------------------------------------------------------
!       In cases where the parcel ascent continues right through the
!       boundary layer, we diagnose cumulus only where the surface
!       buoyancy flux is sufficiently unstable, testing the ratio of
!       the depth of the inversion to the Obukhov length. A value of
!       1 -- 2 is reasonable for this test and 1.6 is selected, but
!       no great precision is attached to this value. Since this is
!       of importance mainly at sea-points, to avoid complications
!       with coastal tiling, the scheme operates only at points
!       where the land fraction is below 0.5.

omp_block = pdims%j_end
!$ omp_block = CEILING(REAL(pdims%j_end)/omp_get_max_threads())

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, jj, ntop, z_scale)

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ! First override the provisional cumulus diagnosis if the
    ! actual surface buoyancy flux is stable
    IF ( fb_surf(i,j) < 0.0 ) THEN
      cumulus(i,j)   = .FALSE.
      l_shallow(i,j) = .FALSE.
      ntml(i,j)      = 1
      zh(i,j)        = z_uv(i,j,2)
    END IF
  END DO
END DO
!$OMP END DO


IF (idyndiag == DynDiag_ZL) THEN

  ! Original version - causes spuriously deep boundary layers if
  ! BL_LEVELS is >> 3km

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      IF ( flandg(i,j) < 0.5 ) THEN
        ntop = MIN(ntpar(i,j),bl_levels-1)
        IF ( -z_uv(i,j,ntop+1) * recip_l_mo_sea(i,j) < 1.6 ) THEN
          cumulus(i,j)   = .FALSE.
          l_shallow(i,j) = .FALSE.
          ntml(i,j)      = ntop
          zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
        END IF
      END IF

    END DO
  END DO
!$OMP END DO

ELSE IF (idyndiag == DynDiag_ZL_corrn) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

        !------------------------------------------------------------
        ! As original code above, except restrict depth scale to <3km
        ! and, if near-neutral, ignore the BL depth diagnosed by the
        ! adiabatic parcel (ie ZH, NTML) completely.
        !------------------------------------------------------------
      IF ( flandg(i,j) < 0.5 ) THEN
        ntop = MIN(ntpar(i,j),bl_levels-1)
        z_scale = MIN( 3000.0, z_uv(i,j,ntop+1) )
        IF ( -z_scale*recip_l_mo_sea(i,j) < 1.6 ) THEN
          cumulus(i,j)   = .FALSE.
          l_shallow(i,j) = .FALSE.
          ntml(i,j)      = 1
          zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
        END IF
      END IF

    END DO
  END DO
!$OMP END DO

ELSE IF (idyndiag == DynDiag_ZL_CuOnly) THEN

!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

        !------------------------------------------------------------
        ! As DynDiag_ZL_corrn but only affects cumulus points and
        ! points that are entirely sea.  Note that DynDiag_ZL_corrn
        ! has been found to switch off non-local mixing in
        ! stratocumulus, where surface fluxes are typically small
        !------------------------------------------------------------
      IF ( cumulus(i,j) .AND. flandg(i,j) < 0.01 ) THEN
        ntop = MIN(ntpar(i,j),bl_levels-1)
        z_scale = MIN( 3000.0, z_uv(i,j,ntop+1) )
        IF ( -z_scale*recip_l_mo_sea(i,j) < 1.6 ) THEN
          ! - ZH/L indicates BL close to neutral
          dynamic_bl_diag(i,j) = .TRUE.
          cumulus(i,j)   = .FALSE.
          l_shallow(i,j) = .FALSE.
          IF (l_fix_dyndiag) THEN
            ntml(i,j) = 1
          ELSE
            ! note that ntop can be very high making a spuriously
            ! deep well-mixed PBL - hence depricated
            ntml(i,j) = ntop
          END IF
          zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
        END IF
      END IF

    END DO
  END DO
!$OMP END DO

ELSE IF (idyndiag == DynDiag_Ribased ) THEN
  !------------------------------------------------------------
  ! As DynDiag_ZL_CuOnly but also allow ZH(Ri) to overrule the
  ! Cumulus diagnosis
  !------------------------------------------------------------
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      topbl(i,j)           = .FALSE.
    END DO
  END DO
!$OMP END DO
  !---------------------------------------------------------------
  !  Loop over levels to find Ri > RiCrit_sharp (=0.25) to find
  !  level to which Ri really is close to neutral
  !---------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
  DO jj = pdims%j_start, pdims%j_end, omp_block
    DO k = 2, bl_levels
      DO j = jj, MIN(jj+omp_block-1,pdims%j_end)
        DO i = pdims%i_start, pdims%i_end
          IF ( .NOT. topbl(i,j) .AND.                                   &
            (ri_ga(i,j,k) >  RiCrit_sharp .OR. k > bl_levels-1) ) THEN
            topbl(i,j) = .TRUE.
            zh_local(i,j) = z_uv(i,j,k)
          END IF
        END DO  ! Loop over points
      END DO  ! Loop over points
    END DO  ! Loop over levels
  END DO
!$OMP END DO
  !---------------------------------------------------------------
  !  Overrule Cumulus flag where close to neutral BL
  !---------------------------------------------------------------
  IF (l_fix_dyndiag) THEN

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        IF ( cumulus(i,j) .AND. flandg(i,j) < 0.01 ) THEN
          ntop = MIN(ntpar(i,j),bl_levels-1)
          z_scale = MIN( 3000.0, z_uv(i,j,ntop+1) )
          IF ( zh_local(i,j)                                            &
                  > zh(i,j)+zhloc_depth_fac*(zhpar(i,j)-zh(i,j)) ) THEN
                ! ZH(Ri>RiCrit) more than zhloc_depth_fac up the
                ! cloud layer, indicating significant shear disruption
            dynamic_bl_diag(i,j) = .TRUE.
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = ntop
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          ELSE IF ( -z_scale*recip_l_mo_sea(i,j) < 1.6 ) THEN
            ! - ZH/L indicates BL close to neutral so overrule
            ! cumulus diagnosis and leave mixing to the local
            ! Ri-based scheme, given no better information on where
            ! unstable bl top should be
            dynamic_bl_diag(i,j) = .TRUE.
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = 1
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          END IF
        END IF ! Cu over sea

      END DO
    END DO
!$OMP END DO

  ELSE  ! old version: not good to set ntml=ntop for small zh/L

!$OMP DO SCHEDULE(STATIC)
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        IF ( cumulus(i,j) .AND. flandg(i,j) < 0.01 ) THEN
          ntop = MIN(ntpar(i,j),bl_levels-1)
          z_scale = MIN( 3000.0, z_uv(i,j,ntop+1) )
          IF ( -z_scale*recip_l_mo_sea(i,j) < 1.6 .OR.                  &
                ! - ZH/L indicates BL close to neutral
             zh_local(i,j) > zh(i,j)+zhloc_depth_fac*(z_scale-zh(i,j))  &
                ! ZH(Ri>RiCrit) more than zhloc_depth_fac up the
                ! cloud layer, indicating significant shear disruption
            ) THEN
            dynamic_bl_diag(i,j) = .TRUE.
            cumulus(i,j)   = .FALSE.
            l_shallow(i,j) = .FALSE.
            ntml(i,j)      = ntop
            zh(i,j)        = z_uv(i,j,ntml(i,j)+1)
          END IF
        END IF ! Cu over sea

      END DO
    END DO
!$OMP END DO

  END IF  ! l_fix_dyndiag



!$OMP SINGLE
  IF ( l_scmdiags(scmdiag_bl) .AND.                                     &
       model_type == mt_single_column) THEN
    CALL scmoutput(zh_local,'zh_loc_DD',                                &
         'ZH found from Ri under DynDiag','m',                          &
         t_avg,d_sl,default_streams,'',routinename)
  END IF
!$OMP END SINGLE

END IF  ! tests on iDynDiag

!$OMP END PARALLEL
!-----------------------------------------------------------------------
! 5.  Turbulent exchange coefficients and "explicit" fluxes between
!     model layers in the boundary layer (P243b, routine KMKH).
!-----------------------------------------------------------------------
! 5.1  Calculate the non-local terms and diffusion coefficients
!-----------------------------------------------------------------------
! Set NTML_NL to NTML as passed in from initial diagnosis routine
!-----------------------------------------------------------------------
zmaxb_for_dsc = 2500.0
zmaxt_for_dsc = 3000.0

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,ntml_nl,ntml)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ntml_nl(i,j) = ntml(i,j)
  END DO
END DO
!$OMP END PARALLEL DO

IF (nl_bl_levels < bl_levels) THEN
      ! Set to huge value to make if-test in KMKHZ redundent
  zmaxb_for_dsc = 1.0e10
  zmaxt_for_dsc = zmaxb_for_dsc

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,ntml_nl,nl_bl_levels,zh,z_uv)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( ntml_nl(i,j) > nl_bl_levels-1 ) THEN
        ntml_nl(i,j) = nl_bl_levels-1
        zh(i,j)      = z_uv(i,j,ntml_nl(i,j)+1)
      END IF
    END DO
  END DO
!$OMP END PARALLEL DO

END IF
!-----------------------------------------------------------------------
! Initialise non-local K and fluxes to zero; necessary for levels
! above NL_BL_LEVELS
!-----------------------------------------------------------------------
!$OMP  PARALLEL DEFAULT(NONE)                                          &
!$OMP& SHARED( bl_levels, pdims, ftl, fqw, rhokmz, rhokhz, rhokm_top,  &
!$OMP&         rhokh_top, f_ngstress, rhof2, rhofsc, ft_nt, fq_nt,     &
!$OMP&         tothf_zh, tothf_zhsc, totqf_zh, totqf_zhsc, ft_nt_dscb, &
!$OMP&         zhnl, zh, fq_nt_dscb ) PRIVATE( i, j, k )
!$OMP  DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      ftl(i,j,k) = 0.0
      fqw(i,j,k) = 0.0
      rhokmz(i,j,k) = 0.0
      rhokhz(i,j,k) = 0.0
      rhokm_top(i,j,k) = 0.0
      rhokh_top(i,j,k) = 0.0
      f_ngstress(i,j,k) = 0.0
      rhof2(i,j,k)  = 0.0
      rhofsc(i,j,k) = 0.0
      ! Initialise Lock-Whelan non-gradient terms to zero
      ! Only calculated in KMKHZ9C
      ft_nt(i,j,k)  = 0.0
      fq_nt(i,j,k)  = 0.0
    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ft_nt(i,j,1)  = 0.0
    fq_nt(i,j,1)  = 0.0
    tothf_zh(i,j)   = 0.0
    tothf_zhsc(i,j) = 0.0
    totqf_zh(i,j)   = 0.0
    totqf_zhsc(i,j) = 0.0
    ft_nt_dscb(i,j) = 0.0
    fq_nt_dscb(i,j) = 0.0
    zhnl(i,j) = zh(i,j)  ! initialise non-local PBL depth
                         ! to that diagnosed in conv_diag
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! switching off non_local_bl based on vertical Smagorinsky has been moved
! to readlsta/readlsta_4a/scm_shell so that it is only executed once

IF (non_local_bl == on) THEN

  IF (i_bl_vn == i_bl_vn_9b) THEN

    CALL kmkhz_9b (                                                     &
    !     IN levels/switches
             nl_bl_levels,BL_diag, nSCMDpkgs,L_SCMDiags,                &
    !     IN fields
             p_theta_levels,rho_wet_tq,rho_mix,rho_mix_tq,t,q,qcl,qcf,  &
             cf_bulk, qw,tl, dzl_charney,rdz_charney_grid,z_tq,z_uv,    &
             rad_hr,                                                    &
             bt,bq,btm,bqm,dqsdt,btm_cld,bqm_cld,a_qs,a_qsm,a_dqsdtm,   &
             u_s,fb_surf,rhostar,ntpar,zh_prev,                         &
             zhpar,zmaxb_for_dsc,zmaxt_for_dsc,l_shallow,               &
    !     INOUT fields
             ftl,fqw,zhnl,cumulus,ntml_nl,w,t1_sd,q1_sd,                &
    !     OUT fields
             rhokmz,rhokhz,rhokm_top,rhokh_top,zhsc,                    &
             unstable,dsc,coupled,sml_disc_inv,dsc_disc_inv,            &
             ntdsc,nbdsc,f_ngstress,                                    &
             grad_t_adj, grad_q_adj,                                    &
             kent, we_lim, t_frac, zrzi,                                &
             kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
             )

  ELSE IF (i_bl_vn == i_bl_vn_9c) THEN

    CALL kmkhz_9c (                                                     &
    !     IN levels/switches
             nl_bl_levels,BL_diag, nSCMDpkgs,L_SCMDiags,                &
    !     IN fields
             p_theta_levels,rho_wet_tq,rho_mix,rho_mix_tq,t,q,qcl,qcf,  &
             cf_bulk, qw,tl, dzl_charney,rdz_charney_grid,z_tq,z_uv,    &
             rad_hr,micro_tends,                                        &
             bt,bq,btm,bqm,dqsdt,btm_cld,bqm_cld,a_qs,a_qsm,a_dqsdtm,   &
             u_s,fb_surf,rhostar,ntpar,zh_prev,                         &
             zhpar,z_lcl,zmaxb_for_dsc,zmaxt_for_dsc,l_shallow,         &
    !     INOUT fields
             ftl,fqw,zhnl,dzh,cumulus,ntml_nl,w,etadot,t1_sd,q1_sd,     &
    !     OUT fields
             rhokmz,rhokhz,rhokm_top,rhokh_top,zhsc,                    &
             unstable,dsc,coupled,sml_disc_inv,dsc_disc_inv,            &
             ntdsc,nbdsc,f_ngstress,                                    &
             grad_t_adj, grad_q_adj,                                    &
             rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb,       &
             tothf_zh, tothf_zhsc, totqf_zh, totqf_zhsc,                &
             kent, we_lim, t_frac, zrzi,                                &
             kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                 &
             )

  END IF   ! test on i_bl_vn

ELSE   ! not NON_LOCAL_BL

       !-------------------------------------------------------------
       ! Set all variables from the non-local scheme to zero or "off"
       !  - reset all fluxes and K's arising from the non-local scheme
       !-------------------------------------------------------------

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i,ient)                                                       &
!$OMP SHARED(pdims,unstable,fb_surf,cumulus,l_shallow,sml_disc_inv,ntpar,     &
!$OMP        ntml_nl,zhnl,grad_t_adj,grad_q_adj,dsc,dsc_disc_inv,ntdsc,nbdsc, &
!$OMP        zhsc,coupled,kent,kent_dsc,t_frac,zrzi,we_lim,t_frac_dsc,        &
!$OMP        zrzi_dsc,we_lim_dsc)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
        ! surface mixed layer
      unstable(i,j) = (fb_surf(i,j) >  0.0)
      cumulus(i,j) = .FALSE.
      l_shallow(i,j) = .FALSE.
      sml_disc_inv(i,j) = 0
      ntpar(i,j)   = 0
      ntml_nl(i,j) = -1    ! to ensure correct diagnostics
      zhnl(i,j)      = 0.0
      grad_t_adj(i,j) = 0.0
      grad_q_adj(i,j) = 0.0
        ! decoupled mixed layer
      dsc(i,j)     = .FALSE.
      dsc_disc_inv(i,j) = 0
      ntdsc(i,j)   = 0
      nbdsc(i,j)   = 0
      zhsc(i,j)    = 0.0
      coupled(i,j) = .FALSE.
        ! entrainment variables for non-local tracer mixing
      kent(i,j) = 2
      kent_dsc(i,j) = 2
      DO ient = 1, 3
        t_frac(i,j,ient) = 0.0
        zrzi(i,j,ient)   = 0.0
        we_lim(i,j,ient) = 0.0
        t_frac_dsc(i,j,ient) = 0.0
        zrzi_dsc(i,j,ient)   = 0.0
        we_lim_dsc(i,j,ient) = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF  ! test on NON_LOCAL_BL

!-----------------------------------------------------------------------
! 5.1  Call local coeff calculation for levels 2 to bl_levels
!-----------------------------------------------------------------------
CALL ex_coef (                                                          &
! IN levels/logicals
   bl_levels,k_log_layr,nSCMDpkgs,L_SCMDiags,BL_diag,                   &
! IN fields
   sigma_h,flandg,dbdz,dvdzm,ri,rho_wet_tq,z_uv,z_tq,z0m_eff_gb,        &
   h_blend_orog,zhpar,ntpar,ntml_nl,ntdsc,nbdsc,u_p,v_p,u_s,fb_surf,    &
   qw,tl,                                                               &
! IN/OUT fields
   cumulus,weight_1dbl,                                                 &
! OUT fields
   lambda_min,zh_local,ntml_local,elm,elh,elh_rho,rhokm,rhokh_th,       &
   fm_3d,fh_3d                                                          &
   )

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_bl) .AND.                                       &
     model_type == mt_single_column ) THEN

  DO k=2, bl_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
            ! so RIOUT(K) is on theta-level K
        riout(i,j,k-1) = ri(i,j,k)
      END DO
    END DO
  END DO
  DO j=pdims%j_start, pdims%j_end
    DO i=pdims%i_start, pdims%i_end
      riout(i,j,bl_levels) = 0.0
    END DO
  END DO

  CALL scmoutput(riout,'Ri_bl',                                         &
       'Richardson number',' ',                                         &
       t_avg,d_bl,default_streams,'',routinename)

  DO k=2, bl_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
            ! so OUT(K) is on theta-level K
        riout(i,j,k-1) = dbdz(i,j,k)
      END DO
    END DO
  END DO
  CALL scmoutput(riout,'dbdz',                                          &
       'Buoyancy gradient in Ri',' ',                                   &
       t_avg,d_bl,default_streams,'',routinename)

  DO k=2, bl_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        ! so OUT(K) is on theta-level K
        riout(i,j,k-1) = dvdzm(i,j,k)
      END DO
    END DO
  END DO
  CALL scmoutput(riout,'dvdz',                                          &
       'Shear in Ri',' ',                                               &
       t_avg,d_bl,default_streams,'',routinename)

  DO k=2, bl_levels
    DO j=pdims%j_start, pdims%j_end
      DO i=pdims%i_start, pdims%i_end
        ! so OUT(K) is on theta-level K
        riout(i,j,k-1) = elm(i,j,k)
      END DO
    END DO
  END DO
  CALL scmoutput(riout,'elm',                                           &
       'Mixing length for momentum','m',                                &
       t_avg,d_bl,default_streams,'',routinename)

  CALL scmoutput(zhpar,'zhpar',                                         &
       'Height of top of parcel ascent','m',                            &
       t_avg,d_sl,default_streams,'',routinename)

  CALL scmoutput(u_s,'ustar',                                           &
       'Explicit surface friction velocity','m/s',                      &
       t_avg,d_sl,default_streams,'',routinename)

END IF ! scmdiag_bl /model_type

! interpolate rhokh_th to rho levels 2 to bl_levels

!$OMP  PARALLEL DEFAULT(SHARED)                                         &
!$OMP& PRIVATE(i, j, k, weight1, weight2, weight3, lambdah, z_scale,    &
!$OMP&         vkz, f_log)
!$OMP  DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      weight1 = r_theta_levels(i,j,k) -                                 &
              r_theta_levels(i,j, k-1)
      weight2 = r_theta_levels(i,j,k) -                                 &
              r_rho_levels(i,j,k)
      weight3 = r_rho_levels(i,j,k) -                                   &
              r_theta_levels(i,j,k-1)
      IF ( k  ==  bl_levels ) THEN
              ! assume rhokh_th(BL_LEVELS+1) is zero
        rhokh(i,j,k) = ( weight2/weight1 ) * rhokh_th(i,j,k)
        IF (blending_option /= off) weight_1dbl_rho(i,j,k) =            &
                               (weight2/weight1) * weight_1dbl(i,j,k)
      ELSE
        rhokh(i,j,k) = (weight3/weight1) * rhokh_th(i,j,k+1)            &
                     + (weight2/weight1) * rhokh_th(i,j,k)
        IF (blending_option /= off) weight_1dbl_rho(i,j,k) =            &
                              (weight3/weight1)*weight_1dbl(i,j,k+1)    &
                            + (weight2/weight1)*weight_1dbl(i,j,k)
      END IF

      IF (local_fa == free_trop_layers) THEN
        ! elh already included in rhokh_th so no need to calculate
        ! here, but interpolate elh separately for diagnostic
        IF (BL_diag%l_elh3d) THEN
          IF ( k  ==  bl_levels ) THEN
            ! assume rhokh_th(BL_LEVELS+1) is zero
            elh_rho(i,j,k) = ( weight2/weight1 ) * elh(i,j,k)
          ELSE
            elh_rho(i,j,k) =                                            &
              weight3/weight1 *                                         &
                      elh(i,j,k+1)                                      &
             +weight2/weight1 *                                         &
                      elh(i,j,k)
          END IF
          BL_diag%elh3d(i,j,k)=elh_rho(i,j,k)
        END IF
      ELSE
        IF ((sbl_op/=equilibrium_sbl) .OR. (fb_surf(i,j) >  0.0)) THEN
          !------------------------------------------------------------
          !  Include mixing length, ELH, in RHOKH.
          !  Code moved from EX_COEF to avoid interpolation
          !------------------------------------------------------------
          IF ( k >= ntml_local(i,j)+2 .AND. l_full_lambdas .AND.        &
               local_fa == to_sharp_across_1km ) THEN
            ! Assuming only LOCAL_FA = "to_sharp_across_1km" option
            ! will have L_FULL_LAMBDAS.
            ! If other LOCAL_FA options are coded here then
            ! changes must be included in section 2.1 of ex_coef
            IF (l_rp2) THEN
              lambdah = MAX ( lambda_min , par_mezcla*zh_local(i,j) )
            ELSE
              lambdah = MAX ( lambda_min , 0.15*zh_local(i,j) )
            END IF
            z_scale = 1000.0
            weight1 = 0.5*( 1.0 -                                       &
                        TANH(3.0*((z_uv(i,j,k)/z_scale )-1.0) ) )
            lambdah = lambdah * weight1                                 &
                         + lambda_min*( 1.0 -  weight1)
            ! no need to do log profile correction as klog_layr eq 2
            vkz = vkman * ( z_uv(i,j,k) + z0m_eff_gb(i,j) )
            elh_rho(i,j,k) = vkz / (1.0 + vkz/lambdah )
          END IF
          ! Reinstate UKV drainage flow bug here, where lambdah was not
          ! enhanced as intended (and as was done in ex_coef)!
          IF (sg_orog_mixing == sg_shear_enh_lambda) THEN
            IF (l_rp2) THEN
              lambdah = MAX ( lambda_min , par_mezcla*zh_local(i,j) )
            ELSE
              lambdah = MAX ( lambda_min , 0.15*zh_local(i,j) )
            END IF
            IF (k >= ntml_local(i,j)+2 .AND. .NOT. l_full_lambdas) THEN
              lambdah = lambda_min
            END IF
            IF (k <= k_log_layr) THEN
              vkz   = vkman * ( z_tq(i,j,k) - z_tq(i,j,k-1) )
              f_log = LOG( ( z_tq(i,j,k) + z0m_eff_gb(i,j)   ) /        &
                           ( z_tq(i,j,k-1) + z0m_eff_gb(i,j) ) )
              elh_rho(i,j,k) = vkz / ( f_log + vkz / lambdah )
            ELSE
              vkz = vkman * ( z_uv(i,j,k) + z0m_eff_gb(i,j) )
              elh_rho(i,j,k) = vkz / (1.0 + vkz/lambdah )
            END IF
          END IF
          ! End of UKV bug!

          IF (BL_diag%l_elh3d) BL_diag%elh3d(i,j,k)=elh_rho(i,j,k)

          rhokh(i,j,k) = elh_rho(i,j,k) * rhokh(i,j,k)

        END IF   ! test on sbl_op
      END IF   ! test on local_fa = free_trop_layers

          ! Finally multiply RHOKH by dry density
      IF (l_mr_physics)                                                 &
         rhokh(i,j,k) = rho_mix(i,j,k) * rhokh(i,j,k)

    END DO
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
      !----------------------------------------------------------
      ! Use local NTML if significantly higher (to allow for
      ! local averaging) than the non-local or if the non-local
      ! is on the ground (=1)
      !----------------------------------------------------------
    IF ( .NOT. cumulus(i,j) .AND.                                       &
              ( ntml_local(i,j)  >   ntml_nl(i,j)+1                     &
                .OR. ntml_nl(i,j)  ==  1 )            ) THEN
      ntml(i,j) = ntml_local(i,j)
      sml_disc_inv(i,j) = 0   ! reset flag for subgrid inversion
    ELSE
      ntml(i,j) = ntml_nl(i,j)
    END IF
      !----------------------------------------------------------
      ! If local NTML is higher than NTDSC then ignore DSC layer
      ! for diagnostics but keep mixing associated with it
      !----------------------------------------------------------
    IF ( ntml_local(i,j)  >   ntdsc(i,j)+1 ) THEN
      dsc_disc_inv(i,j) = 0
      ntdsc(i,j) = 0
      nbdsc(i,j) = 0
      zhsc(i,j)  = 0.0
      dsc(i,j)   = .FALSE.
      coupled(i,j) = .FALSE.
    END IF
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

! Calculate max of two coeffs
CALL kmkh (                                                             &
! IN data
   bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                              &
   ntml,cumulus,ntdsc,sml_disc_inv,dsc_disc_inv,                        &
   weight_1dbl, weight_1dbl_rho,                                        &
! INOUT data
   rhokm,rhokh,rhokmz,rhokhz,rhokm_top,rhokh_top                        &
   )

IF (BL_diag%l_weight1d) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,pdims,BL_diag,weight_1dbl)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%weight1d(i,j,k)=weight_1dbl(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF (BL_diag%l_rhokm) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,pdims,BL_diag,rhokm)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokm(i,j,k)=rhokm(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

IF (BL_diag%l_rhokh) THEN
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,pdims,BL_diag,rhokh)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        BL_diag%rhokh(i,j,k)=rhokh(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Calculation of TKE diagnostic.
! Stored on theta-levels with TKE(K) on theta-level(k-1),
! consistent with RHOKM(K), RI(K), etc.
! The K=1 value could be set to a diagnosed surface value (eg as a
! function of ustar, wstar) but is currently just set to zero
! BL_diag%tke currently contains (the reciprocal of) the non-local
! (SML and DSC) mixed layer timescale, calculated in excf_nl

IF (BL_diag%l_tke) THEN

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED(BL_diag, rho_wet_tq, rhokm, dbdz, bl_levels, pdims)       &
!$OMP& PRIVATE(i, j, k, recip_time_sbl, recip_time_cbl)

  DO k = 2, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end

        ! stable timescale
        recip_time_sbl = SQRT( MAX(dbdz(i,j,k),1.0e-5) )/0.7
        recip_time_cbl = BL_diag%tke(i,j,k)

        ! TKE diagnostic - taking 5 m2/s2 as a suitable maximum
        BL_diag%tke(i,j,k)=  MIN( max_tke,                              &
              ( rhokm(i,j,k) / rho_wet_tq(i,j,k-1) )*(                  &
                 recip_time_cbl + recip_time_sbl ) )

      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  IF ( l_subgrid_qcl_mp .OR. l_wvar_for_conv ) THEN
    ! Set bl_w_var to mimimum value.
    ! Prevents any unset values in the prognostic
    ! or any very low values being passed through
    ! to the microphysics turbulence call
    bl_w_var(:,:,:) = 1.0e-12

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,pdims,BL_diag,bl_w_var)
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          IF (BL_diag%tke(i,j,k) > 1.0e-12) THEN
            bl_w_var(i,j,k) = BL_diag%tke(i,j,k)
          END IF
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO

  END IF ! l_subgrid_qcl_mp .OR. l_wvar_for_conv

  ! at this point, BL_diag%tke really contains sigma_w^2. To make it look
  ! a bit more like TKE near the surface, we will keep it constant below
  ! the max of rhokm_surf
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED( pdims, rhokmz, BL_diag ) PRIVATE(i, j, k, kmax )
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      kmax = MAXLOC(rhokmz(i,j,:))
      DO k = 2, kmax(1)
        BL_diag%tke(i,j,k) = BL_diag%tke(i,j,kmax(1))
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF    ! BL_diag%L_tke

IF (blending_option /= off .AND. l_blend_isotropic) THEN
  !   ! Blended diffusion coefficients now held in rhokm and rhokh
  !   ! so copy to visc_m,h for horizontal diffusion too.
  !   ! Need to interpolate rhokh back to theta levels for visc_h

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, weight1, weight2,      &
!$OMP& weight3)

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF ( blending_option == blend_except_cu .AND.                   &
             cumulus(i,j) .AND. ntdsc(i,j) == 0) THEN
          ! pure cumulus layer so revert to Smag scheme
          visc_m(i,j,k) = visc_m(i,j,k)                                 &
                           *rneutml_sq(i,j,k)*fm_3d(i,j,k+1)
          visc_h(i,j,k) = visc_h(i,j,k)                                 &
                           *rneutml_sq(i,j,k)*fh_3d(i,j,k+1)
          IF (k == bl_levels-1) THEN
            ! need to do bl_levels too apparently
            ! (see non-blend code below)
            visc_m(i,j,k+1) = visc_m(i,j,k+1)*rneutml_sq(i,j,k+1)
            visc_h(i,j,k+1) = visc_h(i,j,k+1)*rneutml_sq(i,j,k+1)
          END IF
        ELSE
          visc_m(i,j,k) = rhokm(i,j,k+1)/rho_wet_tq(i,j,k)
        END IF
      END DO
    END DO
  END DO
!$OMP END DO
!$OMP DO SCHEDULE(STATIC)
  DO k = 2, bl_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        IF (.NOT. ( blending_option == blend_except_cu .AND.            &
             cumulus(i,j) .AND. ntdsc(i,j) == 0)) THEN
          ! NOT a pure cumulus layer or blending_option NE
          ! blend_except_cu, so use standard blending
          weight1 = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
          weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
          weight3 = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
          visc_h(i,j,k) = ( weight2*(rhokh(i,j,k+1)/rho_mix(i,j,k+1))   &
                          + weight3*(rhokh(i,j,k)  /rho_mix(i,j,k)  ))  &
                                         / weight1
        END IF
      END DO
    END DO
  END DO
!$OMP END DO
  k = 1
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF (.NOT. ( blending_option == blend_except_cu .AND.              &
             cumulus(i,j) .AND. ntdsc(i,j) == 0)) THEN
        ! NOT a pure cumulus layer or blending_option NE
        ! blend_except_cu, so use standard blending
        visc_h(i,j,k) = rhokh_th(i,j,k+1)/rho_wet_tq(i,j,k)
      END IF
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

ELSE IF (l_subfilter_horiz .OR. l_subfilter_vert) THEN

  ! visc_m,h on IN are just S and visc_m,h(k) are co-located with w(k)

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,pdims,visc_m,rneutml_sq,visc_h,fh_3d,fm_3d,max_diff)

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        visc_m(i,j,k) = visc_m(i,j,k)*rneutml_sq(i,j,k)
        visc_h(i,j,k) = visc_h(i,j,k)*rneutml_sq(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels-1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        ! stability functions are indexed with Ri, fm(k) on w(k-1)
        visc_h(i,j,k) = visc_h(i,j,k)*fh_3d(i,j,k+1)
        visc_m(i,j,k) = visc_m(i,j,k)*fm_3d(i,j,k+1)
        ! APL why apply this cap here for implicit vertical diffusion?
        ! (also applied in atm_step_phys_init for horiz diffn,
        ! that actually needs it)?
        visc_h(i,j,k) = MIN(visc_h(i,j,k),max_diff(i,j))
        visc_m(i,j,k) = MIN(visc_m(i,j,k),max_diff(i,j))
      END DO
    END DO
  END DO
!$OMP END DO

!$OMP END PARALLEL

  ! visc_m and visc _h are now lambda^2*S*FM and lambda^2*S*FH

  IF (l_subfilter_vert .AND. blending_option == off) THEN

    ! visc_h_rho(k) is held on rho(k), same as BL's rhokh
    ALLOCATE (visc_h_rho(pdims%i_start:pdims%i_end,                     &
                         pdims%j_start:pdims%j_end, bl_levels))

!$OMP  PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k, weight1, weight2,      &
!$OMP& weight3)

!$OMP DO SCHEDULE(STATIC)
    DO k = 2, bl_levels
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          weight1 = r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)
          weight2 = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
          weight3 = r_rho_levels(i,j,k)   - r_theta_levels(i,j,k-1)
          IF ( k  ==  bl_levels ) THEN
            ! assume visc_h(bl_levels) is zero (Ri and thence f_h not defined)
            visc_h_rho(i,j,k) = ( weight2/weight1 ) * visc_h(i,j,k-1)
          ELSE
            visc_h_rho(i,j,k) = ( weight3/weight1 ) * visc_h(i,j,k)     &
                              + ( weight2/weight1 ) * visc_h(i,j,k-1)
          END IF
        END DO
      END DO
    END DO
!$OMP END DO

!$OMP END PARALLEL

    ! Overwrite the diffusion coefficients from the BL scheme
    !(RHOKM and RHOKH) with those obtained from the Smagorinsky scheme

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(k,j,i)                                                          &
!$OMP SHARED(bl_levels,turb_startlev_vert,turb_endlev_vert,pdims,rhokm,       &
!$OMP        visc_m,rho_wet_tq,rhokh,visc_h_rho,rho_mix)
    DO k = 2, bl_levels
      IF (k >= turb_startlev_vert .AND.                                 &
          k <= turb_endlev_vert) THEN
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            rhokm(i,j,k) = visc_m(i,j,k-1)*rho_wet_tq(i,j,k-1)
            rhokh(i,j,k) = visc_h_rho(i,j,k)*rho_mix(i,j,k)
          END DO
        END DO
      ELSE
        DO j = pdims%j_start, pdims%j_end
          DO i = pdims%i_start, pdims%i_end
            rhokm(i,j,k) = 0.0
            rhokh(i,j,k) = 0.0
          END DO
        END DO
      END IF
    END DO
!$OMP END PARALLEL DO
    DEALLOCATE (visc_h_rho)

  END IF ! L_subfilter_vert but not blend
END IF ! L_subfilter_horiz or L_subfilter_vert or isotropic blending_option

!-----------------------------------------------------------------------
! Diagnose boundary layer type.
!      Seven different types are considered:
!      1 - Stable b.l.
!      2 - Stratocumulus over a stable surface layer.
!      3 - Well mixed buoyancy-driven b.l. (possibly with stratocumulus)
!      4 - Decoupled stratocumulus (not over cumulus).
!      5 - Decoupled stratocumulus over cumulus.
!      6 - Cumulus capped b.l.
!      7 - Shear-dominated unstable b.l.
!-----------------------------------------------------------------------
!      First initialise the type variables and set the depth diagnostics

! Top of surface mixed layer (Ksurf profile)
IF (BL_diag%l_smltop) THEN
  BL_diag%smltop=zhnl
END IF
! Top of decoupled stratocu layer
IF (BL_diag%l_dsctop) THEN
  BL_diag%dsctop=zhsc
END IF
! Height of diagnosis parcel top
IF (BL_diag%l_zhpar) THEN
  BL_diag%zhpar=zhpar
END IF

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)
!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ! Max height of BL turbulent mixing
    zht(i,j) = MAX( zhnl(i,j) , zhsc(i,j) )
    ! PBL depth diagnostic: start with top of non-local BL
    zh(i,j) = zhnl(i,j)
    IF ( ntml(i,j)  >   ntml_nl(i,j) ) THEN
      ! Higher local K allowed so reset ZH, ZHT diagnostics
      zh(i,j)  = MAX( zh(i,j) , zh_local(i,j) )
      zht(i,j) = MAX( zht(i,j), zh_local(i,j) )
      IF (.NOT. l_fix_zh) THEN
        IF (forced_cu == off) THEN
          ! going to use zhnl in tr_mix and ex_flux_uv which
          ! was formerly spuriously set to zh
          zhnl(i,j)=zh(i,j)
        ELSE
          ! spuriously overwriting zh with zhnl
          zh(i,j)=zhnl(i,j)
        END IF
      END IF
    END IF
    bl_type_1(i,j) = 0.0
    bl_type_2(i,j) = 0.0
    bl_type_3(i,j) = 0.0
    bl_type_4(i,j) = 0.0
    bl_type_5(i,j) = 0.0
    bl_type_6(i,j) = 0.0
    bl_type_7(i,j) = 0.0
  END DO
END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF (dynamic_bl_diag(i,j)) THEN
      !       ! shear-dominated, via iDynDiag option
      bl_type_7(i,j) = 1.0
    ELSE
      IF (.NOT. unstable(i,j) .AND. .NOT. dsc(i,j) .AND.                &
          .NOT. cumulus(i,j)) THEN
        !         ! Stable b.l.
        bl_type_1(i,j) = 1.0
      ELSE IF (.NOT. unstable(i,j) .AND. dsc(i,j) .AND.                 &
               .NOT. cumulus(i,j)) THEN
        !         ! Stratocumulus over a stable surface layer
        bl_type_2(i,j) = 1.0
      ELSE IF (unstable(i,j) .AND. .NOT. cumulus(i,j) .AND.             &
              .NOT. dsc(i,j) ) THEN
        !         ! Well mixed b.l. (possibly with stratocumulus)
        IF ( ntml(i,j)  >   ntml_nl(i,j) ) THEN
          ! shear-dominated - currently identified
          ! by local NTML overriding non-local
          bl_type_7(i,j) = 1.0
        ELSE
          ! buoyancy-dominated
          bl_type_3(i,j) = 1.0
        END IF
      ELSE IF (unstable(i,j) .AND. dsc(i,j) .AND.                       &
                                      .NOT. cumulus(i,j)) THEN
        !         ! Decoupled stratocumulus (not over cumulus)
        bl_type_4(i,j) = 1.0
      ELSE IF (dsc(i,j) .AND. cumulus(i,j)) THEN
        !         ! Decoupled stratocumulus over cumulus
        bl_type_5(i,j) = 1.0
      ELSE IF (.NOT. dsc(i,j) .AND. cumulus(i,j)) THEN
        !         ! Cumulus capped b.l.
        bl_type_6(i,j) = 1.0
      END IF
    END IF

  END DO
END DO
!$OMP END DO NOWAIT

!-----------------------------------------------------------------------
! 5.5 Calculation of explicit fluxes of T,Q
!-----------------------------------------------------------------------
!$OMP MASTER
CALL ex_flux_tq (                                                       &
! IN levels etc
    bl_levels,nSCMDpkgs,L_SCMDiags,                                     &
! IN fields
    tl,qw,rdz_charney_grid, rhokh, rhokhz,grad_t_adj,grad_q_adj,        &
    rhof2, rhofsc, ft_nt, fq_nt, ft_nt_dscb, fq_nt_dscb, tothf_zh,      &
    tothf_zhsc, totqf_zh, totqf_zhsc, weight_1dbl_rho,                  &
    ntml_nl, ntdsc, nbdsc,                                              &
! INOUT fields
    ftl,fqw                                                             &
    )
!$OMP END MASTER

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      f_ngstress(i,j,k) = weight_1dbl(i,j,k) * f_ngstress(i,j,k)
    END DO
  END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------

IF ( l_scmdiags(scmdiag_bl) .AND.                                       &
     model_type == mt_single_column ) THEN

  DO j=pdims%j_start, pdims%j_end
    jScm = j - pdims%j_start + 1
    DO i=pdims%i_start, pdims%i_end
      iScm = i - pdims%i_start + 1
      TmpScm2d(iScm,jScm)= REAL(sml_disc_inv(i,j))
    END DO
  END DO
  CALL scmoutput(TmpScm2d,'SML_DISC_INV',                               &
       'Indicator for subgrid SML inversion','Indicator',               &
       t_inst,d_point,default_streams,'',routinename)

  DO j=pdims%j_start, pdims%j_end
    jScm = j - pdims%j_start + 1
    DO i=pdims%i_start, pdims%i_end
      iScm = i - pdims%i_start + 1
      TmpScm2d(iScm,jScm)= REAL(dsc_disc_inv(i,j))
    END DO
  END DO
  CALL scmoutput(TmpScm2d,'DSC_DISC_INV',                               &
       'Indicator for subgrid DSC inversion','Indicator',               &
       t_inst,d_point,default_streams,'',routinename)

END IF ! scmdiag_bl / model_type

!-----------------------------------------------------------------------
! 5.6.1 Calculate explicit surface fluxes of U and V on
!       P-grid for convection scheme
!-----------------------------------------------------------------------
!$OMP PARALLEL                                                                &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,uw0,rhokm,u_p,u_0_px,vw0,v_p,v_0_px,wstar,wthvs,           &
!$OMP        cu_over_orog,cumulus,fb_surf,zh,g,bt,l_param_conv,ntml,ntml_nl,  &
!$OMP        formdrag,tau_fd_x,tau_fd_y,ntdsc,BL_diag)

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    uw0(i,j) = -rhokm(i,j,1) *                                          &
                      ( u_p(i,j,1) - u_0_px(i,j) )
    vw0(i,j) = -rhokm(i,j,1) *                                          &
                      ( v_p(i,j,1) - v_0_px(i,j) )
  END DO
END DO
!$OMP END DO NOWAIT
IF (formdrag ==  explicit_stress) THEN
!$OMP DO SCHEDULE(STATIC)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      uw0(i,j) = uw0(i,j) - tau_fd_x(i,j,1)
      vw0(i,j) = vw0(i,j) - tau_fd_y(i,j,1)
    END DO
  END DO
!$OMP END DO NOWAIT
END IF

!-----------------------------------------------------------------------
! 5.7 Set NTML to max number of turbulently mixed layers
!      Calculate quantities to pass to convection scheme.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    wstar(i,j) = 0.0
    wthvs(i,j) = 0.0
    cu_over_orog(i,j) = 0.0
    IF ( cumulus(i,j) ) THEN
      IF ( fb_surf(i,j)  >   0.0 ) THEN
        wstar(i,j) = ( zh(i,j)*fb_surf(i,j) )**one_third
        wthvs(i,j) = fb_surf(i,j) / ( g * bt(i,j,1) )
      END IF
      wstar(i,j) = MAX( 0.1, wstar(i,j) )
      IF (.NOT. l_param_conv) THEN
        ntml(i,j) = MAX( 2, ntml_nl(i,j) - 1 )
      END IF
    ELSE
      ntml(i,j) = MAX( ntml_nl(i,j) , ntdsc(i,j) )
    END IF
    ! Limit explicitly calculated surface stresses
    ! to a physically plausible level.
    IF ( uw0(i,j)  >=  5.0 ) THEN
      uw0(i,j) =  5.0
    ELSE IF ( uw0(i,j)  <=  -5.0 ) THEN
      uw0(i,j) = -5.0
    END IF
    IF ( vw0(i,j)  >=  5.0 ) THEN
      vw0(i,j) =  5.0
    ELSE IF ( vw0(i,j)  <=  -5.0 ) THEN
      vw0(i,j) = -5.0
    END IF
    IF (BL_diag%l_wstar .AND. (fb_surf(i,j) >0.0)) THEN
      BL_diag%wstar(i,j)= (zh(i,j)*fb_surf(i,j))**one_third
    END IF
  END DO
END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL

IF (l_param_conv) THEN

  ! Check for CUMULUS having been diagnosed over steep orography.
  ! Reset to false but keep NTML at NLCL (though decrease by 2 so that
  ! coupling between BL and convection scheme can be maintained).
  ! Reset type diagnostics.

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(l,j,i)                                                          &
!$OMP SHARED(land_pts,land_index,pdims,cumulus,ho2r2_orog,l_shallow,bl_type_5,&
!$OMP        bl_type_6,cu_over_orog,ntml)
  DO l = 1, land_pts
    j=(land_index(l)-1)/pdims%i_end + 1
    i=land_index(l) - (j-1)*pdims%i_end
    IF (cumulus(i,j) .AND. ho2r2_orog(l)  >   900.0) THEN
      cumulus(i,j)      = .FALSE.
      l_shallow(i,j)    = .FALSE.
      bl_type_5(i,j)    = 0.0
      bl_type_6(i,j)    = 0.0
      cu_over_orog(i,j) = 1.0
      IF (ntml(i,j)  >=  3) ntml(i,j) = ntml(i,j) - 2
    END IF
  END DO
!$OMP END PARALLEL DO

  ! Check that CUMULUS and L_SHALLOW are still consistent
  ! and reset ntml back to where it was originally diagnosed
  ! if cumulus is still true in order to trigger convection
  ! at the correct level

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,cumulus,l_shallow,ntml,ntml_save)
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( .NOT. cumulus(i,j) ) l_shallow(i,j) = .FALSE.
      IF ( cumulus(i,j) ) ntml(i,j) = ntml_save(i,j)
    END DO
  END DO
!$OMP END PARALLEL DO

END IF    ! (l_param_conv)
!-----------------------------------------------------------------------
!     Set shallow convection diagnostic: 1.0 if L_SHALLOW (and CUMULUS)
!                                        0.0 if .NOT. CUMULUS
!-----------------------------------------------------------------------

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(NONE)                                                           &
!$OMP PRIVATE(j,i)                                                            &
!$OMP SHARED(pdims,cumulus,l_shallow,shallowc)
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( cumulus(i,j) .AND. l_shallow(i,j) ) THEN
      shallowc(i,j) = 1.0
    ELSE
      shallowc(i,j) = 0.0
    END IF
  END DO
END DO
!$OMP END PARALLEL DO

! ------rhcrit parametrization---------------------
! a_qs = a in documentation on theta level K
! a_dqsdt = b/exner in documentation on theta level K
! rho_wet_tq = rho on theta level K
! qsw = qsat(Tl) on theta level K-1
! dsldzm is on theta-level K-1
! dqwdzm is on theta-level K-1
! elm is on theta-level K-1
! rhokh is on rho-level K
! tke is on theta-level K-1
! exner is on theta-level K-1
IF (i_rhcpt == rhcpt_tke_based) THEN
  b2 = 15.0
  root6 = SQRT(6.0)

  ! level 2 needs special treatment because of the surface
  k = 2

!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(j,i,sh,exner,sgm,delta_x,qsw_arr)
  DO j = tdims%j_start, tdims%j_end

    IF ( l_new_qsat_bl ) THEN
      IF ( l_mr_physics ) THEN
        CALL qsat_wat_mix_new(qsw_arr,tl(:,j,k-1),p_theta_levels(:,j,k-1),    &
                              tdims%i_len)
      ELSE
        CALL qsat_wat_new(qsw_arr,tl(:,j,k-1),p_theta_levels(:,j,k-1),        &
                              tdims%i_len)
      END IF
    ELSE
    !DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsw_arr,tl(:,j,k-1),p_theta_levels(:,j,k-1),          &
                        tdims%i_len,l_mr_physics)
    END IF

    DO i = tdims%i_start, tdims%i_end
      ! calculate sh, don't interpolate with the surface flux or divide by 0
      IF ( BL_diag%tke(i,j,k) > 1.0e-10 ) THEN
        sh = rhokh(i,j,k)                                               &
             / ( rho_wet_tq(i,j,k-1) * elm(i,j,k)                       &
             * SQRT( 2.0 * BL_diag%tke(i,j,k) ) )
      ELSE
        sh = 0.01
      END IF
      IF ( sh < 0.01 ) sh = 0.01
      ! calculate exner
      exner = ( p_theta_levels(i,j,k-1) / pref )**kappa
      ! calculate the variance, use gradient interpolated between levs 1 and 2
      sgm = a_dqsdt(i,j,k-1)**2 * exner**2 * b2 * sh                    &
            * elm(i,j,k)**2 * dsldz(i,j,k)**2                           &
          + a_qs(i,j,k-1)**2 * b2 * sh                                  &
            * elm(i,j,k)**2 * dqwdz(i,j,k)**2                           &
          - 2.0 * a_qs(i,j,k-1) * a_dqsdt(i,j,k-1) * exner * b2         &
            * sh * elm(i,j,k)**2 * dsldz(i,j,k) * dqwdz(i,j,k)
      ! do this for safety, not sure if it's really needed
      sgm = SQRT ( MAX( sgm, 0.0 ) )

      ! calculate rhcrit, with appropriate limits
      ! calculate grid-box size, just take surface for simplicity
      delta_x = SQRT( r_theta_levels(i,j,k-1) * delta_lambda *          &
                      r_theta_levels(i,j,k-1) * delta_phi *             &
                      cos_theta_latitude(i,j) )
      ! max limit, based on curve fitted to aircraft observations
      max_rhcpt(i,j) = MIN( 0.99, 0.997 - 0.0078 *                      &
                                  LOG( 0.001 * delta_x ) )
      ! min limit, based on curve fitted to aircraft observations
      min_rhcpt(i,j) = MAX( 0.6, 0.846 - 0.065 *                        &
                                 LOG( 0.001 * delta_x ) )
      ! full expression
      rhcpt(i,j,k-1) = MIN( max_rhcpt(i,j), MAX( min_rhcpt(i,j),        &
                        1.0 - root6 * sgm / (a_qs(i,j,k-1) * qsw_arr(i))))
    END DO !i
  END DO   !j
!$OMP END PARALLEL DO


  ! remaining bl levels use proper interpolation

!$OMP PARALLEL                                                                &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(k,j,i,weight1,weight2,weight3,sh,exner,sgm,qsw_arr)

  DO k = 3, bl_levels-1
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end

      IF ( l_new_qsat_bl ) THEN
        IF ( l_mr_physics ) THEN
          CALL qsat_wat_mix_new(qsw_arr,tl(:,j,k-1),p_theta_levels(:,j,k-1),  &
                                tdims%i_len)
        ELSE
          CALL qsat_wat_new(qsw_arr,tl(:,j,k-1),p_theta_levels(:,j,k-1),      &
                                tdims%i_len)
        END IF
      ELSE
      !DEPENDS ON: qsat_wat_mix
        CALL qsat_wat_mix(qsw_arr,tl(:,j,k-1),p_theta_levels(:,j,k-1),        &
                          tdims%i_len,l_mr_physics)
      END IF

      DO i = tdims%i_start, tdims%i_end
        weight1 = r_rho_levels(i,j,k) -                                 &
                  r_rho_levels(i,j,k-1)
        dsldzm(i,j,k) = dsldzm(i,j,k) / weight1
        dqwdzm(i,j,k) = dqwdzm(i,j,k) / weight1
        ! calculate sh, don't divide by 0
        IF ( BL_diag%tke(i,j,k) > 1.0e-10 ) THEN
          weight2 = r_theta_levels(i,j,k-1)-                            &
                    r_rho_levels(i,j,k-1)
          weight3 = r_rho_levels(i,j,k) -                               &
                    r_theta_levels(i,j,k-1)
          sh = ( weight2 * rhokh(i,j,k) + weight3 * rhokh(i,j,k-1) )    &
               / ( weight1 * rho_wet_tq(i,j,k-1) * elm(i,j,k)           &
               * SQRT( 2.0 * BL_diag%tke(i,j,k) ) )
        ELSE
          sh = 0.01
        END IF
        IF ( sh < 0.01 ) sh = 0.01
        ! calculate exner
        exner = ( p_theta_levels(i,j,k-1) / pref )**kappa
        ! calculate the variance
        sgm = a_dqsdt(i,j,k-1)**2 * exner**2 * b2 * sh                  &
              * elm(i,j,k)**2 * dsldzm(i,j,k)**2                        &
            + a_qs(i,j,k-1)**2 * b2 * sh                                &
              * elm(i,j,k)**2 * dqwdzm(i,j,k)**2                        &
            - 2.0 * a_qs(i,j,k-1) * a_dqsdt(i,j,k-1) * exner * b2       &
              * sh * elm(i,j,k)**2 * dsldzm(i,j,k) * dqwdzm(i,j,k)
        ! do this for safety, not sure if it's really needed
        sgm = SQRT ( MAX( sgm, 0.0 ) )
        ! calculate rhcrit, with appropriate limits
        rhcpt(i,j,k-1) = MIN( max_rhcpt(i,j), MAX( min_rhcpt(i,j),      &
                        1.0 - root6 * sgm / (a_qs(i,j,k-1) * qsw_arr(i))))
      END DO !i
    END DO   !j
!$OMP END DO NOWAIT
  END DO     !k
!$OMP END PARALLEL

  ! just use 0.8 above bl_levels
!$OMP PARALLEL DO                                                             &
!$OMP SCHEDULE(STATIC)                                                        &
!$OMP DEFAULT(SHARED)                                                         &
!$OMP PRIVATE(k,j,i)
  DO k = bl_levels-1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        rhcpt(i,j,k) = 0.8
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF !i_rhcpt

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE bdy_expl2
END MODULE bdy_expl2_mod
