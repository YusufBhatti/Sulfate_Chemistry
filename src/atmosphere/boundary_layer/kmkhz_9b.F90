! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate the non-local turbulent mixing
!           coefficients KM and KH

!  Programming standard: UMDP3

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Boundary Layer
!---------------------------------------------------------------------
MODULE kmkhz_9b_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'KMKHZ_9B_MOD'
CONTAINS

SUBROUTINE kmkhz_9b (                                                   &
! IN levels/switches
 bl_levels,BL_diag,nSCMDpkgs,L_SCMDiags,                                &
! IN fields
 p,rho_wet_tq,rho_mix,rho_mix_tq,t,q,qcl,qcf,cf,qw,tl,dzl,rdz,z_tq,z_uv,&
 rad_hr,                                                                &
 bt,bq,btm,bqm,dqsdt,btm_cld,bqm_cld,a_qs,a_qsm,a_dqsdtm,               &
 v_s,fb_surf,rhostar_gb,ntpar,zh_prev,                                  &
 zhpar,zmaxb_for_dsc,zmaxt_for_dsc,l_shallow,                           &
! INOUT fields
 ftl,fqw,zh,cumulus,ntml,w,t1_sd,q1_sd,                                 &
! OUT fields
 rhokm,rhokh,rhokm_top,rhokh_top,zhsc,                                  &
 unstable,dsc,coupled,sml_disc_inv,dsc_disc_inv,                        &
 ntdsc,nbdsc,f_ngstress, grad_t_adj, grad_q_adj,                        &
 kent, we_lim, t_frac, zrzi,                                            &
 kent_dsc, we_lim_dsc, t_frac_dsc, zrzi_dsc                             &
 )

USE atm_fields_bounds_mod, ONLY: pdims, pdims_s, tdims, tdims_l,        &
                                 ScmRowLen, ScmRow
USE bl_diags_mod, ONLY: strnewbldiag
USE bl_option_mod, ONLY: one_third
USE tuning_segments_mod, ONLY: bl_segment_size
USE cv_run_mod, ONLY: l_param_conv
USE gen_phys_inputs_mod, ONLY: l_mr_physics
USE model_domain_mod, ONLY: model_type, mt_single_column
USE missing_data_mod, ONLY: rmdi
USE planet_constants_mod, ONLY: cp, r, c_virtual, g,                    &
     etar, grcp, lcrcp, lsrcp
USE s_scmop_mod,  ONLY: default_streams, t_avg, d_bl, d_sl, scmdiag_bl
USE scmoutput_mod,ONLY: scmoutput
USE timestep_mod, ONLY: timestep
USE water_constants_mod, ONLY: tm

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_mix_new     => qsat_mix,                       &
                    l_new_qsat_bl !Currently defaults to FALSE

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! IN arguments
INTEGER, INTENT(IN) ::                                                  &
 bl_levels
                            ! IN No. of atmospheric levels for
                            !    which boundary layer fluxes are
                            !    calculated.

!     Declaration of new BL diagnostics.
TYPE (strnewbldiag), INTENT(INOUT) :: BL_diag

! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN) ::                                                  &
  nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL, INTENT(IN) ::                                                  &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

INTEGER, INTENT(IN) ::                                                  &
 ntpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! IN Top level of parcel ascent.
                            !    Used in convection scheme.
                            !    NOTE: CAN BE > BL_LEVELS-1

LOGICAL, INTENT(IN) ::                                                  &
 l_shallow(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! IN Flag to indicate shallow
                            !    convection
REAL, INTENT(IN) ::                                                     &
 bq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),     &
                            ! IN A buoyancy parameter for clear air
                            !    on p,T,q-levels (full levels).
 bt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),     &
                            ! IN A buoyancy parameter for clear air
                            !    on p,T,q-levels (full levels).
 bqm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                            ! IN A buoyancy parameter for clear air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 btm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                            ! IN A buoyancy parameter for clear air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 bqm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                            ! IN A buoyancy parameter for cloudy air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 btm_cld(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,           &
         bl_levels),                                                    &
                            ! IN A buoyancy parameter for cloudy air
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 a_qs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                            ! IN Saturated lapse rate factor
                            !    on p,T,q-levels (full levels).
 a_qsm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                            ! IN Saturated lapse rate factor
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 a_dqsdtm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,          &
          bl_levels),                                                   &
                            ! IN Saturated lapse rate factor
                            !    on intermediate levels (half levels):
                            !    (*,K) elements are k+1/2 values.
 p(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,0:bl_levels),    &
                            ! IN P(*,K) is pressure at full level k.
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),     &
                            ! IN Total water content (kg per kg air).
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),     &
                            ! IN Liquid/frozen water temperature (K).
 t(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),      &
                            ! IN Temperature (K).
 qcf(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,       &
     tdims_l%k_start:bl_levels),                                        &
                                   ! IN Cloud ice (kg per kg air)
 qcl(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,       &
     tdims_l%k_start:bl_levels),                                        &
                                   ! IN Cloud liquid water
 q(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,         &
   tdims_l%k_start:bl_levels),                                          &
                                   ! IN specific humidity

                            ! IN specific humidity
 cf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),     &
                            ! IN Cloud fractions for boundary levs.
 dzl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                            ! IN Layer depths (m).  DZL(,K) is the
                            !    distance from layer boundary K-1/2
                            !    to layer boundary K+1/2.  For K=1
                            !    the lower boundary is the surface.
 zh_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                            ! IN boundary layer height (m) from
                            !    previous timestep
 rdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                            ! IN Reciprocal of distance between
                            !    full levels (m-1).  1/RDZ(,K) is
                            !    the vertical distance from level
                            !    K-1 to level K, except that for
                            !    K=1 it is the height of the
                            !    lowest atmospheric full level.
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),   &
                            ! IN Z_tq(*,K) is the height of the
                            !    k-th full level above the surface.
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1), &
                            ! IN Z_uv(*,K) is the height of level
                            !       k-1/2 above the surface (m).
 rad_hr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
        2,bl_levels),                                                   &
                            ! IN (LW,SW) radiative heating rates
                            !    (K/s)
 v_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            ! IN Surface friction velocity (m/s)
 fb_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                             !IN Surface buoyancy flux over density
                            !       (m^2/s^3).
 dqsdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),  &
                            ! IN Partial derivative of QSAT w.r.t.
                            !    temperature.
 rho_mix(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
                            ! IN density on UV (ie. rho) levels,
                            !    used in RHOKH so dry density if
                            !    L_mr_physics is true
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                            ! IN density on TQ (ie. theta) levels,
                            !    used in RHOKM so wet density
 rho_mix_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,        &
            bl_levels),                                                 &
                            ! IN density on TQ (ie. theta) levels,
                            !    used in non-turb flux integration
                            !    so dry density if L_mr_physics is true
 rhostar_gb(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! IN Surface air density in kg per
                            !    cubic metre.
 zmaxb_for_dsc,                                                         &
 zmaxt_for_dsc,                                                         &
                            ! IN Max heights to look for DSC cloud
                            !    base and top
 zhpar(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! IN Height of top of NTPAR
                            !    NOTE: CAN BE ABOVE BL_LEVELS-1
! INOUT arrays
INTEGER, INTENT(INOUT) ::                                               &
 ntml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
    ! INOUT Number of model levels in the turbulently mixed layer.

LOGICAL, INTENT(INOUT) ::                                               &
 cumulus(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! INOUT Flag for Cu in the bl

REAL, INTENT(INOUT) ::                                                  &
 zh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                           ! INOUT Boundary layer height (m).
 fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                            ! INOUT "Explicit" flux of QW (i.e.
                            !        evaporation) from layer below
                            !        on P-grid (kg per sq m per s).
 ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),    &
                            ! INOUT "Explicit" flux of TL = H/CP
                            !       (sensible heat/CP) from layer
                            !       below, on P-grid.
 t1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                            ! INOUT Standard Deviation of level 1
                            !    temperature (K).
 q1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                            ! INOUT Standard Deviation of level 1
                            !    specific humidity (kg/kg).
 w(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,0:bl_levels)
                            ! INOUT Vertical velocity (m/s)

! OUT arrays
INTEGER, INTENT(OUT) ::                                                 &
 ntdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                            ! OUT Top level for turb mixing in
                            !     cloud layer
 nbdsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                            ! OUT Bottom level of any decoupled
                            !     turbulently mixed Sc layer.
 sml_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT Flags for whether discontinuous
 dsc_disc_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                            ! OUT inversions are diagnosed
 kent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                            ! OUT grid-level of SML inversion
 kent_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! OUT grid-level of DSC inversion

LOGICAL, INTENT(OUT) ::                                                 &
 unstable(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                            ! OUT Flag to indicate an unstable
                            !     surface layer.
 dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),              &
                            ! OUT Flag set if decoupled
                            !    stratocumulus layer found
 coupled(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                            ! OUT Flag to indicate Sc layer weakly
                            !     coupled to surface (ie weakly
                            !     decoupled)
REAL, INTENT(OUT) ::                                                    &
 grad_t_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Temperature gradient adjustment
                            !     for non-local mixing in unstable
                            !     turbulent boundary layer.
 grad_q_adj(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                            ! OUT Humidity gradient adjustment
                            !     for non-local mixing in unstable
                            !     turbulent boundary layer.
 rhokm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
       2:bl_levels),                                                    &
                            ! OUT Non-local turbulent mixing
                            !     coefficient for momentum.
 rhokh(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
       2:bl_levels),                                                    &
                            ! OUT Non-local turbulent mixing
                            !     coefficient for scalars.
 rhokm_top(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,         &
       2:bl_levels),                                                    &
                            ! OUT Top-down turbulent mixing
                            !     coefficient for momentum.
 rhokh_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
       2:bl_levels),                                                    &
                            ! OUT Top-down turbulent mixing
                            !     coefficient for scalars.
 f_ngstress(pdims_s%i_start:pdims_s%i_end,pdims_s%j_start:pdims_s%j_end,&
            2:bl_levels),                                               &
                            ! OUT dimensionless function for
                            !     non-gradient stresses
 zhsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),             &
                            ! OUT Cloud layer height (m).
 we_lim(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),         &
                            ! OUT rho*entrainment rate implied by
                            !     placing of subsidence
 zrzi(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),           &
                            ! OUT (z-z_base)/(z_i-z_base)
 t_frac(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),         &
                            ! OUT a fraction of the timestep
 we_lim_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),     &
                            ! OUT rho*entrainment rate implied by
                            !     placing of subsidence
 zrzi_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3),       &
                            ! OUT (z-z_base)/(z_i-z_base)
 t_frac_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,3)
                            ! OUT a fraction of the timestep

!----------------------------------------------------------------------
!    Local and other symbolic constants :-

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'KMKHZ_9B'

REAL :: a_plume,b_plume,a_grad_adj,max_t_grad,max_svl_grad,sc_cftol,    &
 ct_resid,svl_coup,dec_svl_grad,fgf
PARAMETER (                                                             &
 a_plume=0.2,                                                           &
 b_plume=3.26,                                                          &
 a_grad_adj=3.26,                                                       &
 max_t_grad=1.0e-3,                                                     &
 max_svl_grad=1.0e-3,                                                   &
                          ! maximum SVL gradient in a mixed layer
 dec_svl_grad=1.0e-3,                                                   &
                          ! SVL gradient required for weak decoupling
 sc_cftol=0.1,                                                          &
                          ! CF required for a Sc layer to be diagnosed
 ct_resid=200.0,                                                        &
                      ! Approx parcel cloud-top residence time (in s)
 svl_coup=0.5,                                                          &
                      ! Parameter controlling positioning of
                      ! surface-driven entrainment
 fgf=0.0)         ! Adiabatic gradient factor for ice

!  Define local storage.

!  (a) Workspace.

LOGICAL ::                                                              &
 cloud_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! Flag set when cloud base
                                 ! is reached.
 dsc_save(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
        ! Copy of DSC needed to indicate decoupling diagnosed in EXCF_NL
REAL ::                                                                 &
 qs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),     &
                                  ! Saturated sp humidity at pressure
                                 ! and temperature of sucessive
                                 ! levels.
 qcl_ic_bot(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                                 ! In-cloud liquid water content
                                 ! at the bottom of the model layer
 qcf_ic_bot(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                                 ! In-cloud frozen water content
                                 ! at the bottom of the model layer
 qcl_ic_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                                 ! In-cloud liquid water content
                                 ! at the top of the model layer.
 qcf_ic_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                                 ! In-cloud frozen water content
                                 ! at the top of the model layer.
 cfl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                                    ! Liquid cloud fraction.
 cff(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                                    ! Frozen cloud fraction.
 dqcldz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
                                 ! Vertical gradient of
                                 ! in-cloud liquid cloud water
                                 ! in a well-mixed layer.
 dqcfdz(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
                                 ! Vertical gradient of in-cloud
                                 ! frozen cloud water in a
                                 ! well-mixed layer.
 chi_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,             &
        2:bl_levels),                                                   &
                                 ! Mixing fraction of just saturate
                                 ! mixture.
 bflux_surf(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                 ! Buoyancy flux at the surface.
 bflux_surf_sat(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                 ! Saturated-air surface buoyancy flux
 db_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                 ! Buoyancy jump at the top of the
                                 ! boundary layer.
 db_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        2:bl_levels),                                                   &
                                 ! In-cloud buoyancy jump across
                                 ! layer interface.
 db(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,2:bl_levels),   &
                                     ! Buoyancy jump across layer
                                     !  interface.
 db_ksurf_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              2:bl_levels),                                             &
                                               ! Dry buoyancy jump in
                                 ! flux integral calculation (m/s2)
 db_ktop_dry(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             2:bl_levels),                                              &
                                              ! Sat. buoyancy jump in
                                 ! flux integral calculation (m/s2)
 db_ksurf_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              2:bl_levels),                                             &
                                               ! Dry buoyancy jump in
                                 ! flux integral calculation (m/s2)
 db_ktop_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,       &
             2:bl_levels),                                              &
                                              ! Sat. buoyancy jump in
                                 ! flux integral calculation (m/s2)
 df_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,        &
            bl_levels),                                                 &
                                 ! Radiative flux change over layer
                                 ! divided by c_P.
 dflw_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels),                                               &
                                 ! LW radiative flux change over layer
                                 ! divided by c_P.
 dfsw_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,      &
              bl_levels),                                               &
                                 ! SW radiative flux change over layer
                                 ! divided by c_P.
 tls_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
                                 ! Static energy increment due to
                                 !    large-scale vertical
                                 !    advection (K s^-1)
 qls_inc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,           &
         bl_levels),                                                    &
                                 ! Specific humidity increment
                                 !    due to large-scale
                                 !    vertical advection (s^-1)
 df_top_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                    ! Radiative flux change at cloud
                                    ! top divided by c_P.
 svl_plume(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    ! Liquid/frozen water virtual
                                    ! static energy over CP for a
                                    ! plume rising without dilution
                                    ! from level 1.
 sl_plume(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                    ! Liquid/frozen water
                                    ! static energy over CP for a
                                    ! plume rising without dilution
                                    ! from level 1.
 qw_plume(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                    ! QW for a plume rising without
                                    ! dilution from level 1.
 env_svl_km1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                    ! Density potential temperature
                                    ! for last layer considered
 z_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),            &
                                    ! Cloud fraction weighted
                                    ! thickness of b.l. cloud.
 bt_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! Buoyancy parameter at the top of
                                    ! the b.l.
 btt_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! In-cloud buoyancy parameter at
                                    ! the top of the b.l.
 btc_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! Cloud fraction weighted buoyancy
                                    ! parameter at the top of the b.l.
 db_top_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! In-cloud buoyancy jump at the
                                    ! top of the b.l.
 cld_factor(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! Fraction of grid box potentially
                                    ! giving evaporative entrainment.
 chi_s_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    ! Mixing fraction of just saturate
                                    ! mixture at top of the b.l.
 zeta_s(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! Non-cloudy fraction of mixing
                                    ! layer for surface forced
                                    ! entrainment term.
 zeta_r(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! Non-cloudy fraction of mixing
                                    ! layer for cloud top radiative
                                    ! cooling entrainment term.
 zc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),               &
                                    ! Cloud depth (not cloud fraction
                                    ! weighted).
 dscdepth(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                    ! Depth of cloud-layer (m)
 db_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! Buoyancy jump at the top of the
                                    ! surface mixed-layer.
 svl(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),    &
                                    ! Liquid/frozen water virtual
                                    ! temperature over CP.
 df_dsct_over_cp(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),  &
                                     ! Radiative flux change at DSC
                                    !  top divided by c_P.
 bt_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! Buoyancy parameter at the top of
                                    ! the DSC
 btt_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                    ! In-cloud buoyancy parameter at
                                    ! the top of the DSC
 btc_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),         &
                                    ! Cloud fraction weighted buoyancy
                                    ! parameter at the top of the DSC
 db_dsct_cld(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                    ! In-cloud buoyancy jump at the
                                    ! top of the DSC
 chi_s_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! Mixing fraction of just saturate
                                    ! mixture at top of the DSC
 zeta_r_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! Non-cloudy fraction of DSC
                                    ! for cloud top radiative
                                    ! cooling entrainment term.
 zc_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! Cloud depth (not cloud fraction
                                    ! weighted).
 z_cld_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    ! Cloud fraction weighted
                                    ! thickness of DSC cloud
 cld_factor_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                    ! Fraction of grid box potentially
                                    ! giving evaporative entrainment
 d_siems(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! Siems (1990) et al. cloud-top
                                    !  entr.t instab. parm
 d_siems_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                    ! Siems (1990) et al. cloud-top
                                    ! entr.t instab. parm for DSC
                                    ! layer
 z_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),  &
                                    ! Z_TOP(*,K) is the height of
                                    ! level k+1/2 above the surface.
 w_grad(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,            &
        bl_levels),                                                     &
                                    ! gradient of w
 tv1_sd(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! Standard Deviation of level 1
                                    ! virtual temperature (K).
 f_net_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
 f_net_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels),                                                  &
                                    ! Net radiative fluxes relative
                                    ! to SML and DSC base
 df_net_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! Net radiative divergences over
 df_net_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),       &
                                    ! various layer depths
 cf_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! cloud fraction of SML
 cf_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! cloud fraction of DSC layer
 z_cf_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    ! cloud base height from cld sch
 z_ctop(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! cloud top height
 dqw_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! QW change across SML disc inv
 dtl_sml(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! TL change across SML disc inv
 dqw_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! QW change across DSC disc inv
 dtl_dsc(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! TL change across DSC disc inv
 tls_inv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! TL subs increment across inv
 rhokh_surf_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                    ! SML surf-driven entr. KH
 rhokh_top_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),    &
                                    ! SML top-driven entr. KH
 rhokh_dsct_ent(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),   &
                                    ! DSC top-driven entr. KH
 zdsc_base(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    ! Height of base of K_top in DSC
 we_parm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! param.d entrainment rates (m/s)
 we_dsc_parm(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                    ! for surf and DSC layers
 zh_np1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),           &
                                    ! estimate of ZH at end of timeste
 zhsc_np1(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! estimate of ZHSC at end of times

REAL :: TmpScm2d(ScmRowLen,ScmRow)  ! For SCM diagnostics

INTEGER ::                                                              &
 k_cloud_top(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),      &
                                    ! Level number of top of b.l.
                                    ! cloud.
 k_cloud_dsct(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),     &
                                    ! Level number of top of dec.
                                    ! cloud.
 ntml_save(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    ! Copy of NTML
 ntml_prev(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),        &
                                    ! NTML from previous timestep
 k_plume(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! Start grid-level for
                                    ! surface-driven plume
 k_cbase(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! grid-level above cloud-base

INTEGER ::                                                              &
 w_nonmono(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,         &
           bl_levels)             ! 0/1 flag for w being non-monotonic

! NEC vectorization
INTEGER ::                                                              &
 k_level(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),          &
                                    ! array to store level selection
 k_cff(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                    ! level counter for CFF

!  (b) Scalars.

REAL ::                                                                 &
 virt_factor,                                                           &
                       ! Temporary in calculation of buoyancy
                       ! parameters.
 dtldz,                                                                 &
                       ! Vertical gradient of TL in a well-mixed
                       ! layer.
 dqw,                                                                   &
                       ! Total water content change across layer
                       ! interface.
 dtl,                                                                   &
                       ! Liquid/ice static energy change across
                       ! layer interface.
 dtl_ga,                                                                &
                       ! As DTL but inc gradient adjustment
 dqw_ga,                                                                &
                       ! As DQW but inc gradient adjustment
 dqcl,                                                                  &
                       ! Cloud liquid water change across layer
                       ! interface.
 dqcf,                                                                  &
                       ! Cloud frozen water change across layer
                       ! interface.
 q_vap_parc,                                                            &
                       ! Vapour content of parcel
 q_liq_parc,                                                            &
                       ! Condensed water content of parcel
 q_liq_env,                                                             &
                       ! Condensed water content of environment
 t_parc,                                                                &
                       ! Temperature of parcel
 t_dens_parc,                                                           &
                       ! Density potential temperature of parcel
 t_dens_env,                                                            &
                       ! Density potential temperature of
                       ! environment
 denv_bydz,                                                             &
                       ! Gradient of density potential
                       ! temperature in environment
 dpar_bydz,                                                             &
                       ! Gradient of density potential
                       ! temperature of parcel
 rho_dz,                                                                &
                       ! rho*dz
 svl_lapse,                                                             &
                    ! Lapse rate of SVL above inversion (K/m)
 svl_lapse_base,                                                        &
                    ! Lapse rate of SVL above inversion (K/m)
 dsvl_top,                                                              &
                 ! s_VL jump across inversion grid layer (K)
 tothf_inv,                                                             &
                 ! total heat flux at discontinous inversion height
 df_inv,                                                                &
                 ! temporary in rad divergence calculation
 dz_disc_min,                                                           &
                 ! smallest allowed DZ_DISC
 tl_k,                                                                  &
                 ! TL (full static energy) on level K
 tl_kp1,                                                                &
                 ! TL (full static energy) on level K+1
 tl_kp2,                                                                &
                 ! TL (full static energy) on level K+2
 db_disc,                                                               &
                 ! Temporary disc inversion buoyancy jump
 w_s_ent,                                                               &
                 ! numerical (subsidence) entrainment rate
 w_ls,                                                                  &
                 ! large-scale (subs) velocity
 w_ls_dsc,                                                              &
                 ! large-scale (subs) velocity
 dz_disc,                                                               &
                 ! height of ZH below Z_uv(NTML+2)
 z_surf,                                                                &
                 ! approx height of top of surface layer
 quad_a,                                                                &
                 ! term `a' in quadratic solver for DZ_DISC
 quad_bm,                                                               &
                 ! term `-b'in quadratic solver for DZ_DISC
 quad_c,                                                                &
                 ! term `c' in quadratic solver for DZ_DISC
 w_m,                                                                   &
                 ! scaling velocity for standard deviations
 w_s_cubed,                                                             &
                 ! convective velocity scale
 z_cbase,                                                               &
                 ! cloud base height (m)
 zdsc_cbase,                                                            &
                 ! DSC cloud base height (m)
 z_int,                                                                 &
                 ! depth of wb integral
 z_int_top,                                                             &
                 ! top of wb-integration depth
 cf_for_wb,                                                             &
                 ! CF for use in wb calculation for decoupling
 dfsw_top,                                                              &
                 ! SW radiative flux change assoc with cloud-top
 w_curv,                                                                &
                 ! curvature of w
 w_curv_nm,                                                             &
 w_del_nm        ! terms used to find where w is non-monotonic

INTEGER ::                                                              &
 i,                                                                     &
             ! Loop counter (horizontal field index).
 j,                                                                     &
             ! Offset counter in certain I loops.
 k,                                                                     &
             ! Loop counter (vertical level index).
 kl,                                                                    &
             ! K
 km1,                                                                   &
             ! K-1
 kp1,                                                                   &
             ! K+1
 kp2,                                                                   &
             ! K+2
 mbl,                                                                   &
             ! Maximum number of model layers allowed in the rapidly
             ! mixing layer; set to BL_LEVELS-1.
 k_rad_smlt,                                                            &
                ! highest SML level for radiative
                !   divergence calculation
 k_rad_lim,                                                             &

 ient,iScm,jScm ! loop counters

!Variables for cache-blocking
INTEGER            :: jj          ! Block index

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!-----------------------------------------------------------------------
!  0.  Check that the scalars input to define the grid are consistent.
!      See comments to routine SF_EXCH for details.
!-----------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  Set MBL, "maximum number of boundary levels" for the purposes of
!  boundary layer height calculation.

mbl = bl_levels - 1

!-----------------------------------------------------------------------
! Calculate Z_TOP (top of levels) and NTML from previous timestep
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ntml_prev(i,j) = 1
  END DO
END DO
DO k = 1, bl_levels-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      z_top(i,j,k) = z_uv(i,j,k+1)
        !------------------------------------------------------------
        !find NTML from previous TS (for accurate gradient adjustment
        !of profiles - also note that NTML LE BL_LEVELS-1)
        !------------------------------------------------------------
      IF ( zh_prev(i,j)  >=  z_uv(i,j,k+1) ) ntml_prev(i,j)=k
    END DO
  END DO
END DO
k = bl_levels
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    z_top(i,j,k) = z_uv(i,j,k) + dzl(i,j,k)
  END DO
END DO
!-----------------------------------------------------------------------
! Calculate SVL: conserved variable used to test for well mixed layers
!-----------------------------------------------------------------------
DO k = 1, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      svl(i,j,k) = ( tl(i,j,k) + grcp * z_tq(i,j,k) )                   &
                 * ( 1.0 + c_virtual*qw(i,j,k) )
    END DO
  END DO
END DO

!No halos
IF ( l_new_qsat_bl ) THEN
  IF ( l_mr_physics ) THEN
    CALL qsat_mix_new(qs,t,p,pdims%i_end,pdims%j_end,bl_levels)
  ELSE
    CALL qsat_new(qs,t,p,pdims%i_end,pdims%j_end,bl_levels)
  END IF
ELSE
  DO k = 1, bl_levels
    ! DEPENDS ON: qsat_mix
    CALL qsat_mix(qs(1,1,k),t(1,1,k),p(1,1,k),pdims%i_end*pdims%j_end,    &
      l_mr_physics)
  END DO
END IF


!--------------------------------------------------------------------
!..Calculate subsidence increments:
!--------------------------------------------------------------------
DO k = 1, bl_levels
  km1 = MAX( 1, k-1 )
  kp1 = MIN( bl_levels, k+1 )
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      w_grad(i,j,k) = (w(i,j,k)-w(i,j,km1))*rdz(i,j,k)
      w_curv_nm = w(i,j,kp1)-2.0*w(i,j,k)+w(i,j,km1)
      w_del_nm  = w(i,j,kp1)-w(i,j,km1)
      w_nonmono(i,j,k) = 0
      IF ( ABS(w_curv_nm) > ABS(w_del_nm) ) w_nonmono(i,j,k) = 1
      tls_inc(i,j,k) = 0.0
      qls_inc(i,j,k) = 0.0
    END DO
  END DO
END DO

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED(bl_levels, pdims, w, tls_inc, rdz, tl,                    &
!$OMP& dzl,w_nonmono,w_grad,qls_inc, qw,bl_segment_size, cp, g)         &
!$OMP& PRIVATE(jj,k,j,i,kp1,w_curv)
DO jj = pdims%j_start, pdims%j_end, bl_segment_size
  DO k = 2, bl_levels-1
    DO j = jj, MIN((jj+bl_segment_size)-1,pdims%j_end)
      DO i = pdims%i_start, pdims%i_end

        IF ( w(i,j,k)  < - TINY(1.0) .AND.                              &
             w(i,j,k-1)< - TINY(1.0) ) THEN
            !-----------------------------------------------------------
            ! Only needed in subsidence regions
            ! Also don't attempt coupling with dynamics if w has
            ! significant vertical structure
            !-----------------------------------------------------------
          w_curv = (w_grad(i,j,k+1)-w_grad(i,j,k))/dzl(i,j,k)
          IF (ABS(w_curv) > 1.0e-6 .AND. w_nonmono(i,j,k) == 1) THEN
                 ! large curvature at a turning point
            tls_inc(i,j,k-1) = 0.0  ! need to make sure increments in
            qls_inc(i,j,k-1) = 0.0  ! level below are also set to zero
            w(i,j,k-1) = 0.0
            w(i,j,k)   = 0.0
            w(i,j,k+1) = 0.0
          ELSE

            kp1 = k+1
            tls_inc(i,j,k) = - w(i,j,k) *rdz(i,j,kp1)                   &
                          * ( tl(i,j,kp1) - tl(i,j,k) +                 &
                              g/( cp*rdz(i,j,kp1) ) )
            qls_inc(i,j,k) = - w(i,j,k) *rdz(i,j,kp1)                   &
                               * ( qw(i,j,kp1) - qw(i,j,k) )
          END IF  ! safe to calculate increments
        END IF

      END DO
    END DO
  END DO
END DO !jj
!$OMP END PARALLEL DO

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    unstable(i,j) = (fb_surf(i,j) >  0.0)
    k_plume(i,j)  = -1
  END DO
END DO
        !------------------------------------------------------------
        ! Find grid-level above top of surface layer, taken
        ! to be at a height, z_surf, given by:
        !       Z_SURF = 0.1*ZH_PREV
        ! Use ZH_prev since that will have determined the shape
        ! of the time-level n profiles.
        !------------------------------------------------------------
DO k = 1, bl_levels-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      IF ( unstable(i,j) ) THEN

        z_surf = 0.1 * zh_prev(i,j)
        IF ( z_tq(i,j,k) >= z_surf .AND. k_plume(i,j) == -1 ) THEN
               !reached z_surf
          k_plume(i,j)=k
        END IF
        IF ( svl(i,j,k+1) >= svl(i,j,k)                                 &
                .AND. k_plume(i,j) == -1 ) THEN
               !reached inversion
          k_plume(i,j)=k
        END IF

      END IF
    END DO
  END DO
END DO

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF (k_plume(i,j) == -1) k_plume(i,j)=1
  END DO
END DO
!-----------------------------------------------------------------------
! 0.3a IF NOT CUMULUS:
!      Look for decoupled cloudy mixed-layer above b.l. top at ZH
!      (starting from level 3 and below 2.5km):
!      find cloud-base above SML inversion, ie. above NTML+1,
!      then cloud-top (ie. CF < SC_CFTOL)
!      and finally check that cloud is well-mixed.
!-----------------------------------------------------------------------
!      Initialise variables

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    cloud_base(i,j) = .FALSE.
    dsc(i,j) = .FALSE.
    coupled(i,j) = .FALSE.
    zhsc(i,j)    = 0.0
    ntdsc(i,j)   = 0
  END DO
END DO

DO k = 3, mbl
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      !----------------------------------------------------------------------
      !..Find cloud-base (where cloud here means CF > SC_CFTOL)
      !----------------------------------------------------------------------
      IF ( .NOT. cumulus(i,j) .AND.                                     &
           z_tq(i,j,k) < zmaxb_for_dsc .AND.                            &
           k  >   ntml(i,j)+1 .AND. cf(i,j,k)  >   sc_cftol             &
                            .AND. .NOT. cloud_base(i,j)                 &
        !                                  not yet found cloud-base
                                      .AND. .NOT. dsc(i,j) ) THEN
        !                                  not yet found a Sc layer
        cloud_base(i,j) = .TRUE.
      END IF
      IF ( cloud_base(i,j) .AND. .NOT. dsc(i,j) .AND.                   &
        !                  found cloud-base but not yet reached cloud-top
                     cf(i,j,k+1) < sc_cftol .AND.                       &
                     z_tq(i,j,k) < zmaxt_for_dsc                        &
        !                  got to cloud-top below ZMAXT_FOR_DSC
                   ) THEN
        cloud_base(i,j) = .FALSE.         ! reset CLOUD_BASE
          !-----------------------------------------------------------
          ! Look to see if at least top of cloud is well mixed:
          ! test SVL-gradient for top 2 pairs of levels, in case
          ! cloud top extends into the inversion.
          ! Parcel descent in Section 4.0 below will determine depth
          ! of mixed layer.
          !----------------------------------------------------------
        IF ( (svl(i,j,k)-svl(i,j,k-1))                                  &
                     /(z_tq(i,j,k)-z_tq(i,j,k-1))                       &
                                            <   max_svl_grad ) THEN
          dsc(i,j) = .TRUE.
          ntdsc(i,j) = k
          zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
        ELSE IF ( (svl(i,j,k-1)-svl(i,j,k-2))                           &
                      /(z_tq(i,j,k-1)-z_tq(i,j,k-2))                    &
                                          <   max_svl_grad ) THEN
            !---------------------------------------------------------
            ! Well-mixed layer with top at k-1 or k.  Check whether
            ! there is a buoyancy inversion between levels k-1 and k
            ! in a manner similar to the surface-driven plume: compare
            ! the buoyancy gradient between levels K-1 and K for an
            ! undiluted parcel and the environment
            !---------------------------------------------------------
          sl_plume(i,j) = tl(i,j,k-1) + grcp * z_tq(i,j,k-1)
          qw_plume(i,j) = qw(i,j,k-1)
          ! -------------------------------------------------------------------
          ! calculate parcel water by linearising qsat about the environmental
          ! temperature.
          ! -------------------------------------------------------------------
          IF (t(i,j,k) >  tm) THEN
            q_liq_parc = MAX( 0.0, ( qw_plume(i,j) - qs(i,j,k) -        &
              dqsdt(i,j,k)*                                             &
              ( sl_plume(i,j)-grcp*z_tq(i,j,k)-t(i,j,k) )               &
                                   ) *a_qs(i,j,k) )
            q_liq_env = MAX( 0.0, ( qw(i,j,k) - qs(i,j,k)               &
                        -dqsdt(i,j,k)*( tl(i,j,k) - t(i,j,k) )          &
                                   ) *a_qs(i,j,k) )
            ! add on the difference in the environment's ql as calculated
            ! by the partial condensation scheme (using some RH_CRIT value)
            ! and what it would be if RH_CRIT=1. This then imitates partial
            ! condensation in the parcel.
            q_liq_parc = q_liq_parc + qcl(i,j,k) + qcf(i,j,k)           &
                           - q_liq_env
            t_parc = sl_plume(i,j) - grcp * z_tq(i,j,k) +               &
                             lcrcp*q_liq_parc
          ELSE
            q_liq_parc = MAX( 0.0, ( qw_plume(i,j) - qs(i,j,k) -        &
              dqsdt(i,j,k)*                                             &
                ( sl_plume(i,j)-grcp*z_tq(i,j,k)-t(i,j,k) )             &
                                   ) *a_qs(i,j,k) )
            q_liq_env = MAX( 0.0, ( qw(i,j,k) - qs(i,j,k)               &
               -dqsdt(i,j,k)*( tl(i,j,k) - t(i,j,k) )                   &
                                   ) *a_qs(i,j,k) )
            ! add on difference in environment's ql between RH_CRIT
            ! and RH_CRIT=1
            q_liq_parc = q_liq_parc + qcl(i,j,k) + qcf(i,j,k)           &
                           - q_liq_env
            t_parc = sl_plume(i,j) - grcp * z_tq(i,j,k) +               &
                             lsrcp*q_liq_parc
          END IF
          q_vap_parc=qw_plume(i,j)-q_liq_parc

          t_dens_parc=t_parc*(1.0+c_virtual*q_vap_parc-q_liq_parc)
          t_dens_env=t(i,j,k)*                                          &
                     (1.0+c_virtual*q(i,j,k)-qcl(i,j,k)-qcf(i,j,k))
          ! find vertical gradients in parcel and environment SVL (using values
          ! from level below (K-1))
          env_svl_km1(i,j) = t(i,j,k-1) * ( 1.0+c_virtual*q(i,j,k-1)    &
               -qcl(i,j,k-1)-qcf(i,j,k-1) ) + grcp*z_tq(i,j,k-1)
          dpar_bydz=(t_dens_parc+grcp*z_tq(i,j,k)-                      &
                      env_svl_km1(i,j)) /                               &
                  (z_tq(i,j,k)-z_tq(i,j,k-1))
          denv_bydz=(t_dens_env+grcp*z_tq(i,j,k)-                       &
                      env_svl_km1(i,j))/                                &
                  (z_tq(i,j,k)-z_tq(i,j,k-1))

          IF ( denv_bydz >  1.25*dpar_bydz ) THEN
              ! there is an inversion between levels K-1 and K
            IF ( k  >=  ntml(i,j)+3 ) THEN
                ! if NTDSC EQ NTML+1 then assume we're looking
                ! at the same inversion and so don't set DSC
              ntdsc(i,j) = k-1
              zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
              dsc(i,j) = .TRUE.
            END IF
          ELSE
              ! no inversion between levels K-1 and K, assume there
              ! is an inversion between K and K+1 because of CF change
            ntdsc(i,j) = k
            zhsc(i,j)  = z_uv(i,j,ntdsc(i,j)+1)
            dsc(i,j) = .TRUE.
          END IF
        END IF
      END IF
    END DO
  END DO
END DO
!-----------------------------------------------------------------------
! 0.4a If the layer to ZHPAR is a cumulus layer capped by cloud and
!      an inversion, declare this layer a decoupled cloud layer and
!      set ZHSC and NTDSC accordingly.
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( (l_param_conv .AND.                                            &
            l_shallow(i,j) .AND. ntpar(i,j)  <   bl_levels )            &
              ! shallow cumulus layer within BL_LEVELS
         .OR. (.NOT. l_param_conv .AND.                                 &
              cumulus(i,j) .AND. ntpar(i,j)  <   bl_levels ) ) THEN
              ! cumulus layer and inversion found
      IF ( cf(i,j,ntpar(i,j))  >   sc_cftol  .OR.                       &
           cf(i,j,ntpar(i,j)+1)  >   sc_cftol ) THEN
           ! cloudy
        dsc(i,j)  = .TRUE.
        zhsc(i,j) = zhpar(i,j)
        ntdsc(i,j)= ntpar(i,j)
      END IF
    END IF
  END DO
END DO
!-----------------------------------------------------------------------
! 0.4b Calculate the radiative flux changes across cloud top for the
!      stratocumulus layer and thence a first guess for the top-down
!      mixing depth of this layer, DSCDEPTH.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    k_cloud_dsct(i,j) = 0
    df_dsct_over_cp(i,j) = 0.0
  END DO
END DO

DO k = 1, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      df_over_cp(i,j,k) = - ( rad_hr(i,j,1,k) + rad_hr(i,j,2,k) )       &
                          * rho_mix_tq(i,j,k) * dzl(i,j,k)
      dflw_over_cp(i,j,k) = - rad_hr(i,j,1,k)                           &
                          * rho_mix_tq(i,j,k) * dzl(i,j,k)
      dfsw_over_cp(i,j,k) = - rad_hr(i,j,2,k)                           &
                          * rho_mix_tq(i,j,k) * dzl(i,j,k)
        !-------------------------------------------------------------
        ! Find the layer with the greatest LW radiative flux jump and
        ! assume that this is the top DSC level.  Limit the
        ! search to above the SML.
        !-------------------------------------------------------------
      k_rad_lim = ntml(i,j)+2

      IF ( dsc(i,j) .AND. k  >=  k_rad_lim .AND. k  <=  ntdsc(i,j)+2    &
          .AND. dflw_over_cp(i,j,k)  >   df_dsct_over_cp(i,j) ) THEN
        k_cloud_dsct(i,j) = k
        df_dsct_over_cp(i,j) = dflw_over_cp(i,j,k)
      END IF

    END DO
  END DO
END DO
    !-----------------------------------------------------------------
    !  In case the LW radiative divergence is spread over 2 or 3 levs,
    !  add on any cooling either side of the level of maximum cooling.
    !-----------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( dsc(i,j) .AND. k_cloud_dsct(i,j)  >   0 ) THEN
      k=k_cloud_dsct(i,j)
      dfsw_top = dfsw_over_cp(i,j,k)
      IF ( k >  1 ) THEN
        df_dsct_over_cp(i,j) = df_dsct_over_cp(i,j)                     &
                                + MAX(0.0, dflw_over_cp(i,j,k-1) )
        dfsw_top = dfsw_top + MIN(0.0, dfsw_over_cp(i,j,k-1) )
      END IF
      IF ( k <  bl_levels ) THEN
        df_dsct_over_cp(i,j) = df_dsct_over_cp(i,j)                     &
                                + MAX(0.0, dflw_over_cp(i,j,k+1) )
        dfsw_top = dfsw_top + MIN(0.0, dfsw_over_cp(i,j,k+1) )
      END IF
          !-----------------------------------------------------------
          ! Combine SW and LW cloud-top divergences into a net
          ! divergence by estimating SW flux divergence at a given
          ! LW divergence = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
          ! Empirically (from LEM data) a reasonable fit is found
          ! with A small and (1-exp{-A*kappa_sw/kappa_lw}) = 0.35
          !-----------------------------------------------------------
      df_dsct_over_cp(i,j) = MAX( 0.0,                                  &
                df_dsct_over_cp(i,j) + 0.35 * dfsw_top )
    END IF
  END DO
END DO
!-----------------------------------------------------------------------
!     Set NBDSC, the bottom level of the DSC layer.
!     Note that this will only be used to give an estimate of the layer
!     depth, DSCDEPTH, used to calculate the entrainment
!     rate (the dependence is only weak), and that a more accurate
!     algorithm is subsequently used to determine the depth over which
!     the top-down mixing profiles will be applied.  If DSC is FALSE,
!     DSCDEPTH = 0.  The plume descent here uses a radiative
!     perturbation to the cloud-layer SVL (use level NTDSC-1 in case
!     SVL is not yet well-mixed to NTDSC), based roughly
!     on a typical cloud-top residence time.  If the plume does not sink
!     and the cloud is decoupled from the surface (ie. above Alan
!     Grant's ZH), then it is assumed to be stable, ie. St rather than
!     Sc, and no mixing or entrainment is applied to it.
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    nbdsc(i,j) = ntdsc(i,j)+1
    IF (dsc(i,j)) THEN
      ! The depth of the radiatively-cooled layer tends to be less than O(50m)
      ! and so RAD_HR will be an underestimate of the cooling tendency there.
      ! Compensate by multiplying by DZL/50. (~4)
      ! Recall that DF_OVER_CP(I,j,K) = RAD_HR * RHO_MIX_TQ * DZL
      ! Thus use cloud-top radiative forcing as follows:

      k = ntdsc(i,j)
      rho_dz = rho_mix_tq(i,j,k) * dzl(i,j,k)
      svl_plume(i,j)=svl(i,j,k-1)                                       &
         - ct_resid * dzl(i,j,k)*df_dsct_over_cp(i,j) / ( 50.0*rho_dz )

    ELSE
      svl_plume(i,j)=0.0
    END IF
  END DO
END DO

DO k = bl_levels-1, 1, -1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( k <  ntdsc(i,j) .AND. svl_plume(i,j)  <   svl(i,j,k) ) THEN
        nbdsc(i,j) = k+1     ! marks lowest level within ML
      END IF
    END DO
  END DO
END DO
!----------------------------------------------------------------------
!  0.4e Tidy up variables associated with decoupled layer
!       NOTE that NTDSC GE 3 if non-zero
!----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ! Note that ZHSC-Z_UV(NTML+2) may = 0, so this test comes first!
    IF (cumulus(i,j) .AND. dsc(i,j))                                    &
                     nbdsc(i,j) = MAX( nbdsc(i,j), ntml(i,j)+2 )
    IF ( ntdsc(i,j)  >=  1 ) THEN
      IF ( nbdsc(i,j) <   ntdsc(i,j)+1 ) THEN
        dscdepth(i,j) =                                                 &
                    z_uv(i,j,ntdsc(i,j)+1) - z_uv(i,j,nbdsc(i,j))
      ELSE
        !----------------------------------------------------------
        ! Indicates a layer of zero depth
        !----------------------------------------------------------
        IF (ntdsc(i,j) == ntpar(i,j)) THEN
          !----------------------------------------------------------
          ! Indicates a Sc layer at the top of Cu: force mixing
          ! over single layer.
          !----------------------------------------------------------
          dscdepth(i,j) = dzl(i,j,ntdsc(i,j))
        ELSE
          dsc(i,j)=.FALSE.
          ntdsc(i,j)=0
          zhsc(i,j)=0.0
          df_dsct_over_cp(i,j) = 0.0
          dscdepth(i,j) = 0.0
        END IF
      END IF
    ELSE  ! ntdsc eq 0, just to make sure
      dscdepth(i,j)=0.0
      dsc(i,j)=.FALSE.
      zhsc(i,j)=0.0
      df_dsct_over_cp(i,j) = 0.0
    END IF
  END DO
END DO
!----------------------------------------------------------------------
!..If decoupled cloud-layer found test to see if it is, in fact,
!..only weakly decoupled from the surface mixed-layer:
!..if SVL difference between NTML and NTDSC is less than SVL_COUP Kelvin
!..then assume there is still some coupling.  This will mean that
!..the surface-driven entrainment term will be applied at ZHSC, no
!..subgrid inversion or entrainment will be calculated for ZH and
!..ZHSC will be the length scale used in the entrainment inputs.
!..Note that for CUMULUS surface-driven entrainment will be done
!..by the convection scheme.
!----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    coupled(i,j) = .FALSE.
    IF ( dsc(i,j) .AND. .NOT. cumulus(i,j) ) THEN
        !------------------------------------------------------------
        ! Note this IF test structure is required because if DSC is
        ! false then NTDSC = 0 and cannot be used to index SVL.
        !------------------------------------------------------------
      IF ( svl(i,j,ntdsc(i,j)) - svl(i,j,ntml(i,j))  <   svl_coup )     &
         coupled(i,j) = .TRUE.
    END IF
  END DO
END DO
!-----------------------------------------------------------------------
! Assuming a discontinuous inversion structure.
!    - to this point in the code, ZH and ZHSC mark the half-level at
!      the base of the inversion
!    - now they will be interpolated into the level above assuming
!      SVL(NTML+1) is a volume average over a subgrid discontinuous
!      inversion structure
!    - discontinuous jumps of TL and QW (and thence buoyancy) can be
!      calculated and used to determine the entrainment rate
!    - parametrized grid-level fluxes at NTML,NTDSC can then be made
!      consistent with this assumed inversion structure
!-----------------------------------------------------------------------
       ! If any `problems' are encountered with this interpolation of ZH
       ! (such as ZH diagnosed at or below Z_UV(NTML+1)), then NTML
       ! is lowered a level and ZH is set fractionally below what has
       ! become Z_UV(NTML+2).  This distance is such that for a net
       ! dZH/dt of 1.E-4 m/s, ZH will be diagnosed as spending at least
       ! half the timestep in level NTML+2, leaving the growth only
       ! marginally affected.  Conversely, it allows a subsiding
       ! inversion to fall more readily.

dz_disc_min = 0.5 * timestep * 1.0e-4
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    sml_disc_inv(i,j) = 0  ! initialise flags to indicate whether a
    dsc_disc_inv(i,j) = 0  ! discontinuous inversion is diagnosed
    !..First interpolate to find ZH

    k = ntml(i,j)
    !..by default, keep ZH at the half-level where it was diagnosed
    !..initially and use grid-level jumps

    dtl_sml(i,j) = tl(i,j,k+1) - tl(i,j,k) + grcp/rdz(i,j,k+1)
    dqw_sml(i,j) = qw(i,j,k+1) - qw(i,j,k)

    IF ( .NOT. cumulus(i,j) .AND. .NOT. coupled(i,j) .AND.              &
                  k  >   1 .AND. k  <=  bl_levels-2 ) THEN
      IF ( svl(i,j,k+2)  >   svl(i,j,k+1)                               &
                    .AND. svl(i,j,k+1)  >   svl(i,j,k) ) THEN

        tl_k = tl(i,j,k) + grcp * z_tq(i,j,k)
        tl_kp1 = tl(i,j,k+1) + grcp * z_tq(i,j,k+1)
        tl_kp2 = tl(i,j,k+2) + grcp * z_tq(i,j,k+2)

        IF ( k  <=  bl_levels-3 ) THEN
          ! need to test for K+1 to K+2 gradient in case profile is concave
          ! (would mess up the inversion diagnosis)
          svl_lapse = MAX(0.0,                                          &
            MIN( ( svl(i,j,k+2) - svl(i,j,k+1) ) * rdz(i,j,k+2),        &
                 ( svl(i,j,k+3) - svl(i,j,k+2) ) * rdz(i,j,k+3)         &
               )         )
        ELSE
          svl_lapse = 0.0
        END IF
        IF ( k  >=  k_plume(i,j)+2 ) THEN
              ! Use mean mixed layer gradient (if resolved) to allow
              ! for stablisation by gradient-adjustment
              ! Ignore level K in case inversion is dropping
          svl_lapse_base = ( svl(i,j,k-1)-svl(i,j,k_plume(i,j)) )/      &
                        (z_tq(i,j,k-1)-z_tq(i,j,k_plume(i,j)))
          svl_lapse_base = MAX( 0.0, svl_lapse_base )
        ELSE
          svl_lapse_base = 0.0
        END IF

        quad_a  = 0.5*( svl_lapse - svl_lapse_base )
        quad_bm = svl(i,j,k+2) - svl(i,j,k)                             &
            - svl_lapse * ( z_tq(i,j,k+2)-z_uv(i,j,k+2) )               &
            - svl_lapse_base * ( z_uv(i,j,k+1)-z_tq(i,j,k) +            &
                                                    dzl(i,j,k+1) )
        quad_c  = dzl(i,j,k+1)*( svl(i,j,k+1) - svl(i,j,k) -            &
            svl_lapse_base * (                                          &
              z_uv(i,j,k+1)-z_tq(i,j,k) + 0.5*dzl(i,j,k+1) ) )

        IF ( quad_bm  >   0.0 ) THEN
          IF ( quad_c  <=  0.0) THEN
                ! SVL extrapolated from K to K+1 is greater than
                ! the level K+1 value - inversion needs to rise so
                ! place it as high as possible
            dz_disc = dz_disc_min
          ELSE IF ( quad_bm*quad_bm  >=  4.0*quad_a*quad_c ) THEN
                ! solve equation for DZ_DISC...
            IF ( quad_a  /=  0.0 ) THEN
                  !   ...quadratic if QUAD_A NE 0
              dz_disc = ( quad_bm - SQRT( quad_bm*quad_bm               &
                                       - 4.0*quad_a*quad_c )            &
                              ) / (2.0*quad_a)
            ELSE
                  !   ...linear if QUAD_A EQ 0
              dz_disc = quad_c / quad_bm
            END IF
          ELSE
            dz_disc = 99999.9  ! large dummy value
          END IF

          IF ( dz_disc  >   0.9 * dzl(i,j,k+1) ) THEN
              ! ZH diagnosed very close to or below Z_UV(K+1):
            IF ( svl(i,j,k)-svl(i,j,k-1)  >   0.0) THEN
                  ! top of ML stably stratified so lower NTML but
                  ! set ZH only fractionally (DZ_DISC_MIN)
                  ! below the top of the inversion level.
              ntml(i,j) = ntml(i,j) - 1
              k=ntml(i,j)
              tl_k = tl(i,j,k) + grcp * z_tq(i,j,k)
              tl_kp1 = tl(i,j,k+1) + grcp * z_tq(i,j,k+1)
              tl_kp2 = tl(i,j,k+2) + grcp * z_tq(i,j,k+2)
              dz_disc = dz_disc_min
            ELSE
                  ! top of ML well-mixed so don't lower the inversion
                  ! level but set ZH just (DZ_DISC_MIN) above the
                  ! half-level to allow the inversion to subside if
                  ! necessary.
              dz_disc = dzl(i,j,k+1) - dz_disc_min
            END IF
          END IF

        ELSE
          !.. ignoring lapse rates
          dsvl_top = svl(i,j,k+2) - svl(i,j,k)
          dz_disc = dzl(i,j,k+1) *                                      &
                          (svl(i,j,k+1)-svl(i,j,k)) / dsvl_top
        END IF

        zh(i,j) = z_uv(i,j,k+2) - dz_disc
        sml_disc_inv(i,j) = 1 ! set flag to indicate disc inv found

        !..Calculate discontinuous jumps of TL and QW:
        !..only accurate enough for large enough DZ_DISC

        IF (dz_disc/dzl(i,j,k+1)  >   0.1) THEN
          dtl_sml(i,j) = (tl_kp1-tl_k) * dzl(i,j,k+1) / dz_disc

          IF ( tl_kp2  >   tl_kp1 .AND.                                 &
                             tl_kp1  >   tl_k ) THEN
            dtl_sml(i,j) = MIN( tl_kp2-tl_k, dtl_sml(i,j) )
          ELSE IF ( tl_kp2  <   tl_kp1 .AND.                            &
                              tl_kp1  <   tl_k ) THEN
            dtl_sml(i,j) = MAX( tl_kp2-tl_k, dtl_sml(i,j) )
          ELSE  ! TL non-monotonic
            dtl_sml(i,j) = tl_kp1 - tl_k
          END IF

          dqw_sml(i,j) = (qw(i,j,k+1)-qw(i,j,k))                        &
                                       * dzl(i,j,k+1) / dz_disc

          IF ( qw(i,j,k+2)  >   qw(i,j,k+1) .AND.                       &
                 qw(i,j,k+1)  >   qw(i,j,k) ) THEN
            dqw_sml(i,j) = MIN( qw(i,j,k+2)-qw(i,j,k), dqw_sml(i,j))
          ELSE IF ( qw(i,j,k+2)  <   qw(i,j,k+1) .AND.                  &
                  qw(i,j,k+1)  <   qw(i,j,k) ) THEN
            dqw_sml(i,j) = MAX( qw(i,j,k+2)-qw(i,j,k), dqw_sml(i,j))
          ELSE  ! QW non-monotonic
            dqw_sml(i,j) = qw(i,j,k+1)-qw(i,j,k)
          END IF

        ELSE
          !DZ_DISC small, ZH close to level above: use double grid-level jumps
          dtl_sml(i,j) = tl_kp2 - tl_k
          dqw_sml(i,j) = qw(i,j,k+2) - qw(i,j,k)
        END IF

      END IF  ! SVL increasing
    END IF ! not cumulus and not at top of bl_levels
    !-----------------------------------------------------------------------
    !..Second interpolate to find ZHSC
    !-----------------------------------------------------------------------
    IF ( dsc(i,j) ) THEN
      k = ntdsc(i,j)
      !..by default, keep ZHSC at the half-level where it was diagnosed
      !..initially and use grid-level jumps
      dtl_dsc(i,j) = tl(i,j,k+1) - tl(i,j,k) + grcp/rdz(i,j,k+1)
      dqw_dsc(i,j) = qw(i,j,k+1) - qw(i,j,k)
      IF ( k  <=  bl_levels-2 ) THEN
        IF ( svl(i,j,k+2)  >   svl(i,j,k+1)                             &
                 .AND. svl(i,j,k+1)  >   svl(i,j,k) ) THEN

          tl_k = tl(i,j,k) + grcp * z_tq(i,j,k)
          tl_kp1 = tl(i,j,k+1) + grcp * z_tq(i,j,k+1)
          tl_kp2 = tl(i,j,k+2) + grcp * z_tq(i,j,k+2)

          IF ( k  <=  bl_levels-3 ) THEN
            ! need to test for K+1 to K+2 gradient in case profile is concave
            ! (would mess up the inversion diagnosis)
            svl_lapse = MAX( 0.0,                                       &
              MIN( ( svl(i,j,k+2) - svl(i,j,k+1) ) * rdz(i,j,k+2),      &
                   ( svl(i,j,k+3) - svl(i,j,k+2) ) * rdz(i,j,k+3)       &
                 )         )
          ELSE
            svl_lapse = 0.0
          END IF
          IF ( k  >=  nbdsc(i,j)+2 ) THEN
              ! Use mean mixed layer gradient (if resolved) to allow
              ! for stablisation by gradient-adjustment
              ! Ignore level K in case inversion is dropping
            svl_lapse_base = ( svl(i,j,k-1)-svl(i,j,nbdsc(i,j)) )/      &
                          (z_tq(i,j,k-1)-z_tq(i,j,nbdsc(i,j)))
            svl_lapse_base = MAX( 0.0, svl_lapse_base )
          ELSE
            svl_lapse_base = 0.0
          END IF

          quad_a  = 0.5*( svl_lapse - svl_lapse_base )
          quad_bm = svl(i,j,k+2) - svl(i,j,k)                           &
               - svl_lapse * ( z_tq(i,j,k+2)-z_uv(i,j,k+2) )            &
               - svl_lapse_base * ( z_uv(i,j,k+1)-z_tq(i,j,k) +         &
                                                      dzl(i,j,k+1) )
          quad_c  = dzl(i,j,k+1)*( svl(i,j,k+1) - svl(i,j,k) -          &
               svl_lapse_base * (                                       &
                z_uv(i,j,k+1)-z_tq(i,j,k) + 0.5*dzl(i,j,k+1) ) )

          IF ( quad_bm  >   0.0 ) THEN
            IF ( quad_c  <=  0.0) THEN
                ! SVL extrapolated from K to K+1 is greater than
                ! the level K+1 value - inversion needs to rise
              dz_disc = dz_disc_min
            ELSE IF ( quad_bm*quad_bm  >=  4.0*quad_a*quad_c ) THEN
                ! solve equation for DZ_DISC...
              IF ( quad_a  /=  0.0 ) THEN
                  !   ...quadratic if QUAD_A NE 0
                dz_disc = ( quad_bm - SQRT( quad_bm*quad_bm             &
                                         - 4.0*quad_a*quad_c )          &
                                ) / (2.0*quad_a)
              ELSE
                  !   ...linear if QUAD_A EQ 0
                dz_disc = quad_c / quad_bm
              END IF
            ELSE
              dz_disc = 99999.9  ! large dummy value
            END IF

            IF ( dz_disc  >   0.9 * dzl(i,j,k+1) ) THEN
              IF ( ntdsc(i,j) == 2 ) THEN
                dz_disc = dzl(i,j,k+1)
              ELSE
                ! ZHSC diagnosed very close to or below Z_UV(K+1):
                IF ( svl(i,j,k)-svl(i,j,k-1)  >   0.0) THEN
                  ! top of ML stably stratified so lower NTDSC but
                  ! set ZHSC only fractionally (DZ_DISC_MIN)
                  ! below the top of the inversion level.
                  ntdsc(i,j) = ntdsc(i,j) - 1
                  k=ntdsc(i,j)
                  tl_k = tl(i,j,k) + grcp * z_tq(i,j,k)
                  tl_kp1 = tl(i,j,k+1) + grcp * z_tq(i,j,k+1)
                  tl_kp2 = tl(i,j,k+2) + grcp * z_tq(i,j,k+2)
                  dz_disc = dz_disc_min
                  dscdepth(i,j) = dscdepth(i,j) - dzl(i,j,k+1)
                  ! Note that all but DZ_DISC_MIN of this layer will
                  ! be added back on to DSCDEPTH a few lines below
                ELSE
                  ! top of ML well-mixed so don't lower the inversion
                  ! level but set ZHSC just (DZ_DISC_MIN) above the
                  ! half-level to allow the inversion to subside if
                  ! necessary.
                  dz_disc = dzl(i,j,k+1) - dz_disc_min
                END IF
              END IF
            END IF

          ELSE  ! QUAD_BM le 0
            !.. ignoring lapse rates
            dsvl_top = svl(i,j,k+2) - svl(i,j,k)
            dz_disc = dzl(i,j,k+1) *                                    &
                            (svl(i,j,k+1)-svl(i,j,k)) / dsvl_top
          END IF

          zhsc(i,j) = z_uv(i,j,k+2) - dz_disc
          dscdepth(i,j) = dscdepth(i,j) + zhsc(i,j) - z_uv(i,j,k+1)
          dsc_disc_inv(i,j) = 1  ! set flag to indicate disc inv found

          !..Calculate discontinuous jumps of TL and QW:
          !..only accurate enough for large enough DZ_DISC
          IF (dz_disc/dzl(i,j,k+1)  >   0.1) THEN

            dtl_dsc(i,j) = (tl_kp1-tl_k) * dzl(i,j,k+1) / dz_disc

            IF ( tl_kp2  >   tl_kp1 .AND.                               &
                               tl_kp1  >   tl_k ) THEN
              dtl_dsc(i,j) = MIN( tl_kp2-tl_k, dtl_dsc(i,j) )
            ELSE IF ( tl_kp2  <   tl_kp1 .AND.                          &
                                tl_kp1  <   tl_k ) THEN
              dtl_dsc(i,j) = MAX( tl_kp2-tl_k, dtl_dsc(i,j) )
            ELSE  ! TL non-monotonic
              dtl_dsc(i,j) = tl_kp1 - tl_k
            END IF

            dqw_dsc(i,j) = (qw(i,j,k+1)-qw(i,j,k))                      &
                                         * dzl(i,j,k+1) / dz_disc

            IF ( qw(i,j,k+2)  >   qw(i,j,k+1) .AND.                     &
                   qw(i,j,k+1)  >   qw(i,j,k) ) THEN
              dqw_dsc(i,j) = MIN( qw(i,j,k+2)-qw(i,j,k), dqw_dsc(i,j))
            ELSE IF ( qw(i,j,k+2)  <   qw(i,j,k+1) .AND.                &
                    qw(i,j,k+1)  <   qw(i,j,k) ) THEN
              dqw_dsc(i,j) = MAX( qw(i,j,k+2)-qw(i,j,k), dqw_dsc(i,j))
            ELSE  ! QW non-monotonic
              dqw_dsc(i,j) = qw(i,j,k+1)-qw(i,j,k)
            END IF

          ELSE
            !..DZ_DISC small, ZHSC close to level above: use double jumps
            dtl_dsc(i,j) = tl_kp2 - tl_k
            dqw_dsc(i,j) = qw(i,j,k+2) - qw(i,j,k)
          END IF

        END IF  ! SVL increasing
      END IF ! test on K LT BL_LEVELS-2
    END IF ! test on DSC
  END DO
END DO
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    !-----------------------------------------------------------------------
    ! 0.5 Calculate the within-layer vertical gradients of cloud liquid
    !     and frozen water for the layer 1
    !-----------------------------------------------------------------------

    virt_factor = 1.0 + c_virtual*q(i,j,1) - qcl(i,j,1) - qcf(i,j,1)

      ! Here this is an estimate of the gradient adjustment applied
      ! the previous timestep (assumes T1_SD has not changed much,
      ! which in turn assumes RHOKH(1) has not)
    grad_t_adj(i,j) = MIN( max_t_grad ,                                 &
                         a_grad_adj * t1_sd(i,j) / zh_prev(i,j) )
    !        IF (T1_SD(I,j)  >   0.0) THEN
    !          GRAD_Q_ADJ(I,j) = (Q1_SD(I,j) / T1_SD(I,j)) * GRAD_T_ADJ(I,j)
    !        ELSE
    grad_q_adj(i,j) = 0.0
    !        END IF
    dtldz = -grcp + grad_t_adj(i,j)
    dqcldz(i,j,1) = -( dtldz*dqsdt(i,j,1) +                             &
                     g*qs(i,j,1)/(r*t(i,j,1)*virt_factor) )             &
                    / (1.0 + lcrcp*dqsdt(i,j,1))
    dqcfdz(i,j,1) = -( dtldz*dqsdt(i,j,1) +                             &
                     g*qs(i,j,1)/(r*t(i,j,1)*virt_factor) ) * fgf       &
                    / (1.0 + lsrcp*dqsdt(i,j,1))

    !-----------------------------------------------------------------------
    ! 0.6 Calculate the cloud liquid and frozen water contents at the
    !     top and bottom of layer 1
    !-----------------------------------------------------------------------

    !MHM limit calculation to greater than a small cloud fraction
    IF ( qcl(i,j,1) + qcf(i,j,1)  >   0.0                               &
         .AND. cf(i,j,1)  >   1.0e-3 ) THEN
      cfl(i,j,1) = cf(i,j,1) * qcl(i,j,1)/(qcl(i,j,1)+qcf(i,j,1))
      cff(i,j,1) = cf(i,j,1) * qcf(i,j,1)/(qcl(i,j,1)+qcf(i,j,1))
    ELSE
      cfl(i,j,1) = 0.0
      cff(i,j,1) = 0.0
    END IF

    IF (cfl(i,j,1)  >   0.0) THEN
      qcl_ic_top(i,j,1) = qcl(i,j,1) / cfl(i,j,1) +                     &
                         0.5*dzl(i,j,1)*dqcldz(i,j,1)
    ELSE
      qcl_ic_top(i,j,1) = 0.0
    END IF

    IF (cff(i,j,1)  >   0.0) THEN
      qcf_ic_top(i,j,1) = qcf(i,j,1) / cff(i,j,1) +                     &
                         0.5*dzl(i,j,1)*dqcfdz(i,j,1)
    ELSE
      qcf_ic_top(i,j,1) = 0.0
    END IF

    qcl_ic_bot(i,j,1) = 0.0
    qcf_ic_bot(i,j,1) = 0.0

  END DO
END DO

!-----------------------------------------------------------------------
! 1.  First loop round boundary layer levels.
!-----------------------------------------------------------------------
!Parameters used: fgf
!$OMP  PARALLEL DEFAULT(NONE)                                           &
!$OMP& SHARED(bl_levels, pdims, ntml_prev, grad_t_adj,                  &
!$OMP& q, qcl, qcf, dqsdt, dqcfdz, t, cf, cfl, cff, dqcldz, dzl,        &
!$OMP& qcl_ic_top, qcl_ic_bot, qcf_ic_top, qcf_ic_bot, qs, qw, tl, rdz, &
!$OMP& btm, bqm, bqm_cld, btm_cld, a_dqsdtm, a_qsm, db,                 &
!$OMP& db_cld, chi_s, zc, zc_dsc, r, c_virtual, g,                      &
!$OMP& etar, grcp, lcrcp, lsrcp ) &
!$OMP& PRIVATE(k,j,i,dtldz,virt_factor,km1,dqw,dtl,dqcl,dqcf)

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels

  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      !-----------------------------------------------------------------
      ! 1.4 Calculate the within-layer vertical gradients of cloud liquid
      !     and frozen water for the current layer
      !-----------------------------------------------------------------

      IF (k  <=  ntml_prev(i,j)) THEN
        dtldz = -grcp + grad_t_adj(i,j)
      ELSE
        dtldz = -grcp
      END IF

      ! RNBS correction  28/11/97
      virt_factor = 1.0 + c_virtual*q(i,j,k) - qcl(i,j,k) -             &
                          qcf(i,j,k)

      dqcldz(i,j,k) = -( dtldz*dqsdt(i,j,k)                             &
                     + g*qs(i,j,k)/(r*t(i,j,k)*virt_factor) )           &
                      / ( 1.0 + lcrcp*dqsdt(i,j,k) )
      dqcfdz(i,j,k) = -( dtldz*dqsdt(i,j,k)                             &
                     + g*qs(i,j,k)/(r*t(i,j,k)*virt_factor) ) * fgf     &
                      / ( 1.0 + lsrcp*dqsdt(i,j,k) )

      !-----------------------------------------------------------------------
      ! 1.5 Calculate the cloud liquid and frozen water contents at the
      !     top and bottom of the current layer
      !-----------------------------------------------------------------------

      !MHM limit calculation to greater than a small cloud fraction
      IF ( qcl(i,j,k) + qcf(i,j,k)  >   0.0                             &
           .AND. cf(i,j,k)  >   1.0e-3 ) THEN
        cfl(i,j,k) = cf(i,j,k) * qcl(i,j,k)                             &
                                  /( qcl(i,j,k) + qcf(i,j,k) )
        cff(i,j,k) = cf(i,j,k) * qcf(i,j,k)                             &
                                  /( qcl(i,j,k) + qcf(i,j,k) )
      ELSE
        cfl(i,j,k) = 0.0
        cff(i,j,k) = 0.0
      END IF

      IF (cfl(i,j,k)  >   0.0) THEN
        qcl_ic_top(i,j,k) = qcl(i,j,k) / cfl(i,j,k) +                   &
                           0.5*dzl(i,j,k)*dqcldz(i,j,k)
        qcl_ic_bot(i,j,k) = MAX( 0.0 , qcl(i,j,k) / cfl(i,j,k) -        &
                                      0.5*dzl(i,j,k)*dqcldz(i,j,k) )
      ELSE
        qcl_ic_top(i,j,k) = 0.0
        qcl_ic_bot(i,j,k) = 0.0
      END IF

      IF (cff(i,j,k)  >   0.0) THEN
        qcf_ic_top(i,j,k) = qcf(i,j,k) / cff(i,j,k) +                   &
                           0.5*dzl(i,j,k)*dqcfdz(i,j,k)
        qcf_ic_bot(i,j,k) = MAX( 0.0 , qcf(i,j,k) / cff(i,j,k) -        &
                                      0.5*dzl(i,j,k)*dqcfdz(i,j,k) )
      ELSE
        qcf_ic_top(i,j,k) = 0.0
        qcf_ic_bot(i,j,k) = 0.0
      END IF
    END DO
  END DO
END DO
!$OMP END DO

!-----------------------------------------------------------------------
! 2.  Second loop round boundary layer levels.
!-----------------------------------------------------------------------

!$OMP DO SCHEDULE(STATIC)
DO k = 2, bl_levels
  km1 = k-1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end

      !-----------------------------------------------------------------
      ! 2.1 Calculate the jumps of QW and TL across the layer interface
      !     at level k-1/2.
      !-----------------------------------------------------------------

      dqw = qw(i,j,k) - qw(i,j,km1)
      !-----------------------------------------------------------------
      ! TL_K_BOT   = TL(I,j,K)   - 0.5*DZL(I,j,K)  *(-GRCP)
      ! TL_KM1_TOP = TL(I,j,KM1) + 0.5*DZL(I,j,KM1)*(-GRCP)
      ! DTL = TL_K_BOT - TL_KM1_TOP   so therefore
      !-----------------------------------------------------------------
      dtl = tl(i,j,k) - tl(i,j,km1) + grcp/rdz(i,j,k)

      dqcl = cfl(i,j,k)*qcl_ic_bot(i,j,k) -                             &
               cfl(i,j,km1)*qcl_ic_top(i,j,km1)
      dqcf = cff(i,j,k)*qcf_ic_bot(i,j,k) -                             &
               cff(i,j,km1)*qcf_ic_top(i,j,km1)

      !--------------------------------------------------------------------
      ! 2.3 Calculate the buoyancy jumps across the interface between layers
      !     k and k-1
      !--------------------------------------------------------------------

      db(i,j,k) = g * ( btm(i,j,km1)*dtl + bqm(i,j,km1)*dqw +           &
                 (lcrcp*btm(i,j,km1) - etar*bqm(i,j,km1)) * dqcl +      &
                  (lsrcp*btm(i,j,km1) - etar*bqm(i,j,km1)) * dqcf )

      db_cld(i,j,k) =g * ( btm_cld(i,j,km1)*dtl                         &
           + bqm_cld(i,j,km1)*dqw )

      !ajm   + changed to - in chi_s calculation at v2p9
      chi_s(i,j,k) = -qcl_ic_top(i,j,km1) /                             &
                    (a_qsm(i,j,km1)*dqw - a_dqsdtm(i,j,km1)*dtl)
    END DO
  END DO
END DO
!$OMP END DO

!$OMP END PARALLEL

!----------------------------------------------------------------------
!..For levels NTML+1 and NTDSC+1, replace grid-level jumps calculated
!..above with jumps calculated assuming a discontinuous subgrid
!..inversion structure.
!----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    IF ( sml_disc_inv(i,j)  ==  1 ) THEN
      k = ntml(i,j)+1

      !..here _top means top of mixed layer, _bot means bottom of free atmos
      IF (cfl(i,j,k-1)  >   0.0) THEN
        qcl_ic_top(i,j,k) = qcl(i,j,k-1)/cfl(i,j,k-1) +                 &
                   ( zh(i,j)-z_tq(i,j,k-1) )*dqcldz(i,j,k-1)
      ELSE
        qcl_ic_top(i,j,k) = 0.0
      END IF

      IF (cfl(i,j,k+1)  >   0.0) THEN
        qcl_ic_bot(i,j,k) = MAX( 0.0 , qcl(i,j,k+1)/cfl(i,j,k+1)        &
            - ( z_tq(i,j,k+1)-zh(i,j) )*dqcldz(i,j,k+1) )
      ELSE
        qcl_ic_bot(i,j,k) = 0.0
      END IF

      IF (cff(i,j,k-1)  >   0.0) THEN
        qcf_ic_top(i,j,k) = qcf(i,j,k-1)/cff(i,j,k-1) +                 &
                  ( zh(i,j)-z_tq(i,j,k-1) )*dqcfdz(i,j,k-1)
      ELSE
        qcf_ic_top(i,j,k) = 0.0
      END IF

      IF (cff(i,j,k+1)  >   0.0) THEN
        qcf_ic_bot(i,j,k) = MAX( 0.0 , qcf(i,j,k+1)/cff(i,j,k+1)        &
            - ( z_tq(i,j,k+1)-zh(i,j) )*dqcfdz(i,j,k+1) )
      ELSE
        qcf_ic_bot(i,j,k) = 0.0
      END IF

      dqcl = cfl(i,j,k+1)*qcl_ic_bot(i,j,k) -                           &
             cfl(i,j,k-1)*qcl_ic_top(i,j,k)
      dqcf = cff(i,j,k+1)*qcf_ic_bot(i,j,k) -                           &
             cff(i,j,k-1)*qcf_ic_top(i,j,k)

      db_disc = g * ( btm(i,j,k-1)*dtl_sml(i,j)                         &
                      + bqm(i,j,k-1)*dqw_sml(i,j) +                     &
             (lcrcp*btm(i,j,k-1) - etar*bqm(i,j,k-1)) * dqcl +          &
             (lsrcp*btm(i,j,k-1) - etar*bqm(i,j,k-1)) * dqcf   )

      IF ( db_disc  <   0.0 ) THEN
            ! Diagnosed inversion statically unstable:
            ! use entrainment K (rather than fluxes) but also
            ! ensure DB>0 so that entrainment is non-zero
        sml_disc_inv(i,j) = 0
        zh(i,j) = z_uv(i,j,ntml(i,j)+1)
      ELSE
        db(i,j,k) = db_disc

        db_cld(i,j,k) = g * ( btm_cld(i,j,k-1)*dtl_sml(i,j) +           &
                            bqm_cld(i,j,k-1)*dqw_sml(i,j) )

        chi_s(i,j,k) = -qcl_ic_top(i,j,k) /                             &
            ( a_qsm(i,j,k-1)*dqw_sml(i,j)                               &
                  - a_dqsdtm(i,j,k-1)*dtl_sml(i,j) )
      END IF

    END IF  ! disc inversion diagnosed

    IF ( dsc_disc_inv(i,j)  ==  1 ) THEN
      k = ntdsc(i,j)+1

      !..here _top means top of mixed layer, _bot means bottom of free atmos
      IF (cfl(i,j,k-1)  >   0.0) THEN
        qcl_ic_top(i,j,k) = qcl(i,j,k-1)/cfl(i,j,k-1) +                 &
            ( zhsc(i,j)-z_tq(i,j,k-1) )*dqcldz(i,j,k-1)
      ELSE
        qcl_ic_top(i,j,k) = 0.0
      END IF

      IF (cfl(i,j,k+1)  >   0.0) THEN
        qcl_ic_bot(i,j,k) = MAX( 0.0 , qcl(i,j,k+1)/cfl(i,j,k+1)        &
             - ( z_tq(i,j,k+1)-zhsc(i,j) )*dqcldz(i,j,k+1) )
      ELSE
        qcl_ic_bot(i,j,k) = 0.0
      END IF

      IF (cff(i,j,k-1)  >   0.0) THEN
        qcf_ic_top(i,j,k) = qcf(i,j,k-1)/cff(i,j,k-1) +                 &
            ( zhsc(i,j)-z_tq(i,j,k-1) )*dqcfdz(i,j,k-1)
      ELSE
        qcf_ic_top(i,j,k) = 0.0
      END IF

      IF (cff(i,j,k+1)  >   0.0) THEN
        qcf_ic_bot(i,j,k) = MAX( 0.0 , qcf(i,j,k+1)/cff(i,j,k+1)        &
             - ( z_tq(i,j,k+1)-zhsc(i,j) )*dqcfdz(i,j,k+1) )
      ELSE
        qcf_ic_bot(i,j,k) = 0.0
      END IF

      dqcl = cfl(i,j,k+1)*qcl_ic_bot(i,j,k) -                           &
             cfl(i,j,k-1)*qcl_ic_top(i,j,k)
      dqcf = cff(i,j,k+1)*qcf_ic_bot(i,j,k) -                           &
             cff(i,j,k-1)*qcf_ic_top(i,j,k)

      db_disc = g * (                                                   &
         btm(i,j,k-1)*dtl_dsc(i,j) + bqm(i,j,k-1)*dqw_dsc(i,j) +        &
               (lcrcp*btm(i,j,k-1) - etar*bqm(i,j,k-1)) * dqcl +        &
               (lsrcp*btm(i,j,k-1) - etar*bqm(i,j,k-1)) * dqcf   )

      IF ( db_disc  <   0.0 ) THEN
            ! Diagnosed inversion statically unstable:
            ! use entrainment K (rather than fluxes) but also
            ! ensure DB>0 so that entrainment is non-zero
        dsc_disc_inv(i,j) = 0
        zhsc(i,j) = z_uv(i,j,ntdsc(i,j)+1)
      ELSE
        db(i,j,k) = db_disc

        db_cld(i,j,k) = g * ( btm_cld(i,j,k-1)*dtl_dsc(i,j)             &
                          + bqm_cld(i,j,k-1)*dqw_dsc(i,j) )

        chi_s(i,j,k) = -qcl_ic_top(i,j,k) /                             &
                 ( a_qsm(i,j,k-1)*dqw_dsc(i,j)                          &
                         - a_dqsdtm(i,j,k-1)*dtl_dsc(i,j) )
      END IF

    END IF  ! disc inversion diagnosed

  END DO
END DO
!-----------------------------------------------------------------------
!  3.0a Calculate surface buoyancy flux
!-----------------------------------------------------------------------

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    ! use mixed-layer average of buoyancy parameters
    bflux_surf(i,j) = 0.5 * g * (                                       &
         (btm(i,j,1)+btm(i,j,ntml(i,j)))*ftl(i,j,1) +                   &
         (bqm(i,j,1)+bqm(i,j,ntml(i,j)))*fqw(i,j,1) )

    IF ( bflux_surf(i,j)   >   0.0 ) THEN
      bflux_surf_sat(i,j) = 0.5 * g * (                                 &
         (btm_cld(i,j,1)+btm_cld(i,j,ntml(i,j)))*ftl(i,j,1) +           &
         (bqm_cld(i,j,1)+bqm_cld(i,j,ntml(i,j)))*fqw(i,j,1) )
      IF ( coupled(i,j) ) bflux_surf_sat(i,j) = 0.5 * g * (             &
         (btm_cld(i,j,1)+btm_cld(i,j,ntdsc(i,j)))*ftl(i,j,1) +          &
         (bqm_cld(i,j,1)+bqm_cld(i,j,ntdsc(i,j)))*fqw(i,j,1) )
    ELSE
      bflux_surf_sat(i,j) = 0.0
    END IF

  END DO
END DO
!-----------------------------------------------------------------------
! 3.0aa Calculate uniform mixed-layer cloud fraction and thence
!       estimate Sc layer cloud depth (not cloud fraction weighted).
!       (If DSC=.FALSE. then NTDSC=0 and ZC_DSC remains equal to 0.)
!-----------------------------------------------------------------------
! First the SML
!--------------
! First find cloud-base as seen by the cloud scheme, at grid-level
! K_LEVEL and height Z_CF_BASE, to use as first guess or lower limit
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    k_level(i,j) = ntml(i,j)
    DO WHILE ( cf(i,j,k_level(i,j))  >   sc_cftol                       &
                .AND. k_level(i,j)  >=  2 )
      k_level(i,j) = k_level(i,j) - 1
    END DO
  END DO
END DO
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    cloud_base(i,j)= .FALSE.
    zc(i,j)        = 0.0
    k_cbase(i,j)   = 0
    z_cf_base(i,j) = zh(i,j)
    z_ctop(i,j)    = zh(i,j)
      ! Use a single CF for whole mixed-layer (more realistic).
      ! Include NTML+1 if a subgrid inversion has been diagnosed
    IF ( coupled(i,j) .OR. cumulus(i,j) ) THEN
      cf_sml(i,j)=0.0
    ELSE
      k = ntml(i,j)
      cf_sml(i,j) = MAX( cf(i,j,k), cf(i,j,k+sml_disc_inv(i,j)) )
    END IF
  END DO
END DO

! initialisation of z_cf_base and first guess for cloud depth, ZC

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( cf_sml(i,j)  >   sc_cftol ) THEN
      k = ntml(i,j)
      IF ( sml_disc_inv(i,j)  ==  0 .AND.                               &
           cf(i,j,k+1)  >   sc_cftol) THEN
          ! if no subgrid inversion and level NTML+1 is cloudy
        z_ctop(i,j) = z_top(i,j,k+1)
      END IF
      IF ( k_level(i,j)  ==  1 .AND.                                    &
           cf(i,j,k_level(i,j))  >   sc_cftol) THEN
        z_cf_base(i,j) = 0.0
      ELSE
        z_cf_base(i,j) = z_uv(i,j,k_level(i,j)+1)
      END IF
      zc(i,j) = z_ctop(i,j) - z_cf_base(i,j)
    END IF
  END DO
END DO
!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------
DO k = bl_levels, 1, -1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( .NOT. cloud_base(i,j) .AND. k  <=  ntml(i,j)+1 .AND.         &
           cf_sml(i,j)  >   sc_cftol ) THEN
             ! within cloudy boundary layer
        IF ( k  ==  1) THEN
          cloud_base(i,j) = .TRUE.
        ELSE
          IF ( cf(i,j,k-1)  <   cf(i,j,k) ) cloud_base(i,j) = .TRUE.
        END IF
        k_cbase(i,j) = k
      END IF
    END DO
  END DO
END DO

!Initialise K_CFF = lowest level with ice cloud
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    k_cff(i,j) = k_cbase(i,j)
    IF (k_cff(i,j) > 1) THEN
      DO WHILE ( cff(i,j,k_cff(i,j))  >   sc_cftol                      &
               .AND. k_cff(i,j)  >   1)
        k_cff(i,j) = k_cff(i,j) - 1
      END DO
    END IF
  END DO
END DO

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

        !--------------------------------------------------
        ! Use adiabatic qcl gradient to estimate cloud-base
        ! from in-cloud qcl in level K_CBASE
        ! If k_cbase = 0 then it hasn't been initialised
        !--------------------------------------------------

    IF ( cloud_base(i,j) .AND. k_cbase(i,j) /= 0 ) THEN
      z_cbase = z_tq(i,j,k_cbase(i,j)) -                                &
                qcl(i,j,k_cbase(i,j)) /                                 &
                ( cf(i,j,k_cbase(i,j))*dqcldz(i,j,k_cbase(i,j)) )
      IF ( dqcfdz(i,j,k_cbase(i,j))  >   0.0 ) THEN
        z_cbase = MIN( z_cbase, z_tq(i,j,k_cbase(i,j)) -                &
               qcf(i,j,k_cbase(i,j)) /                                  &
               ( cf(i,j,k_cbase(i,j))*dqcfdz(i,j,k_cbase(i,j)) )        &
                     )
      ELSE
            !---------------------------------------------------------
            ! No adiabatic QCF gradient so find lowest level, K_CFF,
            ! with CFF>SC_CFTOL and assume cloud-base within that leve
            !---------------------------------------------------------
        IF ( cff(i,j,k_cff(i,j))  <=  sc_cftol .AND.                    &
                      k_cff(i,j)  <   k_cbase(i,j) )                    &
             k_cff(i,j) = k_cff(i,j) + 1
                 ! will want to raise K_CFF back up one level unless
                 ! level 1 is cloudy or no sig frozen cloud at all
        z_cbase = MIN( z_cbase, z_top(i,j,k_cff(i,j)) -                 &
                  dzl(i,j,k_cff(i,j))                                   &
                * cff(i,j,k_cff(i,j))/cf(i,j,k_cff(i,j)) )
      END IF
          !------------------------------------------------------
          ! use cloud-base as seen by cloud scheme as lower limit
          ! and base of level NTML+1 as upper limit
          !------------------------------------------------------
      z_cbase = MIN( z_uv(i,j,ntml(i,j)+1),                             &
                     MAX( z_cf_base(i,j), z_cbase) )

      zc(i,j) = z_ctop(i,j) - z_cbase
    END IF

  END DO
END DO
!-----------------------------------------------------------------------
! Second DSC layer
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    k_level(i,j) = ntdsc(i,j)
    IF ( k_level(i,j) >= 2) THEN
      DO WHILE ( cf(i,j,k_level(i,j))  >   sc_cftol                     &
               .AND. k_level(i,j)  >=  2 )
        k_level(i,j) = k_level(i,j) - 1
      END DO
    END IF
  END DO
END DO

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    cloud_base(i,j) = .FALSE.
    zc_dsc(i,j) = 0.0
    cf_dsc(i,j) = 0.0
    k_cbase(i,j) = 0
    z_cf_base(i,j) = zhsc(i,j)
    z_ctop(i,j)    = zhsc(i,j)

    IF ( dsc(i,j) ) THEN
      k = ntdsc(i,j)
      cf_dsc(i,j) = MAX( cf(i,j,k), cf(i,j,k+1) )
        !-------------------------------------------------------------
        ! Find cloud-base as seen by cloud scheme, Z_CF_BASE,
        ! to use as first guess or lower limit and find cloud top.
        !-------------------------------------------------------------
      IF ( cf_dsc(i,j)  >   sc_cftol ) THEN
        k = ntdsc(i,j)
        IF ( dsc_disc_inv(i,j)  ==  0 .AND.                             &
             cf(i,j,k+1)  >   sc_cftol) THEN
            ! if no subgrid inversion and level NTML+1 is cloudy
            ! then include this layer in cloud-depth
          z_ctop(i,j) = z_top(i,j,k+1)
        END IF
        IF ( k_level(i,j)  ==  1 .AND.                                  &
              cf(i,j,k_level(i,j))  >   sc_cftol) THEN
          z_cf_base(i,j) = 0.0
        ELSE
          z_cf_base(i,j) = z_uv(i,j,k_level(i,j)+1)
        END IF
        zc_dsc(i,j) = z_ctop(i,j) - z_cf_base(i,j)   ! first guess
      END IF
    END IF
  END DO
END DO
!--------------------------------------------------
! Find lowest level within ML with max CF
!--------------------------------------------------
DO k = bl_levels, 1, -1
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( .NOT. cloud_base(i,j) .AND. k  <=  ntdsc(i,j)+1 .AND.        &
           cf_dsc(i,j)  >   sc_cftol ) THEN
             ! within cloudy boundary layer
        IF ( k  ==  1) THEN
          cloud_base(i,j) = .TRUE.
        ELSE
          IF ( cf(i,j,k-1)  <   cf(i,j,k) ) cloud_base(i,j) = .TRUE.
        END IF
        k_cbase(i,j) = k
      END IF
    END DO ! I
  END DO ! J
END DO ! K

! Initialise K_CFF
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    k_cff(i,j) = k_cbase(i,j)
    IF (k_cff(i,j) > 1) THEN
      DO WHILE ( cff(i,j,k_cff(i,j))  >   sc_cftol                      &
                    .AND. k_cff(i,j)  >   1)
        k_cff(i,j) = k_cff(i,j) - 1
      END DO
    END IF
  END DO
END DO

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

      !--------------------------------------------------
      ! use adiabatic qcl gradient to estimate cloud-base
      ! from in-cloud qcl in level K_CBASE
      !--------------------------------------------------
    IF ( cloud_base(i,j) .AND. k_cbase(i,j)  /= 0 ) THEN
      z_cbase = z_tq(i,j,k_cbase(i,j)) -                                &
                qcl(i,j,k_cbase(i,j)) /                                 &
                ( cf(i,j,k_cbase(i,j))*dqcldz(i,j,k_cbase(i,j)) )
      IF ( dqcfdz(i,j,k_cbase(i,j))  >   0.0 ) THEN
        z_cbase = MIN( z_cbase, z_tq(i,j,k_cbase(i,j)) -                &
              qcf(i,j,k_cbase(i,j)) /                                   &
              ( cf(i,j,k_cbase(i,j))*dqcfdz(i,j,k_cbase(i,j)) )         &
                     )
      ELSE
          !----------------------------------------------------------
          ! No adiabatic QCF gradient so find lowest level, K_CFF,
          ! with CFF>SC_CFTOL and assume cloud-base within that level
          !----------------------------------------------------------
        IF ( cff(i,j,k_cff(i,j))  <=  sc_cftol .AND.                    &
                      k_cff(i,j)  <   k_cbase(i,j) )                    &
             k_cff(i,j) = k_cff(i,j) + 1
               ! will want to raise K_CFF back up one level unless
               ! level 1 is cloudy or no sig frozen cloud at all
        z_cbase = MIN( z_cbase, z_top(i,j,k_cff(i,j)) -                 &
                  dzl(i,j,k_cff(i,j))                                   &
                 * cff(i,j,k_cff(i,j))/cf(i,j,k_cff(i,j)) )
      END IF
        !------------------------------------------------------
        ! use cloud-base as seen by cloud scheme as lower limit
        ! and base of level NTDSC+1 as upper limit
        !------------------------------------------------------
      z_cbase = MIN( z_uv(i,j,ntdsc(i,j)+1),                            &
                     MAX( z_cf_base(i,j) , z_cbase) )

      zc_dsc(i,j) = z_ctop(i,j) - z_cbase
    END IF

  END DO !I
END DO !J
    !-----------------------------------------------------------------
    !  Layer cloud depth cannot be > the layer depth itself.
    !-----------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    zc_dsc(i,j) = MIN( zc_dsc(i,j), dscdepth(i,j) )
  END DO
END DO
!-----------------------------------------------------------------------
!- Calculate buoyancy flux factor used in the diagnosis of decoupling:
!-----------------------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE)                       &
!$OMP& SHARED(bl_levels,pdims, z_tq, cumulus, coupled,                  &
!$OMP& ntml, zh, ntdsc, zhsc, qw, tl, rdz, ntml_prev,                   &
!$OMP& grad_t_adj, grad_q_adj, cf_sml, zc, zc_dsc,                      &
!$OMP& cf_dsc, btm, bqm, bqm_cld, btm_cld,                              &
!$OMP& db_ksurf_dry, db_ksurf_cld, db_ktop_dry, db_ktop_cld, g, grcp)   &
!$OMP& PRIVATE(k,j,i,kl,                                                &
!$OMP& z_int_top,km1,z_int,dqw,dtl,dqw_ga,dtl_ga,cf_for_wb,             &
!$OMP& z_cbase,zdsc_cbase)
DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      kl = k

      z_int_top = z_tq(i,j,k)
      IF (.NOT. cumulus(i,j) .AND. .NOT. coupled(i,j) .AND.             &
          k  ==  ntml(i,j)+1 .AND. ntml(i,j)  >=  2) THEN
        kl = ntml(i,j)
        z_int_top = zh(i,j)
      ELSE IF (k  ==  ntdsc(i,j)+1) THEN
        kl = ntdsc(i,j)
        z_int_top = zhsc(i,j)
      END IF
      km1 = kl-1
      z_int = z_int_top - z_tq(i,j,k-1)  ! integration depth

      dqw = qw(i,j,kl) - qw(i,j,km1)
      dtl = tl(i,j,kl) - tl(i,j,km1) + grcp/rdz(i,j,kl)
      IF (kl  <=  ntml_prev(i,j)) THEN
        dtl_ga = dtl - grad_t_adj(i,j)/rdz(i,j,kl)
        dqw_ga = dqw - grad_q_adj(i,j)/rdz(i,j,kl)
      ELSE
        dtl_ga = dtl
        dqw_ga = dqw
      END IF
        !----------------------------------------------------------
        ! CF_FOR_WB is uniform `bl' CF for use within cloud layers
        !----------------------------------------------------------
      cf_for_wb = 0.0
      z_cbase = zh(i,j)-zc(i,j)
      zdsc_cbase = zhsc(i,j)-zc_dsc(i,j)
      IF ( k  <=  ntml(i,j)+1 .AND.                                     &
           z_tq(i,j,k)  >=  z_cbase) cf_for_wb = cf_sml(i,j)
      IF ( k  <=  ntdsc(i,j)+1 .AND.                                    &
           z_tq(i,j,k)  >=  zdsc_cbase) cf_for_wb = cf_dsc(i,j)
        !----------------------------------------------------------
        ! WB = -K_SURF*(DB/DZ - gamma_buoy) - K_TOP*DB/DZ
        ! This is integrated in EXCF_NL, iterating the K profiles.
        ! Here the relevant integrated DB/DZ factors are calculated
        !----------------------------------------------------------
      db_ksurf_dry(i,j,k) = - g*rdz(i,j,kl)*z_int*                      &
               ( btm(i,j,km1)*dtl_ga + bqm(i,j,km1)*dqw_ga )
      db_ktop_dry(i,j,k)  = - g*rdz(i,j,kl)*z_int*                      &
               ( btm(i,j,km1)*dtl + bqm(i,j,km1)*dqw )
      db_ksurf_cld(i,j,k) = - g*rdz(i,j,kl)*z_int*                      &
               ( btm_cld(i,j,km1)*dtl_ga + bqm_cld(i,j,km1)*dqw_ga )
      db_ktop_cld(i,j,k)  = - g*rdz(i,j,kl)*z_int*                      &
               ( btm_cld(i,j,km1)*dtl + bqm_cld(i,j,km1)*dqw )
        !-------------------------------------------------------
        ! Weight cloud layer factors with cloud fraction
        !-------------------------------------------------------
      db_ksurf_cld(i,j,k) = db_ksurf_dry(i,j,k)*(1.0-cf_for_wb) +       &
                            db_ksurf_cld(i,j,k)*cf_for_wb
      db_ktop_cld(i,j,k)  = db_ktop_dry(i,j,k)*(1.0-cf_for_wb) +        &
                            db_ktop_cld(i,j,k)*cf_for_wb

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

!-----------------------------------------------------------------------
!  3.1 Calculate inputs for the top of b.l. entrainment parametrization
!-----------------------------------------------------------------------

DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    zeta_r_dsc(i,j) = 0.0
    chi_s_dsct(i,j) = 0.0
    cld_factor_dsc(i,j) = 0.0
    bt_dsct(i,j) = 0.0
    btt_dsct(i,j) = 0.0
    btc_dsct(i,j) = 0.0
    db_dsct(i,j) = 0.0
    db_dsct_cld(i,j) = 0.0
    chi_s_top(i,j) = 0.0
    cld_factor(i,j) = 0.0
    bt_top(i,j) = 0.0
    btt_top(i,j) = 0.0
    btc_top(i,j) = 0.0
    db_top(i,j) = 0.0
    db_top_cld(i,j) = 0.0    ! default required if COUPLED
    z_cld(i,j) = 0.0
    z_cld_dsc(i,j) = 0.0
  END DO
END DO

DO k = 1, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( k  <=  ntml(i,j)+1 ) THEN
          !---------------------------------------------------
          ! Calculation of cloud fraction weighted
          ! thickness of cloud in the surface mixed layer
          !---------------------------------------------------
        z_cld(i,j) = z_cld(i,j) +                                       &
                   cf(i,j,k) * 0.5 * dzl(i,j,k) +                       &
                    MIN( cfl(i,j,k) * 0.5 * dzl(i,j,k) ,                &
                            qcl(i,j,k) / dqcldz(i,j,k) )
        IF ( dqcfdz(i,j,k)  >   0.0) THEN
          z_cld(i,j) = z_cld(i,j) +                                     &
                    MIN( cff(i,j,k) * 0.5 * dzl(i,j,k) ,                &
                            qcf(i,j,k) / dqcfdz(i,j,k) )
        ELSE
          z_cld(i,j) = z_cld(i,j) + cff(i,j,k) * 0.5 * dzl(i,j,k)
        END IF
      END IF

      IF ( dsc(i,j) .AND. k <= ntdsc(i,j)+1 .AND.                       &
           ( coupled(i,j) .OR.                                          &
                 z_top(i,j,k) >= zhsc(i,j)-zc_dsc(i,j) ) ) THEN
          !-----------------------------------------------------------
          ! Calculation of cloud fraction weighted thickness of
          ! cloud in the DSC layer (or to the surface if COUPLED)
          !-----------------------------------------------------------
        z_cld_dsc(i,j) = z_cld_dsc(i,j) +                               &
                   cf(i,j,k) * 0.5 * dzl(i,j,k) +                       &
                    MIN( cfl(i,j,k) * 0.5 * dzl(i,j,k) ,                &
                            qcl(i,j,k) / dqcldz(i,j,k) )
        IF ( dqcfdz(i,j,k)  >   0.0) THEN
          z_cld_dsc(i,j) = z_cld_dsc(i,j) +                             &
                    MIN( cff(i,j,k) * 0.5 * dzl(i,j,k) ,                &
                            qcf(i,j,k) / dqcfdz(i,j,k) )
        ELSE
          z_cld_dsc(i,j) = z_cld_dsc(i,j) +                             &
                                  cff(i,j,k) * 0.5 * dzl(i,j,k)
        END IF
      END IF

      IF (k  ==  ntml(i,j)  .AND. .NOT. coupled(i,j) ) THEN
          !------------------------------------------------------
          ! Calculation of SML inputs.  If COUPLED then these are
          ! not used (as no entrainment is then applied at ZH)
          !------------------------------------------------------
        chi_s_top(i,j) = MAX( 0.0, MIN( chi_s(i,j,k+1), 1.0) )
        kp2=MIN(k+1+sml_disc_inv(i,j),bl_levels)
        cld_factor(i,j) = MAX( 0.0 , cf(i,j,k)-cf(i,j,kp2) )
        bt_top(i,j) = g * btm(i,j,k)
        btt_top(i,j) = g * btm_cld(i,j,k)
        btc_top(i,j) = btt_top(i,j)
        db_top(i,j) = db(i,j,k+1)
        db_top_cld(i,j) = db_cld(i,j,k+1)
        IF ( db_top(i,j)  <   0.001 ) THEN
            ! Diagnosed inversion statically unstable:
            ! ensure DB>0 so that entrainment is non-zero and
            ! instability can be removed.
          db_top(i,j) = 0.001
          db_top_cld(i,j) = 0.0  ! set buoyancy reversal
          chi_s_top(i,j) = 0.0   ! term to zero
        END IF
      END IF
      IF (k  ==  ntdsc(i,j)) THEN
          !---------------------------------------------------
          ! Calculation of DSC inputs
          ! (if DSC=.FALSE. then K never equals NTDSC(=0))
          !---------------------------------------------------
        chi_s_dsct(i,j) = MAX( 0.0, MIN( chi_s(i,j,k+1), 1.0) )
        kp2=MIN(k+1+dsc_disc_inv(i,j),bl_levels)
        cld_factor_dsc(i,j) = MAX( 0.0 , cf(i,j,k)-cf(i,j,kp2) )
        bt_dsct(i,j) = g * btm(i,j,k)
        btt_dsct(i,j) = g * btm_cld(i,j,k)
        btc_dsct(i,j) = btt_dsct(i,j)
        db_dsct(i,j) = db(i,j,k+1)
        db_dsct_cld(i,j) = db_cld(i,j,k+1)
        IF ( db_dsct(i,j)  <   0.001 ) THEN
            ! Diagnosed inversion statically unstable:
            ! ensure DB>0 so that entrainment is non-zero and
            ! instability can be removed.
          db_dsct(i,j) = 0.001
          db_dsct_cld(i,j) = 0.0  ! set buoyancy reversal
          chi_s_dsct(i,j) = 0.0   ! term to zero
        END IF
      END IF
    END DO
  END DO
END DO
!-----------------------------------------------------------------------
! Next those terms which depend on the presence of buoyancy reversal
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    z_cld(i,j) = MIN( z_cld(i,j), zh(i,j) )
    z_cld_dsc(i,j) = MIN( z_cld_dsc(i,j), zhsc(i,j) )
      !---------------------------------------------------------------
      ! First the surface mixed layer.
      !---------------------------------------------------------------
    IF ( coupled(i,j) ) THEN
      zeta_s(i,j) = 1.0 - z_cld_dsc(i,j) / zhsc(i,j)
      zeta_r(i,j) = 1.0 - zc_dsc(i,j) / zhsc(i,j)
    ELSE
      zeta_s(i,j) = 1.0 - z_cld(i,j) / zh(i,j)
      zeta_r(i,j) = 1.0 - zc(i,j) / zh(i,j)
    END IF

    IF (db_top_cld(i,j)  >=  0.0) THEN
        !--------------------------------------------------
        ! i.e. no buoyancy reversal (or default if COUPLED)
        !--------------------------------------------------
      db_top_cld(i,j) = 0.0
      d_siems(i,j) = 0.0
    ELSE
        !----------------------------
        ! IF (DB_TOP_CLD(I,j)  <   0.0)
        ! i.e. buoyancy reversal
        !----------------------------
      db_top_cld(i,j) = -db_top_cld(i,j) * cld_factor(i,j)
      d_siems(i,j) = MAX( 0.0,                                          &
           chi_s_top(i,j) * db_top_cld(i,j) / (db_top(i,j)+1.0e-14) )
      zeta_r(i,j) = MIN( zeta_r(i,j)+10.0*(1.0-zeta_r(i,j))             &
                                      *d_siems(i,j), 1.0 )
    END IF
      !---------------------------------------------------------------
      ! Now the decoupled Sc layer (DSC).
      !---------------------------------------------------------------
    IF (dsc(i,j)) THEN
      IF ( coupled(i,j) ) THEN
        zeta_r_dsc(i,j) = 1.0 - zc_dsc(i,j) / zhsc(i,j)
      ELSE
        zeta_r_dsc(i,j) = 1.0 - zc_dsc(i,j) / dscdepth(i,j)
      END IF

      IF (db_dsct_cld(i,j)  >=  0.0) THEN
          !----------------------------
          ! i.e. no buoyancy reversal
          !----------------------------
        db_dsct_cld(i,j) = 0.0
        d_siems_dsc(i,j) = 0.0
      ELSE
          !----------------------------
          ! IF (DB_DSCT_CLD(I,j)  <   0.0)
          ! i.e. buoyancy reversal
          !----------------------------
        db_dsct_cld(i,j) = -db_dsct_cld(i,j) * cld_factor_dsc(i,j)
        d_siems_dsc(i,j) = MAX( 0.0, chi_s_dsct(i,j)                    &
                        * db_dsct_cld(i,j) / (db_dsct(i,j)+1.0e-14) )
        zeta_r_dsc(i,j) = MIN( zeta_r_dsc(i,j) +                        &
                10.0*(1.0-zeta_r_dsc(i,j)) * d_siems_dsc(i,j), 1.0 )
      END IF
    END IF
  END DO
END DO

!-----------------------------------------------------------------------
! 4. Calculate the radiative flux change across cloud top for mixed-
!    layer to ZH.  Restrict search for maximum divergence to below
!    NTML+2.  This may introduce errors if NTML changes a lot during
!    the radiative timestep but can't be helped.
!-----------------------------------------------------------------------
!     Initialise variables
!------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    k_cloud_top(i,j) = 0
    df_top_over_cp(i,j) = 0.0
  END DO
END DO

DO k = 1, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      ! *APL*restrict search to `close' to ZH
      k_rad_smlt = ntml(i,j)+1
        !-------------------------------------------------------------
        ! Find the layer with the LW greatest rad flux jump below
        ! K_RAD_SMLT and assume that this is the top of the SML.
        !-------------------------------------------------------------
      IF (dflw_over_cp(i,j,k)  >   df_top_over_cp(i,j)                  &
                  .AND. k  <=  k_rad_smlt ) THEN
        k_cloud_top(i,j) = k
        df_top_over_cp(i,j) = dflw_over_cp(i,j,k)
      END IF

    END DO
  END DO
END DO
    !-----------------------------------------------------------------
    !  In case the radiative divergence is spread over 2 or 3 levels,
    !  add on any cooling either side of the level of maximum cooling.
    !-----------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( k_cloud_top(i,j)  >   0 ) THEN
      k=k_cloud_top(i,j)
      dfsw_top = dfsw_over_cp(i,j,k)
      IF ( k >  1 ) THEN
        df_top_over_cp(i,j) = df_top_over_cp(i,j)                       &
                              + MAX(0.0, dflw_over_cp(i,j,k-1) )
        dfsw_top = dfsw_top + MIN(0.0, dfsw_over_cp(i,j,k-1) )
      END IF
      IF ( k <  bl_levels ) THEN
        df_top_over_cp(i,j) = df_top_over_cp(i,j)                       &
                              + MAX(0.0, dflw_over_cp(i,j,k+1) )
        dfsw_top = dfsw_top + MIN(0.0, dfsw_over_cp(i,j,k+1) )
      END IF
          !-----------------------------------------------------------
          ! Combine SW and LW cloud-top divergences into a net
          ! divergence by estimating SW flux divergence at LW
          ! extinction depth = DF_SW * (1-exp{-A*kappa_sw/kappa_lw})
          ! Choose A=3 (to achieve 95% extinction of LW) and
          ! kappa_lw = 3*kappa_sw
          !-----------------------------------------------------------
      df_top_over_cp(i,j) = MAX( 0.0,                                   &
                  df_top_over_cp(i,j) + 0.35 * dfsw_top )
    END IF
  END DO
END DO
!-----------------------------------------------------------------------
! 5.1 Subroutine EXCF_NL.
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ntml_save(i,j) = ntml(i,j)  ! needed to identify changes
    dsc_save(i,j) = dsc(i,j)    !      in excf_nl
  END DO
END DO

! DEPENDS ON: excf_nl_9b
CALL excf_nl_9b (                                                       &
! IN levels
   bl_levels,                                                           &
! IN fields
   rdz,z_uv,z_tq,rho_mix,rho_wet_tq,rhostar_gb,v_s,fb_surf,db_top,      &
   bflux_surf,bflux_surf_sat,zeta_s,bt_top,btt_top,                     &
   df_top_over_cp,zeta_r,btc_top, db_top_cld,chi_s_top,                 &
   bt_dsct,btt_dsct,df_dsct_over_cp,zeta_r_dsc,                         &
   db_dsct_cld,chi_s_dsct,d_siems,d_siems_dsc,                          &
   db_ksurf_dry,db_ktop_dry,db_ksurf_cld,db_ktop_cld,db_dsct,           &
! INOUT fields
   coupled,cumulus,dsc,ntml,ntdsc,zh,zc,zhsc,dscdepth,zc_dsc,BL_diag,   &
! OUT fields
   rhokm,rhokh,rhokm_top,rhokh_top, f_ngstress,zdsc_base,               &
   rhokh_top_ent,rhokh_dsct_ent,rhokh_surf_ent,nbdsc                    &
  )

!-----------------------------------------------------------------------
!-adjust SML/DSC properties depending on diagnoses in EXCF_NL
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    IF ( dsc(i,j) .AND. .NOT. dsc_save(i,j) ) THEN
      !..decoupling diagnosed in EXCF_NL - change parameters around
      dtl_dsc(i,j) = dtl_sml(i,j)
      dqw_dsc(i,j) = dqw_sml(i,j)
      dsc_disc_inv(i,j) = sml_disc_inv(i,j)
      sml_disc_inv(i,j) = 0
      dtl_sml(i,j) = 0.0
      dqw_sml(i,j) = 0.0
      k_cloud_dsct(i,j) = k_cloud_top(i,j)
      k_cloud_top(i,j) = 0
    END IF
    IF ( .NOT. dsc(i,j) .AND. dsc_save(i,j) ) THEN
      !..decoupled layer removed in EXCF_NL; either...
      IF ( ntml_save(i,j)  ==  ntml(i,j) ) THEN
        !...had no turbulence forcing
        dtl_dsc(i,j) = 0.0
        dqw_dsc(i,j) = 0.0
        dsc_disc_inv(i,j) = 0
        k_cloud_dsct(i,j) = 0
      ELSE
        !...recoupled with surface layer
        dtl_sml(i,j) = dtl_dsc(i,j)
        dqw_sml(i,j) = dqw_dsc(i,j)
        dtl_dsc(i,j) = 0.0
        dqw_dsc(i,j) = 0.0
        sml_disc_inv(i,j) = dsc_disc_inv(i,j)
        dsc_disc_inv(i,j) = 0
        k_cloud_top(i,j) = k_cloud_dsct(i,j)
        k_cloud_dsct(i,j) = 0
      END IF
    END IF
  END DO
END DO
!-----------------------------------------------------------------------
!  6.  Calculate "explicit" entrainment fluxes of TL and QW.
!-----------------------------------------------------------------------
!..Integrate radiative divergence over layers:
!..  DF_NET (SML,DSC) are the total net radiative forcings of the
!..                   layers, used to calculate the total heat fluxes
!..                   at the discontinuous inversion levels
!..  F_NET (SML,DSC) are the net radiative divergences from the base
!..                  of the turbulently mixed layers, used to calculate
!..                  the appropriate grid-level turbulent fluxes
!..  All are in units of rho * Km/s
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    df_net_sml(i,j) = 0.0
    df_net_dsc(i,j) = 0.0
    f_net_sml(i,j,1) = 0.0
    f_net_dsc(i,j,1) = 0.0
      ! If no rad cooling in layer, need to set K_CLOUD_TOP to NTML
      ! so that the correct entrainment heat flux is specified,
      ! consistent with a linear total heat flux profile.
    IF ( k_cloud_top(i,j)  ==  0 ) k_cloud_top(i,j) = ntml(i,j)
    IF ( k_cloud_dsct(i,j)  ==  0 ) k_cloud_dsct(i,j) = ntdsc(i,j)
  END DO
END DO
DO k = 1, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      !..First the surface-forced mixed layer
      IF ( sml_disc_inv(i,j)  ==  1 ) THEN
        IF ( k  <=  k_cloud_top(i,j) )                                  &
          df_net_sml(i,j) = df_net_sml(i,j) + df_over_cp(i,j,k)
      END IF
      !..Second the decoupled mixed layer
      IF ( dsc_disc_inv(i,j)  ==  1 ) THEN
        IF ( k  <=  k_cloud_dsct(i,j) .AND. k  >=  nbdsc(i,j) )         &
           df_net_dsc(i,j) = df_net_dsc(i,j) + df_over_cp(i,j,k)
      END IF
    END DO
  END DO
END DO
!-----------------------------------------------------------------------
!..Assume DF_OVER_CP(K_cloud_top+2) is representative of clear-air rad
!..divergence and so extrapolate this down to discontinuous inversion
!..height and subtract this `clear-air' part from the grid-level
!..divergence.
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    !..First the surface-forced mixed layer
    IF ( sml_disc_inv(i,j)  ==  1 ) THEN
      k = k_cloud_top(i,j)+1
      IF ( k  <   bl_levels ) THEN
        df_inv = df_over_cp(i,j,k) - df_over_cp(i,j,k+1) *              &
                       ( z_uv(i,j,ntml(i,j)+2) - zh(i,j) )              &
                       / dzl(i,j,k+1)
      ELSE IF ( k  ==  bl_levels ) THEN
        df_inv = df_over_cp(i,j,k)
      ELSE
        df_inv = 0.0
      END IF
      df_net_sml(i,j) = df_net_sml(i,j) + MAX( df_inv, 0.0 )
    END IF
    !..Second the decoupled mixed layer
    IF ( dsc_disc_inv(i,j)  ==  1 ) THEN
      k = k_cloud_dsct(i,j)+1
      IF ( k  <   bl_levels ) THEN
        df_inv = df_over_cp(i,j,k) - df_over_cp(i,j,k+1) *              &
                       ( z_uv(i,j,ntdsc(i,j)+2) - zhsc(i,j) )           &
                       / dzl(i,j,k+1)
      ELSE IF ( k  ==  bl_levels ) THEN
        df_inv = df_over_cp(i,j,k)
      ELSE
        df_inv = 0.0
      END IF
      df_net_dsc(i,j) = df_net_dsc(i,j) + MAX( df_inv, 0.0 )
    END IF
  END DO
END DO

!..Now the net radiative flux profiles from the layer base

DO k = 2, bl_levels
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      !..First the surface-forced mixed layer
      IF ( sml_disc_inv(i,j)  ==  1 ) THEN
        f_net_sml(i,j,k) = f_net_sml(i,j,k-1) + df_over_cp(i,j,k-1)
      END IF
      !..Second the decoupled mixed layer
      IF ( dsc_disc_inv(i,j)  ==  1 ) THEN
        IF ( k  <   nbdsc(i,j) ) THEN
          f_net_dsc(i,j,k) = 0.0
        ELSE
          f_net_dsc(i,j,k) = f_net_dsc(i,j,k-1) + df_over_cp(i,j,k-1)
        END IF
      END IF
    END DO
  END DO
END DO
!-----------------------------------------------------------------------
!..Specify entrainment fluxes at NTML+1 and NTDSC+1 directly through FTL
!..and FQW (and set the entrainment RHOKH to zero).
!..For QW, assume turbulent flux at subgrid ZH = - w_e * DQW_SML
!..The flux at the half-level below ZH can then be found by linear
!..interpolation between ZH and the surface (or base of the mixed-layer)
!..For TL, the procedure is basically the same, except that it is
!..the total (turb+rad) heat flux, TOTHF, that is linear.  Thus, TOTHF
!..at ZH is calculated, which is dependent on both the entrainment flux
!..and the net radiative divergence over the mixed layer, DF_NET_SML,
!..the linear profile of TOTHF is interpolated onto the half-level below
!..ZH and the radiative flux there, DF_NET_NTML, is sutracted off to
!..give the required turbulent component.  This has the happy
!..coincidence of ensuring that the `correct' entrainment flux is
!..specified for the given radiative flux (it even compensates for the
!..radiative flux having got a grid-level out, since radiation was last
!..called, for example).
!..For momentum, given the horizontal interpolation required, together
!..with the lack of accuracy in assuming a discontinuous inversion,
!..entrainment continues to be specified using the specified RHOKM.
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end

    !..First the surface-based mixed layer (if entraining)
    k=ntml(i,j)+1
    kent(i,j) = 2
    zh_np1(i,j) = 0.0
    DO ient = 1, 3
      t_frac(i,j,ient) = 0.0
      zrzi(i,j,ient)   = 0.0
      we_lim(i,j,ient) = 0.0
    END DO
    we_parm(i,j) = rdz(i,j,k)*                                          &
                       ( rhokh_top_ent(i,j)+rhokh_surf_ent(i,j) )       &
                                                  / rho_mix(i,j,k)

    IF ( sml_disc_inv(i,j)  ==  1 .AND. .NOT. coupled(i,j) .AND.        &
         (rhokh_top_ent(i,j)+rhokh_surf_ent(i,j))  >   0.0 ) THEN

      kent(i,j) = k
      !-----------------------------------------------------------------------
      !..Calculate ZH at end of timestep, ZH_NP1
      !-----------------------------------------------------------------------
      !..linearly interpolate vertical velocity to ZH
      IF ( zh(i,j)  >=  z_tq(i,j,k) ) THEN
        w_ls = w(i,j,k) + ( w(i,j,k+1) - w(i,j,k) )                     &
                    * (zh(i,j)-z_tq(i,j,k)) * rdz(i,j,k+1)
      ELSE
        w_ls = w(i,j,k) + ( w(i,j,k) - w(i,j,k-1) )                     &
                    * (zh(i,j)-z_tq(i,j,k)) * rdz(i,j,k)
      END IF
      w_ls = MIN ( w_ls, 0.0 )
        ! only interested in subsidence

      zh_np1(i,j) = zh(i,j) +                                           &
                      timestep * ( we_parm(i,j) + w_ls )
      zh_np1(i,j) = MAX( zh_np1(i,j), z_uv(i,j,k-1) )
      IF ( zh_np1(i,j)  >   z_top(i,j,k+1) ) THEN
          ! limit ZH and W_e (and therefore the entraiment fluxes)
          ! because the inversion cannot rise more than one level
          ! in a timestep.
        zh_np1(i,j) = z_top(i,j,k+1)
        we_parm(i,j) =                                                  &
                (z_top(i,j,k+1) - zh(i,j))/timestep - w_ls
      END IF

      !..Linearly interpolate between the known total (turb+rad) TL fluxes at
      !..the surface and the mean height of the discontinuous inversion,
      !..taking care over where the inversion is during the timestep.

      IF ( zh_np1(i,j)  >   z_uv(i,j,k+1) ) THEN
        ! ZH risen above level K+1 so specify appropriate flux at this leve

               ! Adjust the entraiment rate to allow for that implied by
               ! the subsidence increment applied below (combine increments
               ! from levels NTML and NTML+1 to give the subsidence
               ! associated with the inversion.
        tls_inv(i,j) = tls_inc(i,j,k) + tls_inc(i,j,k-1)
        w_s_ent = 0.0
        IF ( dtl_sml(i,j)  /=  0.0 ) w_s_ent =                          &
               MIN( 0.0, -tls_inv(i,j) * dzl(i,j,k) /dtl_sml(i,j) )
          ! Only allow w_e to be reduced to zero!
        we_lim(i,j,3) = rho_mix(i,j,k+1) *                              &
                        MAX( 0.0, we_parm(i,j) + w_s_ent )

        tothf_inv = - we_lim(i,j,3)*dtl_sml(i,j) + df_net_sml(i,j)
        ! T_FRAC is fraction of timestep inversion spends above Z_UV(K+1)
        t_frac(i,j,3) = (zh_np1(i,j)-z_uv(i,j,k+1)) /                   &
                        (zh_np1(i,j)-zh(i,j))
        ! Estimate average flux over timestep as flux for average ZH
        zrzi(i,j,3)  = z_uv(i,j,k+1)*2.0/                               &
                       (z_uv(i,j,k+1)+zh_np1(i,j))

        ftl(i,j,k+1) = t_frac(i,j,3) * (                                &
           ftl(i,j,1) + ( tothf_inv - ftl(i,j,1) )*zrzi(i,j,3)          &
                                            - f_net_sml(i,j,k+1) )
        rhokh_top(i,j,k+1) = 0.0   ! apply entrainment explicitly
        rhokh(i,j,k+1) = 0.0       !      "

        ! Layer should now be well-mixed to NTML+1 (by end of timestep)
        ! so raise NTML by one (this means gradient-adjustment is also
        ! applied at half-level old_NTML+1).
        ! Note KH profiles should already be calculated at level NTML+1
        ! because ZH is above this level.
        ntml(i,j) = ntml(i,j) + 1

      ELSE  ! ZH always below half-level K+1

        rhokh_top(i,j,k) = 0.0   ! apply entrainment explicitly
        rhokh(i,j,k) = 0.0       !      "

        IF ( zh_np1(i,j)  >=  z_uv(i,j,k) ) THEN
          !..ZH always above level K so specify full entrainment flux here

               ! Adjust the entraiment rate to allow for that impled by the
               ! subsidence increment applied below.
          w_s_ent = 0.0
          IF ( dtl_sml(i,j)  /=  0.0 ) w_s_ent =                        &
            MIN( 0.0, -tls_inc(i,j,k-1) * dzl(i,j,k-1) /dtl_sml(i,j) )
          we_lim(i,j,2) = rho_mix(i,j,k) *                              &
                          MAX( 0.0, we_parm(i,j) + w_s_ent )

          tothf_inv = - we_lim(i,j,2)*dtl_sml(i,j) + df_net_sml(i,j)
          t_frac(i,j,2) = 1.0
          zrzi(i,j,2) = z_uv(i,j,k)*2.0/(zh(i,j)+zh_np1(i,j))

          ftl(i,j,k) = ftl(i,j,1) + ( tothf_inv - ftl(i,j,1) )          &
                                   *zrzi(i,j,2) - f_net_sml(i,j,k)
        ELSE
          !ZH dropped below level K so specify full entrainment flux at K-1
          IF (k-1  >=  2) THEN     ! ftl(k=1) is surface flux
            ntml(i,j) = ntml(i,j) - 1
            w_s_ent = 0.0
            IF ( dtl_sml(i,j)  /=  0.0 ) w_s_ent =                      &
             MIN( 0.0,                                                  &
                  -tls_inc(i,j,k-2) * dzl(i,j,k-2) /dtl_sml(i,j) )
            we_lim(i,j,1) = rho_mix(i,j,k-1) *                          &
                                   MAX( 0.0, we_parm(i,j)+w_s_ent )
            tothf_inv = - we_lim(i,j,1)*dtl_sml(i,j) + df_net_sml(i,j)
            t_frac(i,j,1) = 1.0
            zrzi(i,j,1) = z_uv(i,j,k-1)*2.0/(zh(i,j)+zh_np1(i,j))
            ftl(i,j,k-1) = ftl(i,j,1)                                   &
                          + (tothf_inv - ftl(i,j,1)) * zrzi(i,j,1)      &
                                            - f_net_sml(i,j,k-1)
            rhokh_top(i,j,k-1) = 0.0   ! apply entrainment explicitly
            rhokh(i,j,k-1) = 0.0      !      "
          END IF
          ! ...and set specified flux at level K allowing for time spent
          ! with ZH above this level

          w_s_ent = 0.0
          IF ( dtl_sml(i,j)  /=  0.0 ) w_s_ent =                        &
            MIN( 0.0, -tls_inc(i,j,k-1) * dzl(i,j,k-1) /dtl_sml(i,j) )
          we_lim(i,j,2) = rho_mix(i,j,k) *                              &
                           MAX( 0.0, we_parm(i,j) + w_s_ent )
          tothf_inv = - we_lim(i,j,2)*dtl_sml(i,j) + df_net_sml(i,j)
          ! T_FRAC is fraction of timestep inversion spends above Z_UV(K)
          t_frac(i,j,2) = (zh(i,j)-z_uv(i,j,k))                         &
                                     / (zh(i,j)-zh_np1(i,j))
          zrzi(i,j,2) = z_uv(i,j,k)*2.0/(z_uv(i,j,k)+zh(i,j))
          ftl(i,j,k) = t_frac(i,j,2) * (                                &
             ftl(i,j,1) + (tothf_inv - ftl(i,j,1)) * zrzi(i,j,2)        &
                                               - f_net_sml(i,j,k) )
        END IF

      END IF  ! ZH gone above half-level K+1.

    ELSE   ! NOT specifying entrainment flux but KH

        ! Add entrainment KH to K-profiles
        ! (for COUPLED layers these will be zero)
      rhokh_top(i,j,k) = rhokh_top(i,j,k) + rhokh_top_ent(i,j)
      rhokh(i,j,k) = rhokh(i,j,k) + rhokh_surf_ent(i,j)

    END IF  ! IF NOT CUMULUS OR COUPLED
    !-------------------------------------------------
    !..Second the decoupled mixed layer, if entraining
    !-------------------------------------------------
    k=ntdsc(i,j)+1
    kent_dsc(i,j) = 2
    zhsc_np1(i,j) = 0.0
    DO ient = 1, 3
      t_frac_dsc(i,j,ient) = 0.0
      zrzi_dsc(i,j,ient)   = 0.0
      we_lim_dsc(i,j,ient) = 0.0
    END DO
    we_dsc_parm(i,j) = rdz(i,j,k)*rhokh_dsct_ent(i,j)                   &
                                               / rho_mix(i,j,k)

    IF ( dsc_disc_inv(i,j)  ==  1                                       &
                  .AND. rhokh_dsct_ent(i,j)  >   0.0 ) THEN

      kent_dsc(i,j) = k
      !-----------------------------------------------------------------------
      !..Calculate ZHSC at end of timestep, ZHSC_NP1
      !-----------------------------------------------------------------------
      !..interpolate vertical velocity to ZH
      IF ( zhsc(i,j)  >=  z_tq(i,j,k) ) THEN
        w_ls_dsc = w(i,j,k) + ( w(i,j,k+1) - w(i,j,k) ) *               &
                         (zhsc(i,j)-z_tq(i,j,k)) * rdz(i,j,k+1)
      ELSE
        w_ls_dsc = w(i,j,k) + ( w(i,j,k) - w(i,j,k-1) ) *               &
                         (zhsc(i,j)-z_tq(i,j,k)) * rdz(i,j,k)
      END IF
      w_ls_dsc = MIN ( w_ls_dsc, 0.0 )
        ! only interested in subsidence

      zhsc_np1(i,j) = zhsc(i,j) +                                       &
            timestep * ( we_dsc_parm(i,j) + w_ls_dsc )
      zhsc_np1(i,j) = MAX( zhsc_np1(i,j), z_uv(i,j,k-1) )
      IF ( zhsc_np1(i,j)  >   z_top(i,j,k+1) ) THEN
          ! limit ZHSC and W_e (and therefore the entrainment fluxes)
          ! because the inversion cannot rise more than one level
          ! in a timestep.
        zhsc_np1(i,j) = z_top(i,j,k+1)
        we_dsc_parm(i,j) =                                              &
           (z_top(i,j,k+1) - zhsc(i,j))/timestep - w_ls_dsc
      END IF

      !..linearly interpolate between the known total (turb+rad) TL fluxes at
      !..the base (assumed zero) and the mean height of the discontinuous
      !..inversion, taking care over where the inversion is during the
      !..timestep.  Assume layer base is fixed during timestep.

      IF ( zhsc_np1(i,j)  >   z_uv(i,j,k+1) ) THEN
        ! ZHSC risen above level K+1 so specify approp. flux at this level

              ! Adjust the entraiment rate to allow for that implied by
              ! the subsidence increment applied below (combine  increments
              ! from levels NTDSC and NTDSC+1)
        tls_inv(i,j) = tls_inc(i,j,k) + tls_inc(i,j,k-1)
        w_s_ent = 0.0
        IF ( dtl_dsc(i,j)  /=  0.0 ) w_s_ent =                          &
                 MIN( 0.0, -tls_inv(i,j) * dzl(i,j,k) /dtl_dsc(i,j) )
         ! Only allow w_e to be reduced to zero!
        we_lim_dsc(i,j,3) = rho_mix(i,j,k+1) *                          &
                        MAX( 0.0, we_dsc_parm(i,j) + w_s_ent )

        tothf_inv = - we_lim_dsc(i,j,3)*dtl_dsc(i,j) +                  &
                                                      df_net_dsc(i,j)
        zrzi_dsc(i,j,3) = ( z_uv(i,j,k+1)-                              &
                                        (zhsc(i,j)-dscdepth(i,j)) )/    &
          ( 0.5*(zhsc_np1(i,j)+z_uv(i,j,k+1))-                          &
                                        (zhsc(i,j)-dscdepth(i,j)) )
        t_frac_dsc(i,j,3) = (zhsc_np1(i,j)-z_uv(i,j,k+1)) /             &
                                          (zhsc_np1(i,j)-zhsc(i,j))
        ftl(i,j,k+1) = t_frac_dsc(i,j,3) * (                            &
                   tothf_inv * zrzi_dsc(i,j,3) - f_net_dsc(i,j,k+1) )
        rhokh_top(i,j,k+1) = 0.0   ! apply entrainment explicitly
        rhokh(i,j,k+1) = 0.0       !      "

        ! Layer should now be well-mixed to NTDSC+1 (by end of timestep)
        ! so raise NTDSC by one:
        ntdsc(i,j) = ntdsc(i,j) + 1
        ! Note KH profiles should already be calculated at level
        ! NTDSC_old+1 because ZHSC is above this level.

      ELSE  ! ZHSC always below half-level K+1

        rhokh(i,j,k) = 0.0       ! apply entrainment explicitly
        rhokh_top(i,j,k) = 0.0   !      "
        IF ( zhsc_np1(i,j)  >=  z_uv(i,j,k) ) THEN
          !..ZH always above level K so specify full entrainment flux here
          w_s_ent = 0.0
          IF ( dtl_dsc(i,j)  /=  0.0 ) w_s_ent =                        &
            MIN( 0.0, -tls_inc(i,j,k-1) * dzl(i,j,k-1) /dtl_dsc(i,j) )
          we_lim_dsc(i,j,2) = rho_mix(i,j,k) *                          &
                          MAX( 0.0, we_dsc_parm(i,j) + w_s_ent )
          tothf_inv = - we_lim_dsc(i,j,2)*dtl_dsc(i,j) +                &
                                                     df_net_dsc(i,j)
          zrzi_dsc(i,j,2) =( z_uv(i,j,k)-(zhsc(i,j)-dscdepth(i,j)) )    &
                      /( dscdepth(i,j)+0.5*(zhsc_np1(i,j)-zhsc(i,j)) )
          ! ZRZI_DSC = (z - z_base) / (zhsc - z_base)
          t_frac_dsc(i,j,2) = 1.0

          ftl(i,j,k) = tothf_inv * zrzi_dsc(i,j,2) - f_net_dsc(i,j,k)

        ELSE
          !..ZH dropped below level K so specify full entrainment flux at K-1...

          ntdsc(i,j) = ntdsc(i,j) - 1  ! could reduce NTDSC to 1
          w_s_ent = 0.0
          IF ( dtl_dsc(i,j)  /=  0.0 ) w_s_ent = MIN( 0.0,              &
                     -tls_inc(i,j,k-2) * dzl(i,j,k-2) /dtl_dsc(i,j) )
          we_lim_dsc(i,j,1) = rho_mix(i,j,k-1) *                        &
                             MAX( 0.0, we_dsc_parm(i,j) + w_s_ent )

          tothf_inv = - we_lim_dsc(i,j,1)*dtl_dsc(i,j) +                &
                                                       df_net_dsc(i,j)
          zrzi_dsc(i,j,1) =                                             &
                   ( z_uv(i,j,k-1)-(zhsc(i,j)-dscdepth(i,j)) )          &
                      /( dscdepth(i,j)+0.5*(zhsc_np1(i,j)-zhsc(i,j)) )
          t_frac_dsc(i,j,1) = 1.0
          ftl(i,j,k-1) = tothf_inv*zrzi_dsc(i,j,1) -f_net_dsc(i,j,k-1)

          rhokh_top(i,j,k-1) = 0.0   ! apply entrainment explicitly
          rhokh(i,j,k-1) = 0.0       !      "
          ! ...and set specified flux at level K for time spent
          ! with ZH above this level
          w_s_ent = 0.0
          IF ( dtl_dsc(i,j)  /=  0.0 ) w_s_ent = MIN( 0.0,              &
                      -tls_inc(i,j,k-1) * dzl(i,j,k-1) /dtl_dsc(i,j) )
          we_lim_dsc(i,j,2) = rho_mix(i,j,k) *                          &
                       MAX( 0.0, we_dsc_parm(i,j) + w_s_ent )
          tothf_inv = - we_lim_dsc(i,j,2)*dtl_dsc(i,j) +                &
                                                       df_net_dsc(i,j)
          t_frac_dsc(i,j,2) = (zhsc(i,j)-z_uv(i,j,k)) /                 &
                                      (zhsc(i,j)-zhsc_np1(i,j))
          zrzi_dsc(i,j,2) =( z_uv(i,j,k)-(zhsc(i,j)-dscdepth(i,j)) )    &
                      /( dscdepth(i,j)+0.5*(zhsc(i,j)-z_uv(i,j,k)) )
          ftl(i,j,k) = t_frac_dsc(i,j,2) * (                            &
                      tothf_inv * zrzi_dsc(i,j,2) - f_net_dsc(i,j,k) )
        END IF

      END IF  ! ZHSC gets above or stays below half-level K+1

    ELSE IF ( dsc(i,j) ) THEN

        ! not specifying entrainment flux but KH
      rhokh_top(i,j,k) = rhokh_dsct_ent(i,j)

    END IF  ! IF NOT DSC

  END DO
END DO
!-----------------------------------------------------------------------
! Specify QW entrainment fluxes
! (must be separate loops to avoid overwriting)
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    DO ient = 1, 3
      IF ( kent(i,j)-2+ient  >   1 )                                    &
        fqw(i,j,kent(i,j)-2+ient) = t_frac(i,j,ient)*( fqw(i,j,1)       &
                 - ( we_lim(i,j,ient)*dqw_sml(i,j) + fqw(i,j,1) )       &
                                     *zrzi(i,j,ient) )

    END DO
    DO ient = 1, 3
      IF ( kent_dsc(i,j) >= 3 .AND. t_frac_dsc(i,j,ient) >  0.0 )       &
        fqw(i,j,kent_dsc(i,j)-2+ient) =  - t_frac_dsc(i,j,ient) *       &
            we_lim_dsc(i,j,ient) * dqw_dsc(i,j) * zrzi_dsc(i,j,ient)
    END DO
  END DO
END DO
!-----------------------------------------------------------------------
!- Update standard deviations and gradient adjustment to use this
!- timestep's ZH (code from SF_EXCH)
!-----------------------------------------------------------------------
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    w_s_cubed = 0.25 * zh(i,j) * fb_surf(i,j)
    IF (w_s_cubed  >   0.0) THEN
      w_m  =                                                            &
     ( w_s_cubed + v_s(i,j) * v_s(i,j) * v_s(i,j) ) ** one_third
      t1_sd(i,j) = 1.93 * ftl(i,j,1) / (rhostar_gb(i,j) * w_m)
      q1_sd(i,j) = 1.93 * fqw(i,j,1) / (rhostar_gb(i,j) * w_m)
      tv1_sd(i,j) = t(i,j,1) *                                          &
          ( 1.0 + c_virtual*q(i,j,1) - qcl(i,j,1) - qcf(i,j,1) ) *      &
          ( bt(i,j,1)*t1_sd(i,j) + bq(i,j,1)*q1_sd(i,j) )
      t1_sd(i,j) = MAX ( 0.0 , t1_sd(i,j) )
      q1_sd(i,j) = MAX ( 0.0 , q1_sd(i,j) )
      IF (tv1_sd(i,j)  <=  0.0) THEN
        tv1_sd(i,j) = 0.0
        t1_sd(i,j) = 0.0
        q1_sd(i,j) = 0.0
      END IF
    END IF
    grad_t_adj(i,j) = MIN( max_t_grad ,                                 &
                         a_grad_adj * t1_sd(i,j) / zh(i,j) )
    !        IF (T1_SD(I,j)  >   0.0) THEN
    !          GRAD_Q_ADJ(I,j) = (Q1_SD(I,j) / T1_SD(I,j)) * GRAD_T_ADJ(I,j)
    !        ELSE
    grad_q_adj(i,j) = 0.0
    !        END IF
  END DO
END DO
!-----------------------------------------------------------------------
!- Save diagnostics
!-----------------------------------------------------------------------
IF (BL_diag%l_dscbase) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( dsc(i,j) ) THEN
        BL_diag%dscbase(i,j)= zhsc(i,j)-dscdepth(i,j)
      ELSE
        BL_diag%dscbase(i,j)= rmdi
      END IF
    END DO
  END DO
END IF
IF (BL_diag%l_cldbase) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( dsc(i,j) ) THEN
        BL_diag%cldbase(i,j)= zhsc(i,j)-zc_dsc(i,j)
      ELSE
        BL_diag%cldbase(i,j)= zh(i,j)-zc(i,j)
      END IF
    END DO
  END DO
END IF
IF (BL_diag%l_weparm_dsc) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      IF ( dsc(i,j) ) THEN
        BL_diag%weparm_dsc(i,j)= we_dsc_parm(i,j)
      ELSE
        BL_diag%weparm_dsc(i,j)= we_parm(i,j)
      END IF
    END DO
  END DO
END IF
IF (BL_diag%l_weparm) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      BL_diag%weparm(i,j)= we_parm(i,j)
    END DO
  END DO
END IF

!-----------------------------------------------------------------------
!     SCM Boundary Layer Diagnostics Package
!-----------------------------------------------------------------------
IF ( l_scmdiags(scmdiag_bl) .AND.                                       &
     model_type == mt_single_column ) THEN

  !------------------------------------------------------
  DO j=pdims%j_start, pdims%j_end
    jScm = j - pdims%j_start + 1
    DO i=pdims%i_start, pdims%i_end
      iScm = i - pdims%i_start + 1
      IF ( dsc(i,j) ) THEN
        TmpScm2d(iScm,jScm)= zhsc(i,j) - dscdepth(i,j)
      ELSE
        TmpScm2d(iScm,jScm)= zh(i,j)
      END IF
    END DO
  END DO
  CALL scmoutput(TmpScm2d,'DSCbase',                                    &
       'Base of decoupled layer','m',                                   &
       t_avg,d_sl,default_streams,'',routinename)
  !------------------------------------------------------
  DO j=pdims%j_start, pdims%j_end
    jScm = j - pdims%j_start + 1
    DO i=pdims%i_start, pdims%i_end
      iScm = i - pdims%i_start + 1
      IF ( dsc(i,j) ) THEN
        TmpScm2d(iScm,jScm) = zhsc(i,j) - zc_dsc(i,j)
      ELSE
        TmpScm2d(iScm,jScm) = zh(i,j)   - zc(i,j)
      END IF
    END DO
  END DO
  CALL scmoutput(TmpScm2d,'ScCldBase',                                  &
       'Stratocumulus cloud base','m',                                  &
       t_avg,d_sl,default_streams,'',routinename)
  !------------------------------------------------------
  CALL scmoutput(we_parm,'Entr_SML',                                    &
       'SML-top entrainment rate','m/s',                                &
       t_avg,d_sl,default_streams,'',routinename)
  !------------------------------------------------------
  DO j=pdims%j_start, pdims%j_end
    jScm = j - pdims%j_start + 1
    DO i=pdims%i_start, pdims%i_end
      iScm = i - pdims%i_start + 1
      IF ( dsc(i,j) ) THEN
        TmpScm2d(iScm,jScm) = we_dsc_parm(i,j)
      ELSE
        TmpScm2d(iScm,jScm) = we_parm(i,j)
      END IF
    END DO
  END DO
  CALL scmoutput(TmpScm2d,'Entr_BL',                                    &
       'BL-top entrainment rate','m/s',                                 &
       t_avg,d_sl,default_streams,'',routinename)
  !------------------------------------------------------

END IF ! scmdiag_bl / model_type


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE kmkhz_9b
END MODULE kmkhz_9b_mod
