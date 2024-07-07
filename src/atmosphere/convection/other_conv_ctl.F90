! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Control routine for calling other alternative convections scheme

MODULE other_conv_ctl_mod

IMPLICIT NONE
!---------------------------------------------------------------------------
! Description: Called from atmos_physic2 in place of ni_conv_ctl
!
!  List of alternative convection schemes: with their i_convection_vn
!  number and code name :
!  
!   Number   :  Scheme                   : code name
!            :                           :
!    10      :  Lambert-Lewis moist-     : llcs
!            :  convective adjustment    :
!            :  scheme.                  :
!            :                           :
!
! method: Passes in all prognostics a convection scheme may need
!         plus some other fields from boundary layer
!         *_n arrays hold the start of time step values
!         *_inc arrays hold the convection increments
!         *_star arrays hold start-of-timestep plus all increments so far
!         Normally the convection scheme has worked off the incoming *_star
!         arrays passing back updated values.
!         Note U and V winds are interpolated in the horizontal by 
!         atmos-physics2 but are are still on rho levels in the vertical.

!         Note U and V winds are only put on the p grid if l_mom = .true.
!         i.e. the convection scheme is expected to produce CMT increments. 
!         U & V tendencies from convection are stored on the p-grid here,
!         in the output arrays dubydt_pout & dvbydt_pout.
!         Interpolation back to the U and V grids is done towards the end of
!         atmos_physics2.

! Moist variables - these may be specific quantities or mixing ratios
! dependent on the setting of l_mr_physics

! Section 5 diagnostics
!         Most convective diagnostics are held in cv_diagnostic_array_mod
!         The convective cloud and surface precipitation diagnostics pass 
!         through the arguement list to this routine.
!         Space for some is allocated by  cv_alloc_diag_array called from 
!         atmos_physics2 before this routine dependent on stashflag settings.  
!         Any new convection scheme must ensure it sets diagnostics from
!         within the scheme by using the module.
!         The call to diagnostics_conv towards the end of atmos_physics2 
!         outputs diagnostics to stash.
!         Note that convective increment diagnostics are set in the call to
!         save_conv_diags.
!         If you want to output any convection diagnostic to STASH, you
!         need to ensure your alternative convection scheme's version number
!         has the version mask set in the STASHmaster entry for that
!         diagnostic.  This includes the convection increment diagnostics.

!  Convection time stepping - this routine assumes no sub-stepping
!  Sub-stepping is best avoided.       



! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection

! code description:
!   language: fortran 90
!   This code is written to umdp3 programming standards version 10.4
!---------------------------------------------------------------------------

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='OTHER_CONV_CTL_MOD'

CONTAINS

SUBROUTINE  other_conv_ctl(                                                  &
                numcycles, cycleno,                                          &
                row_length, rows, n_cca_levels,                              &
                tr_levels, tr_vars, tr_ukca,                                 &
                land_points, nscmdpkgs,                                      &
                ntml,                                                        &
             ! Logical in
                l_calc_dxek, l_q_interact, l_scmdiags,                       &
                land_sea_mask,                                               &
             ! Real in...
                zh, dzh, taux_p, tauy_p, ftl, fqw, rhokm, rhokh, bl_w_var,   &
                rho_wet_rsq, rho_wet, rho_wet_th, rho_dry, rho_dry_th,       &
                p_star, p, exner_rho_levels, exner_theta_levels,             &
                p_layer_boundaries, p_layer_centres,                         &
                exner_layer_boundaries, exner_layer_centres,                 &
             ! in Prognostics at beginning of time step
                theta_n, q_n, qcl_n, qcf_n, qrain_n, qgraup_n, qcf2_n,       &
                cf_liquid_n, cf_frozen_n, bulk_cf_n,                         &
                u_p, v_p, w,                                                 &
                m_v, m_cl, m_cf, m_cf2, m_r, m_gr,                           &
             ! Integer inout (cloud level info) 
                ccb0, cct0, lcbase0, ccb, cct, lcbase, lctop,                &
             ! In/Out: fields updated with all increments so far
                theta_star, q_star, qcl_star, qcf_star,                      & 
                qcf2_star, qrain_star, qgraup_star,                          &
                cf_liquid_star, cf_frozen_star, bulk_cf_star,                &
             ! In: winds updated with increments so far, interpolated to p-grid
                ustar_p, vstar_p,                                            &
             ! Out: convective momentum transport tendencies on p-grid
                dubydt_pout, dvbydt_pout,                                    &
             ! Cloud prognostics & diagostics
                cca0_dp, cca0_md, cca0_sh, cca0, ccw0, cclwp0, cca0_2d,      &
                lcca, cca,  ccw, cclwp, cca_2d,                              &
             ! Other output fields
                conv_rain, conv_snow,                                        &
             ! Tracers
                aerosol, dust_div1, dust_div2,                               &
                dust_div3, dust_div4, dust_div5, dust_div6,                  &
                so2, so4_aitken, so4_accu, so4_diss,                         &
                dms, nh3, soot_new, soot_aged, soot_cld,                     &
                bmass_new, bmass_aged, bmass_cld,                            &
                ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss,          &
                co2,  free_tracers, tracer_ukca, ozone_tracer,               &
             ! Real out - may be expected by surface scheme
                ddmfx                                                        &
                )


!---------------------------------------------------------------------------
! Modules holding information which may be required.
!---------------------------------------------------------------------------
! Required for any openmp
!$ USE omp_lib

USE gen_phys_inputs_mod, ONLY: l_mr_physics

! Definitions of prognostic variable array sizes
USE atm_fields_bounds_mod, ONLY:                                        &
   tdims, pdims, tdims_s, pdims_s, wdims_s,                             &
   tdims_l

USE nlsizes_namelist_mod,   ONLY: model_levels, bl_levels

USE cloud_inputs_mod,       ONLY: i_cld_vn ! used to decide how rain and
                                           ! cloud are incremented
USE pc2_constants_mod,      ONLY: i_cld_pc2
! Convection scheme number
USE cv_run_mod,  ONLY:                                                  &
    i_convection_vn, l_mom, i_cv_llcs,                                  &
    l_mcr_conv

USE llcs, ONLY: llcs_control, l_output_is_cloud
USE dynamics_input_mod, ONLY: l_check_moist_inc
USE check_dmoist_inc_mod, ONLY: check_dmoist_inc
USE cv_stash_flg_mod, ONLY: l_apply_diag
USE model_domain_mod, ONLY: model_type, mt_global, mt_single_column

! Extra moist variable switches
USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup

! Tracer switches

USE run_aerosol_mod, ONLY: l_sulpc_so2, l_sulpc_dms, l_sulpc_nh3,       &
                           l_soot, l_ocff, l_biomass, l_nitrate
USE carbon_options_mod, ONLY: l_co2_interactive
USE dust_parameters_mod, ONLY: l_dust

! Convective diagnostic output arrays
USE cv_diagnostic_array_mod, ONLY: dwn_flux, ep1, ep2, ep3

! Flags controlling convective diagnostic output
!USE cv_stash_flg_mod, ONLY:                                            &

! Add convection tendencies to physics_tendencies module for SPT scheme
USE physics_tendencies_mod,  ONLY:                                     &
    l_retain_conv_tendencies, dt_conv, dq_conv

! Control request for information to be passed to surface exchange
USE jules_surface_mod, ONLY: ISrfExCnvGust, IP_SrfExWithCnv

!---------------------------------------------------------------------------
! Note information a scheme may want - a lot of information is now in modules
! so does not need to be passed down in arguement lists, e.g.
!
! Model time stepping information
USE timestep_mod, ONLY: timestep, recip_timestep
USE cderived_mod, ONLY: delta_lambda, delta_phi
!
! Model level heights from centre of Earth
! USE level_heights_mod, ONLY: &
!   r_theta_levels             &  ! Radii on theta levels (m)
!  ,r_rho_levels                  ! Radii on rho levels (m)
!---------------------------------------------------------------------------
! Subroutines
USE save_conv_diags_mod,   ONLY: save_conv_diags
USE pc2_from_conv_ctl_mod, ONLY: pc2_from_conv_ctl

!---------------------------------------------------------------------------
! Required by most routines - message printing and error reporting
USE umPrintMgr
USE errormessagelength_mod, ONLY: errormessagelength
USE ereport_mod, ONLY: ereport
!---------------------------------------------------------------------------

IMPLICIT NONE

!----------------------------------------------------------------------
! Subroutine arguments

INTEGER, INTENT(IN)  ::      &
  numcycles                  & ! Number of cycles
, cycleno                      ! cycle number

! Model dimensions
INTEGER, INTENT(IN)  ::      &
  row_length                 & ! Row length
, rows                       & ! Number of rows
, n_cca_levels               & ! Number of levels for conv cloud
                               ! amount: 1 for 2D, nlevs for 3D.
, tr_levels                  & ! free tracer levels
, tr_vars                    & ! number of free tracers
, tr_ukca                    & ! number of ukca tracers
, land_points                  ! No.of land points being processed, can be 0.


! Additional variables for SCM diagnostics which are dummy in full UM
INTEGER, INTENT(IN)  ::  &
  nscmdpkgs                ! No of SCM diagnostics packages

INTEGER, INTENT(IN)  ::  &
  ntml(row_length, rows)   ! Top level of surface mixed layer used by BL
                           ! Comes from standard conv_diag

LOGICAL, INTENT(IN)  ::            &
  l_calc_dxek                      & ! Switch for PC2 set by atmos_physics2
 ,l_q_interact                       ! Switch for PC2 set by atmos_physics2

LOGICAL, INTENT(IN)  ::            &
  l_scmdiags(nscmdpkgs)            & ! Logicals for SCM diagnostics packages
 ,land_sea_mask(row_length, rows)    ! land sea land

! Data fields 2d fields 

! From the boundary layer called just before convection
REAL, INTENT(IN) :: zh       ( pdims%i_start:pdims%i_end,           &
                               pdims%j_start:pdims%j_end )
                    ! boundary layer height (m)
REAL, INTENT(IN) :: dzh      ( pdims%i_start:pdims%i_end,           &
                               pdims%j_start:pdims%j_end )
                    ! inverison thickness (m)
REAL, INTENT(IN) :: taux_p   ( pdims%i_start:pdims%i_end,           &
                               pdims%j_start:pdims%j_end,           &
                               0:bl_levels-1 )
                    ! X-comp turbulent stress on p grid (on theta-levels) N/m2
REAL, INTENT(IN) :: tauy_p   ( pdims%i_start:pdims%i_end,           &
                               pdims%j_start:pdims%j_end,           &
                               0:bl_levels-1 )
                    ! Y-comp turbulent stress on p grid (on theta-levels) N/m2
! Note: in order for the above arrays of turbulent stresses on the p-grid
! to actually be set on input to this routine, the flag
! l_calc_tau_at_p needs to be set to true in atmos_physics2_alloc
REAL, INTENT(IN) :: ftl      ( pdims%i_start:pdims%i_end,           &
                               pdims%j_start:pdims%j_end,           &
                               bl_levels )
                    ! Sensible heat flux from BL, on rho-levels
                    ! (W/m2) i.e. cp*rho*w'tl'
REAL, INTENT(IN) :: fqw      ( pdims%i_start:pdims%i_end,           &
                               pdims%j_start:pdims%j_end,           &
                               bl_levels )
                    ! Total water flux from BL, on rho-levels
                    ! (kg/m2/s) i.e. rho*w'qT'
REAL, INTENT(IN) :: rhokm    ( pdims_s%i_start:pdims_s%i_end,       &
                               pdims_s%j_start:pdims_s%j_end,       &
                               0:bl_levels-1 )
                    ! rho * turbulent diffusivity for momentum, on theta-grid.
REAL, INTENT(IN) :: rhokh    ( pdims%i_start:pdims%i_end,           &
                               pdims%j_start:pdims%j_end,           &
                               bl_levels )
                    ! rho * turbulent diffusivity for scalars, on rho-levels.

REAL, INTENT(IN) :: bl_w_var ( tdims%i_start:tdims%i_end,           &
                               tdims%j_start:tdims%j_end,           &
                               1:tdims%k_end )
                    ! Sub-grid turbulent variance in vertical velocity

! Note: some of the above turbulence fields follow a slightly different
! vertical level indexing convention to the rest of the model:
!
! The fields defined on rho-levels (ftl, fqw, rhokh) have k=1 defined to
! be the surface, and don't have a defined value on the actual first rho-level.
!
! The fields defined on theta-levels (taux_p, tauy_p, rhokm, bl_w_var) have
! k=0 defined to be the surface, like theta (bl_w_var not defined at surface).
! This is consistent with how they are indexed in atmos_physics2
! (within the boundary-layer scheme taux, tauy, rhokm are declared and
!  indexed differently, with k=1 at the surface and rhokm(k) at theta(k-1) ).

! Density information (note prognostic is held on the rho grid)
REAL, INTENT(IN) ::                                    &
  rho_wet_rsq(pdims_s%i_start:pdims_s%i_end,           & ! density *r*r (kg/m)
              pdims_s%j_start:pdims_s%j_end,           & !
              pdims_s%k_start:pdims_s%k_end)           & !
, rho_wet(pdims%i_end,pdims%j_end,pdims%k_end)         & ! wet density (kg/m3)
, rho_wet_th(tdims%i_end,tdims%j_end,tdims%k_end-1)    & ! wet on th lev (kg/m3)
, rho_dry(pdims%i_end,pdims%j_end,pdims%k_end)         & ! dry density (kg/m3)
, rho_dry_th(tdims%i_end,tdims%j_end,tdims%k_end-1)      ! dry on th lev (kg/m3)

! Pressure/Exner info
REAL, INTENT(IN) :: p_star(row_length, rows)         ! surface pressure (Pa)
REAL, INTENT(IN) ::                                &
  p(pdims_s%i_start:pdims_s%i_end,                 & ! pressure rho levels (Pa)
    pdims_s%j_start:pdims_s%j_end,                 &
    pdims_s%k_start:pdims_s%k_end)                 &
, exner_rho_levels(pdims_s%i_start:pdims_s%i_end,  & ! Exner on rho level
                   pdims_s%j_start:pdims_s%j_end,  & !
                   pdims_s%k_start:pdims_s%k_end)  &
, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,  & ! Exner on theta level
                     tdims_s%j_start:tdims_s%j_end,  & !
                     tdims_s%k_start:tdims_s%k_end)

! Arrays setup in atmos_physics2 for use by BL & conv
REAL, INTENT(IN) ::                          & ! pressure at layer boundaries.
  p_layer_boundaries(pdims%i_end,            & ! Same as p except at bottom
                     pdims%j_end,            & ! level = pstar, and at top = 0.
                     0:pdims%k_end)          &
, p_layer_centres(tdims%i_end,               & ! pressure at layer centres.
                  tdims%j_end,               & ! Same as p_theta_levels except
                  0:tdims%k_end)             & ! bottom level = pstar, & top =0.
, exner_layer_boundaries(pdims%i_end,        & ! Exner at p_layer_boundaries
                         pdims%j_end,        &
                         0:pdims%k_end)      &
, exner_layer_centres(tdims%i_end,           & ! Exner at p_layer_centres
                      tdims%j_end,           &
                      0:tdims%k_end)


! Prognostics start of timestep values INTENT(IN) as must not be altered

REAL,INTENT(IN) ::                           &
  theta_n(tdims_s%i_start:tdims_s%i_end,     & ! theta (K)
          tdims_s%j_start:tdims_s%j_end,     &
          tdims_s%k_start:tdims_s%k_end)     &
, q_n(tdims_l%i_start:tdims_l%i_end,         & ! q (kg/kg)
      tdims_l%j_start:tdims_l%j_end,         &
      tdims_l%k_start:tdims_l%k_end)         &
, qcl_n(tdims_l%i_start:tdims_l%i_end,       & ! qcl (kg/kg)
        tdims_l%j_start:tdims_l%j_end,       &
        tdims_l%k_start:tdims_l%k_end)       &
, qcf_n(tdims_l%i_start:tdims_l%i_end,       & ! qcf (kg/kg)
        tdims_l%j_start:tdims_l%j_end,       &
        tdims_l%k_start:tdims_l%k_end)       &
, qrain_n(tdims_l%i_start:tdims_l%i_end,     & ! qrain (kg/kg)
          tdims_l%j_start:tdims_l%j_end,     &
          tdims_l%k_start:tdims_l%k_end)     &
, qgraup_n(tdims_l%i_start:tdims_l%i_end,    & ! qraup (kg/kg)
           tdims_l%j_start:tdims_l%j_end,    &
           tdims_l%k_start:tdims_l%k_end)    &
, qcf2_n(tdims_l%i_start:tdims_l%i_end,      & ! qcf2 (kg/kg)
         tdims_l%j_start:tdims_l%j_end,      &
         tdims_l%k_start:tdims_l%k_end)      &
, cf_liquid_n(tdims_l%i_start:tdims_l%i_end, & ! liquid cloud fraction
              tdims_l%j_start:tdims_l%j_end, &
              tdims_l%k_start:tdims_l%k_end) &
, cf_frozen_n(tdims_l%i_start:tdims_l%i_end, & ! ice cloud fraction
              tdims_l%j_start:tdims_l%j_end, &
              tdims_l%k_start:tdims_l%k_end) &
, bulk_cf_n(tdims_l%i_start:tdims_l%i_end,   & ! total cloud fraction
            tdims_l%j_start:tdims_l%j_end,   &
            tdims_l%k_start:tdims_l%k_end)

! Winds at start of timestep but on p-grid and without halos
REAL, INTENT(IN) ::                          &
  u_p(pdims%i_end,pdims%j_end,pdims%k_end)   & ! U (m/s) P-grid points
, v_p(pdims%i_end,pdims%j_end,pdims%k_end)   & ! V (m/s) P-grid points
, w(wdims_s%i_start:wdims_s%i_end,           & ! W (m/s)
    wdims_s%j_start:wdims_s%j_end,           &
    wdims_s%k_start:wdims_s%k_end)

! Start-of-timestep mixing ratios of water species
! Note: in SCM runs, these are not set unless the flag l_mr_conv in
!       cv_run_mod is set to true.
                    ! Water vapour
REAL, INTENT(IN) :: m_v   ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Cloud liquid water
REAL, INTENT(IN) :: m_cl  ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Cloud ice and snow
REAL, INTENT(IN) :: m_cf  ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! 2nd cloud ice category (optional)
REAL, INTENT(IN) :: m_cf2 ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Prognostic rain water (optional)
REAL, INTENT(IN) :: m_r   ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )
                    ! Prognostic graupel (optional)
REAL, INTENT(IN) :: m_gr  ( tdims_s%i_start:tdims_s%i_end,              &
                            tdims_s%j_start:tdims_s%j_end,              &
                            tdims_s%k_start:tdims_s%k_end )


! Variables with 0 are prognostics passed around to radiation, use
! depends on whether using PC2, convective cores etc. Others diagnostic
INTEGER, INTENT(INOUT) ::    &
  ccb0    (row_length,rows)  &! Cnv Cld Base of highest layer on gridpoint
, cct0    (row_length,rows)  &! Cnv Cld Top  of highest layer on gridpoint
, lcbase0 (row_length,rows)   ! Cnv Cld Top  of lowest  layer on gridpoint
INTEGER, INTENT(INOUT) ::  &!
  ccb    (row_length,rows) &! Cnv Cld Base of highest layer on gridpoint
, cct    (row_length,rows) &! Cnv Cld Top  of highest layer on gridpoint
, lcbase (row_length,rows) &! Cnv Cld Base of lowest  layer on gridpoint
, lctop  (row_length,rows)  ! Cnv Cld Top  of lowest  layer on gridpoint

! '_star' = '_n' + all increments up to now in time step
! _star variables held on same levels as prognostics in the vertical for
! ENDGame. Note Convection will not update any variables stored on level 0,
! so the zero-level has not been passed in (by subscripting the arrays in
! the argument list to this routine in atmos_physics2).

REAL,INTENT(INOUT) ::                          &
  theta_star(tdims%i_start:tdims%i_end,        & ! Potential temperature (K)
             tdims%j_start:tdims%j_end,        &
                         1:tdims%k_end)        &
, q_star(tdims%i_start:tdims%i_end,            & ! Water vapour (kg/kg)
         tdims%j_start:tdims%j_end,            &
                     1:tdims%k_end)            &
, qcl_star(tdims%i_start:tdims%i_end,          & ! cloud liquid water (kg/kg)
           tdims%j_start:tdims%j_end,          &
                       1:tdims%k_end)          &
, qcf_star(tdims%i_start:tdims%i_end,          & ! cloud ice water (kg/kg)
           tdims%j_start:tdims%j_end,          &
                       1:tdims%k_end)          &
, qcf2_star(tdims%i_start:tdims%i_end,         & ! 2nd ice cloud type (kg/kg)
            tdims%j_start:tdims%j_end,         &
                        1:tdims%k_end)         &
, qrain_star(tdims%i_start:tdims%i_end,        & ! rain water (kg/kg)
             tdims%j_start:tdims%j_end,        &
                         1:tdims%k_end)        &
, qgraup_star(tdims%i_start:tdims%i_end,       & ! graupel (kg/kg)
              tdims%j_start:tdims%j_end,       &
                          1:tdims%k_end)       &
, cf_liquid_star(tdims%i_start:tdims%i_end,    & ! Fraction of liquid cloud
                 tdims%j_start:tdims%j_end,    &
                             1:tdims%k_end)    &
, cf_frozen_star(tdims%i_start:tdims%i_end,    & ! Fraction of ice cloud
                 tdims%j_start:tdims%j_end,    &
                             1:tdims%k_end)    &
, bulk_cf_star(tdims%i_start:tdims%i_end,      & ! Fraction of cloud
               tdims%j_start:tdims%j_end,      &
                           1:tdims%k_end)
! u and v updated with all increments so far and interpolated to p-grid.
REAL,INTENT(IN) ::                          &
  ustar_p(pdims%i_start:pdims%i_end,        & ! U+du (m/s) P-grid
          pdims%j_start:pdims%j_end,        &
          pdims%k_start:pdims%k_end)        &
, vstar_p(pdims%i_start:pdims%i_end,        & ! V+dv (m/s) P-grid
          pdims%j_start:pdims%j_end,        &
          pdims%k_start:pdims%k_end)
! Note: cannot update the primary fields by modifying the above; u,v are
! updated with convection increments in atmos_physics2, after the
! convective tendencies have been interpolated onto u/v points.

! Out convective momentum transport tendencies on p-grid
REAL, INTENT(OUT) :: dubydt_pout ( pdims_s%i_start:pdims_s%i_end,        &
                                   pdims_s%j_start:pdims_s%j_end,        &
                                   pdims_s%k_start:pdims_s%k_end )
REAL, INTENT(OUT) :: dvbydt_pout ( pdims_s%i_start:pdims_s%i_end,        &
                                   pdims_s%j_start:pdims_s%j_end,        &
                                   pdims_s%k_start:pdims_s%k_end )

! Convective cloud prognostics with 0 in name, those actual in use depend
! on PC2, and various settings controlling whether 3d convective cloud
REAL, INTENT(INOUT) ::                    &
  cca0    (row_length,rows,n_cca_levels)  &! Cnv.Cld Amount (0-1)
, cca0_dp (row_length,rows,n_cca_levels)  &! Deep cnv.Cld Amount (0-1)
, cca0_md (row_length,rows,n_cca_levels)  &! Mid-level cnv.Cld Amount (0-1)
, cca0_sh (row_length,rows,n_cca_levels)  &! Shallow cnv.Cld Amount (0-1)
, ccw0    (tdims%i_end,tdims%j_end,       &! Cnv.Cld Water (kg/kg)
           tdims%k_end)                   &! 
, cclwp0  (row_length,rows)               &! Cloud Condensed water path
, cca0_2d (row_length,rows)                ! Cnv.Cld Amount (2d) with no anvil
                                           ! Only used if: l_3d_cca = .False. 
                                           ! .AND. l_pc2 .AND. l_pc2_diag_sh
! Diagnostics
REAL, INTENT(INOUT) ::                    &
  cca    (row_length,rows,n_cca_levels)   &! Cnv.Cld Amount (0-1)
, ccw    (tdims%i_end,tdims%j_end,        &! Cnv.Cld Water (kg/kg)
          tdims%k_end)                    &! 
, cclwp  (row_length,rows)                &! Cloud Condensed water path
                                           ! (kg/m^2)
, cca_2d(row_length,rows)                 &! Cnv.Cld Amount (2d)
                                           ! with no anvil
, lcca   (row_length,rows)                 ! Conv. Cloud Amount of
                                           ! lowest Convective Cloud Layer

REAL, INTENT(OUT) ::                  &
  conv_rain(row_length, rows)         & ! Convective surface rainfall (kg/m2/s)
, conv_snow(row_length, rows)           ! Convective surface snowfall (kg/m2/s)


! Tracers Prognostics - these are time level n plus all increments up to this
! point.

REAL,INTENT(INOUT) ::                         &
  aerosol(tdims_s%i_start:tdims_s%i_end,      & ! aerosol tracer
          tdims_s%j_start:tdims_s%j_end,      & ! (called "murk" in atm_step).
          tdims_s%k_start:tdims_s%k_end)      &
, dust_div1(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div1
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div2(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div2
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div3(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div3
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div4(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div4
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div5(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div5
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, dust_div6(tdims_s%i_start:tdims_s%i_end,    & !dust mmr in div6
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                         &
  so2(tdims_s%i_start:tdims_s%i_end,          & ! SO2
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, so4_aitken(tdims_s%i_start:tdims_s%i_end,   & ! SO4 aitken mode
             tdims_s%j_start:tdims_s%j_end,   &
             tdims_s%k_start:tdims_s%k_end)   &
, so4_accu(  tdims_s%i_start:tdims_s%i_end,   & ! SO4 accumulation mode
             tdims_s%j_start:tdims_s%j_end,   &
             tdims_s%k_start:tdims_s%k_end)   &
, so4_diss(tdims_s%i_start:tdims_s%i_end,     & ! SO4 dissipation mode
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, dms(tdims_s%i_start:tdims_s%i_end,          & ! DMS
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, nh3(tdims_s%i_start:tdims_s%i_end,          & ! NH3
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)          &
, soot_new(tdims_s%i_start:tdims_s%i_end,     & ! New soot
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, soot_aged(tdims_s%i_start:tdims_s%i_end,    & ! Aged soot
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, soot_cld(tdims_s%i_start:tdims_s%i_end,     & ! soot in cloud
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, bmass_new(tdims_s%i_start:tdims_s%i_end,    & ! New biomass
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, bmass_aged(tdims_s%i_start:tdims_s%i_end,   & ! Aged biomass
             tdims_s%j_start:tdims_s%j_end,   &
             tdims_s%k_start:tdims_s%k_end)   &
, bmass_cld(tdims_s%i_start:tdims_s%i_end,    & ! Biomass in cloud
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                         &
  ocff_new(tdims_s%i_start:tdims_s%i_end,     & ! New ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, ocff_aged(tdims_s%i_start:tdims_s%i_end,    & ! Aged ocff
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, ocff_cld(tdims_s%i_start:tdims_s%i_end,     & ! Cloud ocff
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_acc(tdims_s%i_start:tdims_s%i_end,     & ! nitrate accumulation mode
           tdims_s%j_start:tdims_s%j_end,     &
           tdims_s%k_start:tdims_s%k_end)     &
, nitr_diss(tdims_s%i_start:tdims_s%i_end,    & ! nitrate dissipation mode
            tdims_s%j_start:tdims_s%j_end,    &
            tdims_s%k_start:tdims_s%k_end)    &
, co2(tdims_s%i_start:tdims_s%i_end,          & ! CO2
      tdims_s%j_start:tdims_s%j_end,          &
      tdims_s%k_start:tdims_s%k_end)

REAL,INTENT(INOUT) ::                                   &
  free_tracers(tdims_s%i_start:tdims_s%i_end,           & ! free tracers
               tdims_s%j_start:tdims_s%j_end,           &
               tdims_s%k_start:tdims_s%k_end,           &
               tr_vars)                                 &
, tracer_ukca(tdims_s%i_start:tdims_s%i_end,            & ! UKCA tracers
              tdims_s%j_start:tdims_s%j_end,            &
              tdims_s%k_start:tdims_s%k_end,            &
              tr_ukca)                                  &
, ozone_tracer( tdims_s%i_start:tdims_s%i_end,          & ! Ozone tracer
                tdims_s%j_start:tdims_s%j_end,          &
                tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(OUT) ::             &
  ddmfx(row_length, rows)          ! Convective downdraught mass-flux at 
                                   ! cloud-base, store if required
                                   ! for surface exchange.

!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------

! '_inc' - convective increments to a field
! Note Convection not expected to update or alter level zero, so the
! increment arrays are declared without it.
REAL :: theta_inc ( tdims%i_start:tdims%i_end,          & ! Potential temp.
                    tdims%j_start:tdims%j_end,          &
                    1:tdims%k_end )
REAL :: q_inc     ( tdims%i_start:tdims%i_end,          & ! Water vapour
                    tdims%j_start:tdims%j_end,          &
                    1:tdims%k_end )
REAL :: qcl_inc   ( tdims%i_start:tdims%i_end,          & ! Cloud liquid water
                    tdims%j_start:tdims%j_end,          &
                    1:tdims%k_end )
REAL :: qcf_inc   ( tdims%i_start:tdims%i_end,          & ! Cloud ice
                    tdims%j_start:tdims%j_end,          &
                    1:tdims%k_end )
! Only need increments to extra microphysics fields if they are used
REAL, ALLOCATABLE :: qcf2_inc(:,:,:)                  ! 2nd cloud ice category
REAL, ALLOCATABLE :: qrain_inc(:,:,:)                 ! Prognostic rain
REAL, ALLOCATABLE :: qgraup_inc(:,:,:)                ! Prognostic graupel
! Only need increments to prognostic cloud fields if PC2 is on.
REAL, ALLOCATABLE :: cf_liquid_inc(:,:,:)             ! Liquid cloud fraction
REAL, ALLOCATABLE :: cf_frozen_inc(:,:,:)             ! Ice cloud fraction
REAL, ALLOCATABLE :: bulk_cf_inc(:,:,:)               ! Total cloud fraction

INTEGER ::             &
  i, j, k              & ! Loop counters
 ,n_conv_levels          ! Number of model levels over which convection 
                         ! can act. 

REAL :: phalfllcs(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,  &
                  1:pdims%k_end+1)
REAL :: pfullllcs(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,  &
                  1:tdims%k_end)
REAL, ALLOCATABLE :: dqt(:,:,:)

CHARACTER (LEN=80) :: scheme_name

! Error reporting variables
INTEGER ::                                                        &
  errorstatus

CHARACTER (LEN=errormessagelength) ::                             &
  cmessage

CHARACTER (LEN=*), PARAMETER ::                                   &
  RoutineName='OTHER_CONV_CTL'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!----------------------------------------------------------------------

!----------------------------------------------------------------------
! Any code required before calling all schemes
!----------------------------------------------------------------------

n_conv_levels = tdims%k_end

IF (l_mom) THEN

  ! Limit convection calling levels to maximum of model_levels - 1
  ! This is because CMT increments to u and v exist on n_levels + 1
  ! Note this setting of n_conv_levels is appropriate for the mass flux
  ! scheme call by ni_conv_ctl and the fact the winds are on rho levels.

  IF (n_conv_levels  >   model_levels - 1 ) THEN
    n_conv_levels = model_levels - 1
  END IF
END IF


! Initialise convection increments to zero...

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      theta_inc(i,j,k) = 0.0
      q_inc(i,j,k) = 0.0
      qcl_inc(i,j,k) = 0.0
      qcf_inc(i,j,k) = 0.0
    END DO
  END DO
END DO
!$OMP END PARALLEL DO

! If extra microphysics prognostic hydrometeors are affected by convection,
! allocate and set increments to zero
IF ( l_mcr_conv ) THEN

  ! 2nd cloud ice category
  IF ( l_mcr_qcf2 ) THEN
    ALLOCATE( qcf2_inc ( tdims%i_start:tdims%i_end,               &
                         tdims%j_start:tdims%j_end,               &
                         1:tdims%k_end ) )
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qcf2_inc(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Prognostic rain
  IF ( l_mcr_qrain ) THEN
    ALLOCATE( qrain_inc ( tdims%i_start:tdims%i_end,              &
                          tdims%j_start:tdims%j_end,              &
                          1:tdims%k_end ) )
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qrain_inc(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Prognostic graupel
  IF ( l_mcr_qgraup ) THEN
    ALLOCATE( qgraup_inc ( tdims%i_start:tdims%i_end,             &
                           tdims%j_start:tdims%j_end,             &
                           1:tdims%k_end ) )
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qgraup_inc(i,j,k) = 0.0
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

END IF  ! ( l_mcr_conv )

IF ( .NOT. ALLOCATED(qcf2_inc) )  ALLOCATE( qcf2_inc(1,1,1) )
IF ( .NOT. ALLOCATED(qrain_inc) )  ALLOCATE( qrain_inc(1,1,1) )
IF ( .NOT. ALLOCATED(qgraup_inc) )  ALLOCATE( qgraup_inc(1,1,1) )

! If PC2 cloud scheme is expecting increments from convection
IF ( l_calc_dxek .OR. i_convection_vn==i_cv_llcs ) THEN
                      ! Note: LLCS scheme sets these to zero even if PC2 is
                      !       off, so they need to be allocated regardless.
  ! Allocate increment arrays for cloud fraction fields.
  ALLOCATE( cf_liquid_inc ( tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                            1:tdims%k_end ) )
  ALLOCATE( cf_frozen_inc ( tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                            1:tdims%k_end ) )
  ALLOCATE( bulk_cf_inc   ( tdims%i_start:tdims%i_end,          &
                            tdims%j_start:tdims%j_end,          &
                            1:tdims%k_end ) )
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        cf_liquid_inc(i,j,k) = 0.0
        cf_frozen_inc(i,j,k) = 0.0
        bulk_cf_inc(i,j,k)   = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
ELSE
  ! Allocate to minimal size if not needed
  ALLOCATE( cf_liquid_inc (1,1,1) )
  ALLOCATE( cf_frozen_inc (1,1,1) )
  ALLOCATE( bulk_cf_inc   (1,1,1) )
END IF

! If convective momentum transport is on, initialise wind tendencies to zero
IF (l_mom) THEN
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j, k)
  DO k = pdims_s%k_start, pdims_s%k_end
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        dubydt_pout(i,j,k) = 0.0
        dvbydt_pout(i,j,k) = 0.0
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO
END IF

! Initialise output convective precip fields to zero

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, j)
DO j = 1, rows
  DO i = 1, row_length
    conv_rain(i,j) = 0.0
    conv_snow(i,j) = 0.0
  END DO
END DO
!$OMP END PARALLEL DO


!=====================================================================
! Which alternative convection scheme ? - controlled by i_convection_vn 
!=====================================================================

SELECT CASE (i_convection_vn)

!CASE (i_cv_betts_miller)
  !----------------------------------------------------------------------
  ! Scheme 10 - eg Betts Miller?  i_cv_betts_miller should be added 
  !             to cv_run_mod
  !----------------------------------------------------------------------

!  CALL betts_miller()

!CASE (i_cv_name)

  !----------------------------------------------------------------------
  ! Scheme ?? - Name    (i_cv_name) should be added to cv_run_mod
  !----------------------------------------------------------------------

  ! Either a direct call to the scheme or a call to a glue routine
  ! to interface to the alternative convection scheme.
  !  CALL another_scheme()
  !or  CALL glue_another_scheme()

CASE (i_cv_llcs)

  ! Set appropriate logical here
  IF ( i_cld_vn == i_cld_pc2 ) l_output_is_cloud=.TRUE.

  DO k = 1, pdims%k_end+1
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        phalfllcs(i,j,k) = p_layer_boundaries(i,j,k-1)
      END DO
    END DO
  END DO

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        pfullllcs(i,j,k) = p_layer_centres(i,j,k)
      END DO
    END DO
  END DO

  CALL llcs_control(theta_star, q_star, p_star, phalfllcs, pfullllcs,   &
       theta_inc, q_inc, qcl_inc, cf_liquid_inc, conv_rain)

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        theta_inc(i,j,k) = theta_inc(i,j,k) - theta_star(i,j,k)
        theta_star(i,j,k) = theta_star(i,j,k) + theta_inc(i,j,k)
        q_inc(i,j,k) = q_inc(i,j,k) - q_star(i,j,k)
        q_star(i,j,k) = q_star(i,j,k) + q_inc(i,j,k)
      END DO
    END DO
  END DO

  IF ( i_cld_vn == i_cld_pc2 ) THEN 
    ! increment cloud fields with condensed moisture
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qcl_star(i,j,k) = qcl_star(i,j,k) + qcl_inc(i,j,k)
          cf_liquid_star(i,j,k) = cf_liquid_star(i,j,k) +               &
                                  cf_liquid_inc(i,j,k)
          bulk_cf_inc(i,j,k) = cf_liquid_inc(i,j,k)
          ! assume minimum overlap with ice
          bulk_cf_star(i,j,k) = bulk_cf_star(i,j,k) +                   &
                                bulk_cf_inc(i,j,k)
        END DO
      END DO
    END DO
  END IF

  IF (cycleno == numcycles .AND. l_check_moist_inc .AND.                &
       model_type == mt_global ) THEN

    ALLOCATE(dqt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,   &
                 1:tdims%k_end))

    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          dqt(i,j,k) = ( q_inc(i,j,k) + qcl_inc(i,j,k) )* recip_timestep
        END DO
      END DO
    END DO

    scheme_name = 'LLCS convection'
    CALL check_dmoist_inc(row_length, rows, model_levels,               &
                          delta_lambda, delta_phi, timestep,            &
                          p_layer_boundaries, rho_wet_rsq,              &
                          dqt,                                          &
                          conv_rain, scheme_name,                       &
                          ep1, ep2, ep3)

    DEALLOCATE(dqt)

  END IF

CASE DEFAULT
  !----------------------------------------------------------------------
  ! Should never reach here if namelist checking is working
  !----------------------------------------------------------------------

  ! Will stop the model
  errorstatus = 10
  WRITE (cmessage,'(A,A,I6)') 'Convection scheme version value ',  &
                                  'i_convection_vn = ',i_convection_vn
  CALL Ereport ( routinename, errorstatus, cmessage)

END SELECT

!=====================================================================
!  End of section calling different schemes
!=====================================================================

!--------------------------------------------------------------------------
! Wind increments - interpolation back to uv grid in atmos_physics2
!                   so that swap-bound calls can be grouped with others
!                   to save CPU.
!--------------------------------------------------------------------------

!----------------------------------------------------------------------
! SPT scheme requires dq and dT increments from convection and any PC2
!----------------------------------------------------------------------
IF (l_retain_conv_tendencies) THEN

  ! Add SPT increments for T
  DO k=1,n_conv_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dt_conv(i,j,k) = theta_inc(i,j,k)*exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

  ! Add SPT increments for q
  DO k=1,n_conv_levels
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dq_conv(i,j,k) = q_inc(i,j,k)
      END DO
    END DO
  END DO

END IF

!----------------------------------------------------------------------
! May be expected for surface exchange 
! Note assumes the convection scheme has set ccb and dwn_flux
! if this is not done then the run should not be set up with
! isrfexcnvgust == ip_srfexwithcnv
!----------------------------------------------------------------------

IF (isrfexcnvgust == ip_srfexwithcnv) THEN
  ! Store convective mass fluxes at cloud base
  ! if required for surface exchange.
  DO j = 1, rows
    DO i = 1, row_length
      IF (ccb(i,j) > 0) THEN
        ddmfx(i,j)=dwn_flux(i,j,ccb(i,j))
      ELSE
        ddmfx(i,j)=0.0
      END IF
    END DO
  END DO
END IF

! ----------------------------------------------------------------------
! Copy convection increments into diagnostic arrays
! ----------------------------------------------------------------------
IF ( l_apply_diag ) THEN
  CALL save_conv_diags( theta_inc, q_inc, qcl_inc, qcf_inc,          &
                        bulk_cf_inc, cf_liquid_inc, cf_frozen_inc,   &
                        exner_theta_levels, l_calc_dxek )
END IF

! ----------------------------------------------------------------------
! Section to calculate PC2 scheme increments on top of the convection 
! increments as done in ni_conv_ctl
! Will only be called if the model is running with the PC2 cloud scheme on.
! ----------------------------------------------------------------------

IF (l_calc_dxek) THEN

  CALL pc2_from_conv_ctl(                                                 &
                 row_length, rows, n_conv_levels,                         &
                 ! Model switches
                 l_calc_dxek, l_q_interact, l_mr_physics,               &
                 ! Coordinate info
                 exner_theta_levels, p_layer_centres,                     &  
                 ! in/out Convection increments
                 theta_inc, q_inc, qcl_inc, qcf_inc,                      &
                 cf_liquid_inc, cf_frozen_inc, bulk_cf_inc,               &
                 ! in/out values with all increments added
                 theta_star, q_star, qcl_star, qcf_star,                  &
                 cf_liquid_star, cf_frozen_star, bulk_cf_star )

END IF  ! L_calc_dxek

! Deallocate cloud fraction increment arrays
DEALLOCATE( bulk_cf_inc )
DEALLOCATE( cf_frozen_inc )
DEALLOCATE( cf_liquid_inc )
! Deallocate prognostic hydrometeor increment arrays
DEALLOCATE( qgraup_inc )
DEALLOCATE( qrain_inc )
DEALLOCATE( qcf2_inc )


!----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!----------------------------------------------------------------------
RETURN
END SUBROUTINE other_conv_ctl

END MODULE other_conv_ctl_mod
