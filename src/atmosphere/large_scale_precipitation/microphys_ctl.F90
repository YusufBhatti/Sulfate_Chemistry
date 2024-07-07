! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Subroutine microphys_ctl
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

SUBROUTINE microphys_ctl (                                                    &

! Parallel variables
  halo_i, halo_j, off_x, off_y, global_row_length,                            &
  at_extremity,                                                               &

! model dimensions.
  row_length, rows, rhc_row_length, rhc_rows, n_rows, land_points,            &
  bl_levels,                                                                  &
  lspice_dim1,lspice_dim2,lspice_dim3,                                        &
  salt_dim1, salt_dim2, salt_dim3,                                            &
  cdnc_dim1, cdnc_dim2, cdnc_dim3,                                            &
  n_arcl_species, n_arcl_compnts, i_arcl_compnts,                             &

! Model switches
  l_sulpc_so2, l_sulpc_nh3, l_soot,                                           &
  l_biomass, l_ocff, l_nitrate,                                               &
  l_cosp_lsp,                                                                 &

! Model parameters
  rhcrit,                                                                     &

! Primary fields passed in
  t_n, theta_n, q_n, qcl_n, qcf_n, qcf2_n, qrain_n, qgraup_n,                 &
  cf_n, cfl_n, cff_n,                                                         &
  u_on_p, v_on_p, w,                                                          &
  snow_depth,                                                                 &
  land_sea_mask, ice_fract,                                                   &
  p_layer_centres, p_layer_boundaries,                                        &
  rho_r2,                                                                     &
  aerosol, flash_pot,                                                         &
  ukca_cdnc, easyaerosol_cdnc,                                                &
  dust_div1,dust_div2,dust_div3,dust_div4,dust_div5,dust_div6,                &
  so2, nh3, so4_aitken, so4_accu, so4_diss,                                   &
  aged_soot, cloud_soot, aged_bmass, cloud_bmass,                             &
  aged_ocff, cloud_ocff, nitr_acc, nitr_diss, biogenic,                       &
  sea_salt_film, sea_salt_jet, arcl,                                          &

! Other fields passed in
  ntml, cumulus,                                                              &
  fland, land_index,                                                          &

! For seeder feeder scheme
  sd_orog_land,                                                               &

! Variables for stochastic physics random parameters
  m_ci,                                                                       &

! diagnostic info

  stashwork4, stashwork21,                                                    &

! SCM diagnostics switches (dummy in full UM)
  nSCMDpkgs, L_SCMDiags,                                                      &

! Increment fields passed in/out
  t_inc, q_inc, qcl_inc, qcf_inc,                                             &
  qcf2_inc, qrain_inc, qgraup_inc,                                            &
  cf_inc, cfl_inc, cff_inc,                                                   &

! Fields required elsewhere
  ls_rain, ls_rainfrac, ls_snow, ls_graup, micro_tends,                       &
  cosp_gbx, cosp_frac_agg, n_drop_pot,                                        &
! Field for Rh crit parametrization
  ext_p_layer_centres,                                                        &
  ext_tl,                                                                     &
  ext_ql,                                                                     &
  ext_qcf,                                                                    &
  ext_ice_frac,                                                               &
  ext_land_frac,                                                              &
! BL w-variance required for mixed phase cloud calculation
  bl_w_var,                                                                   &
! error information
  error_code  )

! purpose: Interface to Atmospheric Physics parametrizations.
!          This version interfaces to cloud scheme, boundary layer,
!          hydrology, large scale rain, radiation, convection.
!          Called by atmos_physics1.

! but Not:
!          gravity wave drag, and optionally energy
!          correction. Diagnostics also missing.

!          IN/OUT etc intents to be added later.

! code description: this code is written to umdp3 programming standards.

    ! Microphysics modules

USE mphys_diags_mod,       ONLY: l_aggfr_diag, l_point_diag,                  &
                                 l_psdep_diag, l_psaut_diag,                  &
                                 l_psacw_diag, l_psacr_diag,                  &
                                 l_psaci_diag, l_psmlt_diag,                  &
                                 l_psmltevp_diag, l_praut_diag,               &
                                 l_pracw_diag, l_prevp_diag,                  &
                                 l_pgaut_diag, l_pgacw_diag,                  &
                                 l_pgacs_diag, l_pgmlt_diag,                  &
                                 l_pifrw_diag, l_piprm_diag,                  &
                                 l_pidep_diag, l_piacw_diag,                  &
                                 l_piacr_diag, l_pimlt_diag,                  &
                                 l_pimltevp_diag, l_pifall_diag,              &
                                 l_psfall_diag, l_prfall_diag,                &
                                 l_pgfall_diag, l_plset_diag,                 &
                                 l_plevpset_diag, l_pifrr_diag,               &
                                 l_refl_tot, l_refl_gr, l_refl_ice,           &
                                 l_refl_ice2, l_refl_rain, l_refl_qcl,        &
                                 l_refl_1km, l_refl_surf, l_refl_max,         &
                                 l_sfwater_diag, l_sfrain_diag,               &
                                 l_sfsnow_diag,                               &
                                 frac_agg, mphys_pts,                         &
                                 l_vm_cry_diag, l_vm_agg_diag,                &
                                 l_vm_used_diag,                              &
                                 l_vtbranch_diag,                             &
                                 psdep, psaut,                                &
                                 psacw, psacr, psaci, psmlt,                  &
                                 psmltevp, praut, pracw, prevp,               &
                                 pgaut, pgacw, pgacs, pgmlt, pifrw,           &
                                 piprm, pidep, piacw, piacr, pimlt,           &
                                 pimltevp, pifall, psfall, prfall,            &
                                 pgfall, plset, plevpset, pifrr,              &
                                 dbz_tot, dbz_g, dbz_i, dbz_i2,               &
                                 dbz_l, dbz_r, vm_cry, vm_agg,                &
                                 vtbranch_flag, vm_used,                      &
                                 sfwater, sfrain, sfsnow

USE science_fixes_mod,     ONLY: l_fix_lsp_incs_to_spt

USE mphys_bypass_mod,      ONLY: l_ref_diag

USE cloud_inputs_mod,      ONLY: l_micro_eros, i_cld_vn, i_rhcpt,             &
                                 l_inc_cld_prec
USE pc2_constants_mod,     ONLY: rhcpt_tke_based, rhcpt_horiz_var,            &
                                 i_cld_pc2
USE mphys_inputs_mod,      ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain,       &
                                 l_subgrid_qcl_mp,                            &
                                 l_subgrid_cfl_mp_by_erosion,                 &
                                 l_use_sulphate_autoconv,                     &
                                 l_orograin, l_orograin_block, nsigmasf
USE electric_inputs_mod,   ONLY: l_use_electric
USE electric_main_mod,     ONLY: electric_main

USE mphys_radar_mod,       ONLY: ref_lim

  ! General atmosphere modules
USE arcl_mod,              ONLY: npd_arcl_compnts
USE planet_constants_mod,  ONLY: cp
USE conversions_mod,       ONLY: zerodegc
USE water_constants_mod,   ONLY: lc, lf
USE dust_parameters_mod,   ONLY: ndiv, krain_dust, ksnow_dust,                &
                                 l_twobin_dust, l_dust
USE gen_phys_inputs_mod,   ONLY: l_mr_physics
USE timestep_mod,          ONLY: timestep

  ! Grid bounds module
USE atm_fields_bounds_mod, ONLY: tdims, tdims_s,                              &
                                 pdims, wdims_s
USE level_heights_mod,     ONLY: r_theta_levels, r_rho_levels

  ! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim
!$ USE omp_lib

USE Field_Types
USE UM_ParParams
USE nlstcall_mod, ONLY: ltimer

  ! COSP modules
USE cosp_types_mod,         ONLY: cosp_gridbox

USE diagnostics_lsrain_mod, ONLY: diagnostics_lsrain
USE ls_ppn_mod,             ONLY: ls_ppn
USE lsp_froude_moist_mod,   ONLY: lsp_froude_moist

USE mphys_turb_gen_mixed_phase_mod, ONLY: mphys_turb_gen_mixed_phase

! Modules to save tendencies for SPT or other schemes
USE physics_tendencies_mod,  ONLY:                                            &
    l_retain_mic_tendencies, dt_mic, dq_mic

  ! Aerosols modules
USE mass_calc_mod,      ONLY: mass_calc
USE nh3dwash_mod,       ONLY: nh3dwash
USE rainout_intctl_mod, ONLY: rainout_intctl
USE sl3dwash_mod,       ONLY: sl3dwash
USE stash_array_mod,    ONLY: si, sf
USE s_scmop_mod,        ONLY: default_streams, t_avg, t_mult,                 &
                              d_wet, d_all, scmdiag_lsp,                      &
                              scmdiag_incs
USE scmoutput_mod,      ONLY: scmoutput

USE model_domain_mod, ONLY: model_type, mt_single_column

USE nlsizes_namelist_mod, ONLY: model_levels

USE def_easyaerosol, ONLY: t_easyaerosol_cdnc

USE ls_calc_rhcrit_mod, ONLY: ls_calc_rhcrit

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

! arguments with intent in. ie: input variables.

! Parallel setup variables
INTEGER ::                                                                    &
  halo_i,                                                                     &
                 ! Size of halo in i direction.
  halo_j,                                                                     &
                 ! Size of halo in j direction.
  off_x,                                                                      &
                 ! Size of small halo in i
  off_y,                                                                      &
                 ! Size of small halo in j.
  global_row_length
                 ! number of points on a row

LOGICAL ::                                                                    &
  at_extremity(4)  ! Indicates if this processor is at north,
!                          south, east or west of the processor grid

LOGICAL ::                                                                    &
  l_sulpc_so2,                                                                &
                  ! true for sulphur cycle calcs
  l_sulpc_nh3,                                                                &
                  ! true for sulphur/nh3 calcs
  l_soot,                                                                     &
                  ! true for soot cycles calcs
  l_biomass,                                                                  &
                  ! true for biomass aerosol calcs
  l_ocff,                                                                     &
                  ! true for fossil-fuel organic carbon calcs
  l_nitrate
                  ! true for ammonium nitrate aerosol calcs

! Model dimensions
INTEGER ::                                                                    &
  row_length,                                                                 &
  rows,                                                                       &
  rhc_row_length,                                                             &
             ! = row_length if i_rhcpt = rhcpt_horiz_var, 1 otherwise.
  rhc_rows,                                                                   &
             ! = rows       if i_rhcpt = rhcpt_horiz_var, 1 otherwise.
  n_rows,                                                                     &
  land_points,                                                                &
                  ! IN No.of land points being processed, can be 0.
  bl_levels,                                                                  &
  lspice_dim1,                                                                &
                    ! Dimensions for 3D diagnostic arrays.
  lspice_dim2,                                                                &
                    ! These are set to 1 in order to save
  lspice_dim3,                                                                &
                    ! memory if the diagnostics are not used.
  salt_dim1,                                                                  &

  salt_dim2,                                                                  &
                    ! Dimensions for sea-salt aerosol arrays
  salt_dim3,                                                                  &
  n_arcl_species,                                                             &
                    ! Number of requested species within the climatology
  n_arcl_compnts,                                                             &
                    ! Corresponding number of requested components
  i_arcl_compnts(npd_arcl_compnts)
                    ! Array index of each aerosol clim component


INTEGER, INTENT(IN) :: cdnc_dim1
INTEGER, INTENT(IN) :: cdnc_dim2
INTEGER, INTENT(IN) :: cdnc_dim3
                ! Dimensions of UKCA cloud drop number concentration

! Switch for COSP LS diagnostics
LOGICAL ::  l_cosp_lsp

!  Additional arguments for SCM diagnostics which are dummy in full UM
INTEGER ::                                                                    &
  nSCMDpkgs             ! No of SCM diagnostics packages

LOGICAL ::                                                                    &
  L_SCMDiags(nSCMDpkgs) ! Logicals for SCM diagnostics packages

! Model parameters
REAL ::                                                                       &
  rhcrit(1 : tdims%k_end)  ! IN Critical relative humidity.
                                ! the values need to be tuned
                                ! for the given set of levels.
! Variable used in stochastic physics random parameters
REAL ::                                                                       &
 m_ci            ! Used to modify ice fall speed

! Primary fields passed in
! Note these are copies of prognostics from atmos_physics1
REAL, INTENT (INOUT) ::                                                       &
  t_n(      tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  q_n(      tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  qcl_n(    tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  qcf_n(    tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  qcf2_n(   tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  qrain_n(  tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  qgraup_n( tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  cf_n(     tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  cfl_n(    tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  cff_n(    tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  u_on_p(   pdims%i_start : pdims%i_end,                                      &
            pdims%j_start : pdims%j_end,                                      &
            pdims%k_start : pdims%k_end ),                                    &
  v_on_p(   pdims%i_start : pdims%i_end,                                      &
            pdims%j_start : pdims%j_end,                                      &
            pdims%k_start : pdims%k_end ),                                    &
  snow_depth( tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end ),                                  &
  ice_fract(  tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end ),                                  &
  fland( land_points )

REAL, INTENT(IN) ::    w( wdims_s%i_start : wdims_s%i_end,                    &
                          wdims_s%j_start : wdims_s%j_end,                    &
                          wdims_s%k_start : wdims_s%k_end )

REAL, INTENT(IN) :: bl_w_var( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

REAL, INTENT (IN) ::  theta_n( tdims%i_start : tdims%i_end,                   &
                               tdims%j_start : tdims%j_end,                   &
                               1 : tdims%k_end )

INTEGER ::                                                                    &
  land_index( land_points )

REAL ::                                                                       &
  sd_orog_land( land_points )

! Variables describing the low-level sub-grid orographic blocking
! estimated by the seeder feeder scheme, on landpoints
REAL ::                                                                       &
  zb_land( land_points ),                                                     &
          ! Blocked layer depth (metres)
  hmteff_land( land_points )
          ! Effective mountain height (metres) - this is the height
          ! of the sub-grid mountain rising above the blocked layer

REAL :: u_land(land_points,model_levels),                                     &
        v_land(land_points,model_levels) ,                                    &
        rho_land(land_points,model_levels) ,                                  &
        theta_land(land_points,model_levels) ,                                &
        rho_levels_land(land_points,model_levels),                            &
        theta_levels_land(land_points,model_levels)

LOGICAL ::                                                                    &
  land_sea_mask( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end )

REAL ::                                                                       &
  p_layer_boundaries( tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  0 : tdims%k_end ),                          &
            ! pressure at layer boundaries. Same as p except at
            ! bottom level = pstar, and at top = 0.

  p_layer_centres(    tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  0 : tdims%k_end ),                          &
            ! pressure at layer centres. Same as p_theta_levels
            !except bottom level = pstar, and at top = 0.

  rho_r2( tdims_s%i_start : tdims_s%i_end,                                    &
          tdims_s%j_start : tdims_s%j_end,                                    &
                        1 : tdims_s%k_end)
                             ! Density*planet_radius^2

INTEGER ::                                                                    &
  ntml( tdims%i_start : tdims%i_end,                                          &
        tdims%j_start : tdims%j_end )
                                         ! Height of diagnosed BL top

LOGICAL ::                                                                    &
  cumulus(  tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end )
                                   ! Logical indicator for convection

REAL, INTENT(INOUT) ::                                                        &
                                   ! tracer variables
  aerosol   ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
                         ! murk (total) aerosol
  flash_pot(  tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
              tdims%k_start : tdims%k_end ),                                  &
                         ! flash potential
  dust_div1(  tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
                         !dust mmr in div1
  dust_div2(  tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
                         !dust mmr in div2
  dust_div3(  tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
                         !dust mmr in div3
  dust_div4(  tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
                         !dust mmr in div4
  dust_div5(  tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
                         !dust mmr in div5
  dust_div6(  tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
                         !dust mmr in div6
  so2       ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  nh3       ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  so4_aitken( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  so4_accu  ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  so4_diss  ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  aged_soot ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  cloud_soot( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  aged_bmass( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  cloud_bmass(tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  aged_ocff(  tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  cloud_ocff( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  nitr_acc  ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end ),                              &
  nitr_diss ( tdims_s%i_start : tdims_s%i_end,                                &
              tdims_s%j_start : tdims_s%j_end,                                &
              tdims_s%k_start : tdims_s%k_end )

REAL, INTENT(IN) ::                                                           &
arcl(       tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end,                                      &
                        1 : n_arcl_compnts ),                                 &

biogenic(   tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end   ),                                  &

sea_salt_film( salt_dim1,                                                     &
               salt_dim2,                                                     &
               salt_dim3 ),                                                   &
!          Film-mode sea-salt aerosol number concentration
    sea_salt_jet(  salt_dim1,                                                 &
                   salt_dim2,                                                 &
                   salt_dim3 )
!          Jet-mode sea-salt aerosol number concentration

REAL, INTENT(IN) :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3)
!                     CDNC from UKCA

TYPE (t_easyaerosol_cdnc), INTENT(IN) :: easyaerosol_cdnc
!                     CDNC from EasyAerosol

! Diagnostics info
REAL ::  stashwork4(*)     ! STASH workspace (Microphysics)
REAL ::  stashwork21(*)    ! STASH workspace (Electrification)


! arguments with intent in/out. ie: input variables changed on output.
REAL, INTENT (INOUT) ::                                                       &
  t_inc(      tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  q_inc(      tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  qcl_inc(    tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  qcf_inc(    tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  qcf2_inc(   tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  qrain_inc(  tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  qgraup_inc( tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  cf_inc(     tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  cfl_inc(    tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end ),                                  &
  cff_inc(    tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end )

INTEGER ::                                                                    &
  error_code

! arguments with intent out. ie: output variables.
REAL, INTENT(OUT)::                                                           &
  ls_rain( tdims%i_start : tdims%i_end,                                       &
           tdims%j_start : tdims%j_end ),                                     &
  ls_snow( tdims%i_start : tdims%i_end,                                       &
           tdims%j_start : tdims%j_end ),                                     &
  ls_graup( tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end ),                                    &
  micro_tends( tdims%i_start : tdims%i_end,                                   &
               tdims%j_start : tdims%j_end, 2, bl_levels ),                   &
                      ! Tendencies from microphys within BL levels
                      ! (TL, K/s; QW, kg/kg/s)
                      ! The number 2 refers to the two categories
                      ! required by b.layer. 1 is TL and 2 is QW
 ls_rainfrac(land_points)
                      ! Rain fraction array on land points to be
                      ! passed to atmos_physics2

! Gridbox information. Input for COSP
TYPE(cosp_gridbox),INTENT(INOUT) :: cosp_gbx
REAL,INTENT(INOUT) :: cosp_frac_agg(tdims%i_start : tdims%i_end,              &
                                    tdims%j_start : tdims%j_end,              &
                                    1 : tdims%k_end )

! local variables.

CHARACTER(LEN=*), PARAMETER ::  RoutineName = 'MICROPHYS_CTL'

LOGICAL :: l_pc2_prod_qcl_mp = .FALSE. ! Apply dqcl from mp calc


! loop counters
INTEGER ::                                                                    &
  i, j, k, l,                                                                 &
  idiv, jj !loop counter for dust divisions

! open mp block size
INTEGER :: omp_block

! Diagnostic switches

! Local variables
REAL ::                                                                       &
 land_frac(tdims%i_start : tdims%i_end,                                       &
           tdims%j_start : tdims%j_end)

! Variables describing the low-level sub-grid orographic blocking
! estimated by the seeder feeder scheme - full 2D fields:
! Effective mountain height (metres)
REAL :: hmteff(tdims%i_start : tdims%i_end,                                   &
               tdims%j_start : tdims%j_end)
! Blocked layer depth (metres)
REAL :: zb(tdims%i_start : tdims%i_end,                                       &
           tdims%j_start : tdims%j_end)

REAL ::                                                                       &
  rhcpt(rhc_row_length, rhc_rows, 1 : tdims%k_end),                           &
! Critical relative humidity
    ls_rain3d(lspice_dim1,lspice_dim2,lspice_dim3),                           &
! rainfall rate out of each model level for diagonstic
    ls_snow3d(lspice_dim1,lspice_dim2,lspice_dim3),                           &
! snowfall rate out of each model level for diagonstic
   ls_graup3d(lspice_dim1, lspice_dim2, lspice_dim3),                         &
! Graupel fall rate out of each model level for diagnostic
   rainfrac3d(lspice_dim1,lspice_dim2,lspice_dim3),                           &
! Rain fraction for diagnostic
    rnout_tracer(lspice_dim1,lspice_dim2),                                    &
! Total tracer amount scavenged by rainout (kg m-2)
    lscav_tr(lspice_dim1,lspice_dim2),                                        &
! Total tracer amount scavenged by washout (kg m-2)
    rnout_soot(lspice_dim1,lspice_dim2),                                      &
! Total soot amount scavenged by rainout (kg m-2)
    lscav_soot(lspice_dim1,lspice_dim2),                                      &
! Total soot amount scavenged by washout (kg m-2)
    rnout_bmass(lspice_dim1,lspice_dim2),                                     &
! Total biomass amount scavenged by rainout (kg m-2)
    lscav_bmass(lspice_dim1,lspice_dim2),                                     &
! Total biomass amount scavenged by washout (kg m-2)
    rnout_ocff(lspice_dim1, lspice_dim2),                                     &
! Total fossil-fuel organic carbon amount scavenged by rainout (kg m-2)
    lscav_ocff(lspice_dim1, lspice_dim2),                                     &
! Total fossil-fuel organic carbon amount scavenged by washout (kg m-2)
    rnout_nitrate(lspice_dim1,lspice_dim2),                                   &
! Total ammonium nitrate amount scavenged by rainout (kg[N] m-2)
    lscav_nitrate(lspice_dim1,lspice_dim2)
! Total ammonium nitrate amount scavenged by washout (kg[N] m-2)


REAL ::                                                                       &
  t_work(   tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  q_work(   tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  qcl_work( tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  qcf_work( tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  cf_work(  tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  cfl_work( tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
  cff_work( tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
n_drop_3d(  tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
                 ! 3D calculated droplet number
rhodz_dry(  tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end ),                                    &
                 ! Dry density
rhodz_moist(tdims%i_start : tdims%i_end,                                      &
            tdims%j_start : tdims%j_end,                                      &
                        1 : tdims%k_end )
                 ! Moist density

REAL :: dqcl_mp( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end,                                 &
                             1 : tdims%k_end )



REAL, INTENT(OUT) ::                                                          &
  n_drop_pot( tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end,                                    &
                          1 : tdims%k_end )
!   Potential droplet number (includes values where cloud not present)

REAL, INTENT(IN) ::                                                           &
  ext_p_layer_centres(0:rhc_row_length+1,0:rhc_rows+1,                        &
                                         0:tdims%k_end),                      &
  ext_tl(0:rhc_row_length+1, 0:rhc_rows+1, 1: tdims%k_end),                   &
  ext_ql(0:rhc_row_length+1, 0:rhc_rows+1, 1: tdims%k_end),                   &
  ext_qcf(0:rhc_row_length+1,0:rhc_rows+1, 1: tdims%k_end),                   &
  ext_ice_frac(0:rhc_row_length+1,0:rhc_rows+1),                              &
  ext_land_frac(0:rhc_row_length+1,0:rhc_rows+1)

REAL ::                                                                       &
 dm( tdims%i_start : tdims%i_end,                                             &
     tdims%j_start : tdims%j_end,                                             &
                 1 : tdims%k_end  )   ! mass air p.u. area in lev

    ! Local work arrays for additional microphysics fields if in use

REAL, ALLOCATABLE :: qcf2_work(:,:,:)
REAL, ALLOCATABLE :: qrain_work(:,:,:)
REAL, ALLOCATABLE :: qgraup_work(:,:,:)

REAL ::                                                                       &
  tinc_np1(   tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end ),                                  &
  qinc_np1(   tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end ),                                  &
  qclinc_np1( tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end ),                                  &
  qcfinc_np1( tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end )

! Local work array for increment diagnostics
REAL ::                                                                       &
  work_3d( tdims%i_start : tdims%i_end,                                       &
           tdims%j_start : tdims%j_end,                                       &
                       1 : tdims%k_end )

! Scavenged tracers (in column) for diagnostics
REAL ::                                                                       &
  lscav_so2( tdims%i_start : tdims%i_end,                                     &
             tdims%j_start : tdims%j_end ),                                   &
  lscav_nh3( tdims%i_start : tdims%i_end,                                     &
             tdims%j_start : tdims%j_end ),                                   &
  lscav_so4ait( tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end ),                                &
  lscav_so4acc( tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end ),                                &
  lscav_so4dis( tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end ),                                &
  lscav_dust_all( tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end, ndiv ),                        &
  dust_all(       tdims%i_start : tdims%i_end,                                &
                  tdims%j_start : tdims%j_end,                                &
                              1 : tdims%k_end, ndiv )

REAL ::                                                                       &
rainrate,                                                                     &
               !rate corrected for possible -ive precip
snowrate,                                                                     &
               !rate corrected for possible -ive precip
rate,                                                                         &

delta_dust

REAL :: deltaz (tdims%i_start : tdims%i_end,                                  &
                tdims%j_start : tdims%j_end,                                  &
                            1 : tdims%k_end )

! Turbulent mixed phase diagnostics
REAL, ALLOCATABLE :: qcl_mpt(:,:,:) ! qcl increment by mixed-phase scheme
REAL, ALLOCATABLE :: tau_d(:,:,:)   ! Turbulent decorrelation timescale
REAL, ALLOCATABLE :: inv_prt(:,:,:) ! Inverse phase-relaxation timescale
REAL, ALLOCATABLE :: disprate(:,:,:)! Diagnosed turbulent dissipation rate
REAL, ALLOCATABLE :: inv_mt(:,:,:)  ! Inverse mixing timescale
REAL, ALLOCATABLE :: si_avg(:,:,:)  ! Mean of subgrid supersaturation PDF
REAL, ALLOCATABLE :: dcfl_mp(:,:,:) ! cfl diagnosed by mixed phase scheme
REAL, ALLOCATABLE :: sigma2_s(:,:,:)! Variance of subgrid supersaturation PDF
REAL, ALLOCATABLE :: cfl_copy(:,:,:) ! Copy of cfl_work
REAL, ALLOCATABLE :: cfl_incr_mp(:,:,:) ! cfl increment from mixed phase
!                                           scheme (PC2)
REAL, ALLOCATABLE :: qcl_copy(:,:,:) ! Copy of qcl_work
REAL, ALLOCATABLE :: qcl_incr_mp_pc2(:,:,:) ! qcl increment from mixed phase
!                                               scheme (PC2)

! Work variables for tracer variables

REAL ::                                                                       &
  so4_aitken_work ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                               1  : tdims%k_end ),                            &
  so4_accu_work   ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  so4_diss_work   ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  aged_bmass_work ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  cloud_bmass_work( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  aged_ocff_work  ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  cloud_ocff_work ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  nitr_acc_work   ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  nitr_diss_work  ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ),                            &
  aerosol_work    ( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end )

! Microphysical process rate diagnostics
! Note: These arrays will only increase memory usage and are
!       only referenced if the particular diagnostic is active
!       (i.e. if the associated logical is .true., set from SF
!       STASH flag.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

omp_block = (tdims%j_end - tdims%j_start) + 1

! ----------------------------------------------------------------------
! Section Microphysics. Call microphys_ctl routine
! ----------------------------------------------------------------------
! Call timer for cloud code
IF (ltimer) CALL timer ('AP1M LS Cloud',5)

! Store values in work arrays

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, i_cld_vn,t_work, q_work, t_n, q_n, cf_n, qcl_work,&
!$OMP         qcf_work, cf_work, cfl_work, cff_work, qcl_n, qcf_n,     &
!$OMP         cfl_n, cff_n, i_rhcpt, cp )                              &
!$OMP PRIVATE( i, j, k )
DO k = 1, tdims%k_end

  DO j = tdims%j_start, tdims%j_end

    DO i = tdims%i_start, tdims%i_end

      IF (i_cld_vn == i_cld_pc2 .OR. i_rhcpt == rhcpt_tke_based .OR.   &
          .NOT. l_inc_cld_prec) THEN

        ! No diagnostic cloud scheme is required.
        !Store temperature and
        ! vapour content in T_work and q_work

        t_work(i,j,k)   = t_n(i,j,k)
        q_work(i,j,k)   = q_n(i,j,k)

      ELSE  ! i_cld_pc2
        ! inputs to diagnostic cloud scheme are
        ! Tl_star, qT_star
        ! outputs from diag. cloud scheme
        ! are cloud_fraction, T, q, qcl, qcf
        ! output T is held in theta_star
        ! Calculate Tl, store in T_work
        ! Calculate qT, store in q_work

        t_work(i,j,k) = t_n(i,j,k) -                                          &
                        (lc * qcl_n(i,j,k) ) / cp
        q_work(i,j,k) = q_n(i,j,k) + qcl_n(i,j,k)

      END IF  ! i_cld_pc2

      qcl_work(i,j,k) = qcl_n(i,j,k)
      qcf_work(i,j,k) = qcf_n(i,j,k)
      cf_work(i,j,k)  = cf_n(i,j,k)
      cfl_work(i,j,k) = cfl_n(i,j,k)
      cff_work(i,j,k) = cff_n(i,j,k)

    END DO
  END DO
END DO
!$OMP END PARALLEL DO

IF (l_micro_eros) THEN
  ! DEPENDS ON: pc2_turbulence_ctl
  CALL pc2_turbulence_ctl (                                                   &
  ! Primary fields passed in, updated on exit
       T_work, q_work, qcl_work, cf_work, cfl_work, cff_work,                 &
       p_layer_centres(1,1,1),                                                &
       dqcl_mp,                                                               &
  ! diagnostic info
       stashwork4,                                                            &
  ! SCM diagnostics switches (dummy in full UM)
       nSCMDpkgs, L_SCMDiags,                                                 &
  ! Increment fields passed in, unchanged on exit
       t_inc, q_inc, qcl_inc, cf_inc, cfl_inc, l_pc2_prod_qcl_mp )

END IF  ! l_micro_eros

IF (l_mcr_qcf2) THEN

      ! Second cloud ice variable in use

  ALLOCATE ( qcf2_work( tdims%i_start : tdims%i_end,                          &
                        tdims%j_start : tdims%j_end,                          &
                                    1 : tdims%k_end ) )

      ! Add qcf and qcf2 to get total ice for cloud scheme
      ! qcf_work and qcf2_work are reset after call to ls_cld

  qcf_work(:,:,:) = qcf_n(:,:,:) + qcf2_n(:,:,:)

ELSE

  ALLOCATE ( qcf2_work(1,1,1) )

END IF

IF (l_mcr_qrain) THEN

      ! Prognostic rain in use

  ALLOCATE ( qrain_work( tdims%i_start : tdims%i_end,                         &
                         tdims%j_start : tdims%j_end,                         &
                                     1 : tdims%k_end ) )

!$OMP PARALLEL DO DEFAULT(NONE)  SCHEDULE(STATIC)                     &
!$OMP SHARED( tdims, qrain_work, qrain_n )                            &
!$OMP PRIVATE ( i, j, k )
  DO k = 1,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        qrain_work(i,j,k) = qrain_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE

  ALLOCATE ( qrain_work(1,1,1) )

END IF

IF (l_mcr_qgraup) THEN

      ! Prognostic graupel in use

  ALLOCATE ( qgraup_work(tdims%i_start : tdims%i_end,                         &
                         tdims%j_start : tdims%j_end,                         &
                                     1 : tdims%k_end ) )

!$OMP PARALLEL DO DEFAULT(NONE)  SCHEDULE(STATIC)                     &
!$OMP SHARED( tdims, qgraup_work, qgraup_n )                          &
!$OMP PRIVATE ( i, j, k )
  DO k = 1,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        qgraup_work(i,j,k) = qgraup_n(i,j,k)
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

ELSE

  ALLOCATE ( qgraup_work(1,1,1) )

END IF

IF (model_type /= mt_single_column) THEN
  l_aggfr_diag    = sf(100,4)
  l_point_diag    = sf(101,4)
  l_vm_cry_diag   = sf(102,4)
  l_vm_agg_diag   = sf(103,4)
  l_vtbranch_diag = sf(104,4)
  l_vm_used_diag  = sf(105,4)
  l_pifrw_diag    = sf(240,4)
  l_piprm_diag    = sf(241,4)
  l_pidep_diag    = sf(243,4)
  l_psdep_diag    = sf(245,4)
  l_piacw_diag    = sf(247,4)
  l_psacw_diag    = sf(248,4)
  l_piacr_diag    = sf(249,4)
  l_psacr_diag    = sf(250,4)
  l_pimltevp_diag = sf(251,4)
  l_psmltevp_diag = sf(252,4)
  l_pimlt_diag    = sf(253,4)
  l_psmlt_diag    = sf(254,4)
  l_psaut_diag    = sf(255,4)
  l_psaci_diag    = sf(256,4)
  l_praut_diag    = sf(257,4)
  l_pracw_diag    = sf(258,4)
  l_prevp_diag    = sf(259,4)
  l_pgaut_diag    = sf(260,4)
  l_pgacw_diag    = sf(261,4)
  l_pgacs_diag    = sf(262,4)
  l_pgmlt_diag    = sf(263,4)
  l_pifall_diag   = sf(265,4)
  l_psfall_diag   = sf(266,4)
  l_prfall_diag   = sf(267,4)
  l_pgfall_diag   = sf(268,4)
  l_plset_diag    = sf(269,4)
  l_plevpset_diag = sf(270,4)
  l_pifrr_diag    = sf(271,4)

  l_refl_surf     = sf(110,4)
  l_refl_max      = sf(111,4)
  l_refl_1km      = sf(112,4)
  l_refl_gr       = sf(113,4)
  l_refl_ice      = sf(114,4)
  l_refl_ice2     = sf(115,4)
  l_refl_rain     = sf(116,4)
  l_refl_qcl      = sf(117,4)
  l_refl_tot      = sf(118,4)

  l_sfwater_diag  = sf(400,4)
  l_sfrain_diag   = sf(401,4)
  l_sfsnow_diag   = sf(402,4)
END IF ! model_type

IF (l_refl_surf .OR. l_refl_tot  .OR. l_refl_1km  .OR. l_refl_gr .OR.         &
    l_refl_ice  .OR. l_refl_ice2 .OR. l_refl_rain .OR. l_refl_qcl .OR.        &
    l_refl_max )                                                              &
    l_ref_diag = .TRUE.

!==============================================================
! ALLOCATE arrays for required microphysics diagnostics
!==============================================================

      ! Aggregate Fraction
IF (l_aggfr_diag .OR. l_cosp_lsp) THEN
  ALLOCATE ( frac_agg( tdims%i_start : tdims%i_end,                           &
                       tdims%j_start : tdims%j_end,                           &
                                   1 : tdims%k_end ) )
  frac_agg(:,:,:) = 0.0
ELSE
  ALLOCATE ( frac_agg(1,1,1) )
END IF

    ! Microphysics Points Used
IF (l_point_diag) THEN
  ALLOCATE ( mphys_pts( tdims%i_start : tdims%i_end,                          &
                        tdims%j_start : tdims%j_end,                          &
                                    1 : tdims%k_end ) )
  mphys_pts(:,:,:) = .FALSE.
ELSE
  ALLOCATE ( mphys_pts(1,1,1) )
END IF

    ! Homogeneous freezing nucl.
IF (l_pifrw_diag) THEN
  ALLOCATE ( pifrw( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pifrw(:,:,:) = 0.0
ELSE
  ALLOCATE ( pifrw(1,1,1) )
END IF

    ! Homoegenous freezing of rain
IF (l_pifrr_diag) THEN
  ALLOCATE ( pifrr( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pifrr(:,:,:) = 0.0
ELSE
  ALLOCATE ( pifrr(1,1,1) )
END IF

    ! Heterogeneous nucl.
IF (l_piprm_diag) THEN
  ALLOCATE ( piprm( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  piprm(:,:,:) = 0.0
ELSE
  ALLOCATE ( piprm(1,1,1) )
END IF

    ! Deposition of vapour to ice
IF (l_pidep_diag) THEN
  ALLOCATE ( pidep( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pidep(:,:,:) = 0.0
ELSE
  ALLOCATE ( pidep(1,1,1) )
END IF

    ! Deposition of vapour to snow
IF (l_psdep_diag) THEN
  ALLOCATE ( psdep( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  psdep(:,:,:) = 0.0
ELSE
  ALLOCATE ( psdep(1,1,1) )
END IF

    ! Accretion of liq. by ice
IF (l_piacw_diag) THEN
  ALLOCATE ( piacw( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  piacw(:,:,:) = 0.0
ELSE
  ALLOCATE ( piacw(1,1,1) )
END IF

    ! Accretion of liq. by snow
IF (l_psacw_diag) THEN
  ALLOCATE ( psacw( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  psacw(:,:,:) = 0.0
ELSE
  ALLOCATE ( psacw(1,1,1) )
END IF

    ! Collection of rain by ice
IF (l_piacr_diag) THEN
  ALLOCATE ( piacr( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  piacr(:,:,:) = 0.0
ELSE
  ALLOCATE ( piacr(1,1,1) )
END IF

    ! Collection of rain by snow
IF (l_psacr_diag) THEN
  ALLOCATE ( psacr( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  psacr(:,:,:) = 0.0
ELSE
  ALLOCATE ( psacr(1,1,1) )
END IF

    ! Evaporation of melting ice
IF (l_pimltevp_diag) THEN
  ALLOCATE ( pimltevp( tdims%i_start : tdims%i_end,                           &
                       tdims%j_start : tdims%j_end,                           &
                                   1 : tdims%k_end ) )
  pimltevp(:,:,:) = 0.0
ELSE
  ALLOCATE ( pimltevp(1,1,1) )
END IF

    ! Evap. of melting aggregates
IF (l_psmltevp_diag) THEN
  ALLOCATE ( psmltevp( tdims%i_start : tdims%i_end,                           &
                       tdims%j_start : tdims%j_end,                           &
                                   1 : tdims%k_end ) )
  psmltevp(:,:,:) = 0.0
ELSE
  ALLOCATE ( psmltevp(1,1,1) )
END IF

    ! Melting of ice crystals
IF (l_pimlt_diag) THEN
  ALLOCATE ( pimlt( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pimlt(:,:,:) = 0.0
ELSE
  ALLOCATE ( pimlt(1,1,1) )
END IF

    ! Melting of snow aggregates
IF (l_psmlt_diag) THEN
  ALLOCATE ( psmlt(tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ) )
  psmlt(:,:,:) = 0.0
ELSE
  ALLOCATE ( psmlt(1,1,1) )
END IF

    ! Autoconversion of snow
IF (l_psaut_diag) THEN
  ALLOCATE ( psaut( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  psaut(:,:,:) = 0.0
ELSE
  ALLOCATE ( psaut(1,1,1) )
END IF

    ! Collection of ice crystals
IF (l_psaci_diag) THEN
  ALLOCATE ( psaci( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  psaci(:,:,:) = 0.0
ELSE
  ALLOCATE ( psaci(1,1,1) )
END IF

    ! Autoconversion of cloud
IF (l_praut_diag) THEN
  ALLOCATE ( praut( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  praut(:,:,:) = 0.0
ELSE
  ALLOCATE ( praut(1,1,1) )
END IF

    ! Accretion of liq. by rain
IF (l_pracw_diag) THEN
  ALLOCATE ( pracw( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pracw(:,:,:) = 0.0
ELSE
  ALLOCATE ( pracw(1,1,1) )
END IF

    ! Evaporation of rain
IF (l_prevp_diag) THEN
  ALLOCATE ( prevp( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  prevp(:,:,:) = 0.0
ELSE
  ALLOCATE ( prevp(1,1,1) )
END IF

    ! Autoconversion of graupel
IF (l_pgaut_diag) THEN
  ALLOCATE ( pgaut( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pgaut(:,:,:) = 0.0
ELSE
  ALLOCATE ( pgaut(1,1,1) )
END IF

    ! Accretion of liq. by graup
IF (l_pgacw_diag) THEN
  ALLOCATE ( pgacw( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pgacw(:,:,:) = 0.0
ELSE
  ALLOCATE ( pgacw(1,1,1) )
END IF

    ! Collection of snow by graup
IF (l_pgacs_diag) THEN
  ALLOCATE ( pgacs( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pgacs(:,:,:) = 0.0
ELSE
  ALLOCATE ( pgacs(1,1,1) )
END IF

    ! Melting of graupel
IF (l_pgmlt_diag) THEN
  ALLOCATE ( pgmlt( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  pgmlt(:,:,:) = 0.0
ELSE
  ALLOCATE ( pgmlt(1,1,1) )
END IF

    ! Sedimentation of ice crystals
IF (l_pifall_diag) THEN
  ALLOCATE ( pifall( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  pifall(:,:,:) = 0.0
ELSE
  ALLOCATE ( pifall(1,1,1) )
END IF

    ! Sedimentation of aggregates
IF (l_psfall_diag) THEN
  ALLOCATE ( psfall( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  psfall(:,:,:) = 0.0
ELSE
  ALLOCATE ( psfall(1,1,1) )
END IF

    ! Sedimentation of rain
IF (l_prfall_diag) THEN
  ALLOCATE ( prfall( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  prfall(:,:,:) = 0.0
ELSE
  ALLOCATE ( prfall(1,1,1) )
END IF

    ! Sedimentation of graupel
IF (l_pgfall_diag) THEN
  ALLOCATE ( pgfall( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  pgfall(:,:,:) = 0.0
ELSE
  ALLOCATE ( pgfall(1,1,1) )
END IF

    ! Droplet settling of liquid water
IF (l_plset_diag) THEN
  ALLOCATE ( plset( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  plset(:,:,:) = 0.0
ELSE
  ALLOCATE ( plset(1,1,1) )
END IF

    ! Evaporated settled droplets
IF (l_plevpset_diag) THEN
  ALLOCATE ( plevpset( tdims%i_start : tdims%i_end,                           &
                       tdims%j_start : tdims%j_end,                           &
                                   1 : tdims%k_end ) )
  plevpset(:,:,:) = 0.0
ELSE
  ALLOCATE ( plevpset(1,1,1) )
END IF

    ! Fallspeed branch to use. 0=crystals; 1=aggregates
IF (l_vtbranch_diag) THEN
  ALLOCATE ( vtbranch_flag( tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end,                      &
                                        1 : tdims%k_end ) )
  vtbranch_flag(:,:,:) = -1.0
ELSE
  ALLOCATE ( vtbranch_flag(1,1,1) )
END IF

    ! Mass-weighted fallspeed with crystal parameters
IF (l_vm_cry_diag) THEN
  ALLOCATE ( vm_cry( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  vm_cry(:,:,:) = -1.0
ELSE
  ALLOCATE ( vm_cry(1,1,1) )
END IF

    ! Mass-weighted fallspeed with aggregate parameters
IF (l_vm_agg_diag) THEN
  ALLOCATE ( vm_agg( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  vm_agg(:,:,:) = -1.0
ELSE
  ALLOCATE ( vm_agg(1,1,1) )
END IF

    ! Fallspeed used
IF (l_vm_used_diag) THEN
  ALLOCATE ( vm_used( tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  vm_used(:,:,:) = 0.0
ELSE
  ALLOCATE ( vm_used(1,1,1) )
END IF


! Allocate turbulent mixed phase diagnostic space.
IF (l_subgrid_qcl_mp) THEN
  ALLOCATE ( qcl_mpt( tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  qcl_mpt(:,:,:) = 0.0

  ALLOCATE ( tau_d(   tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  tau_d(:,:,:) = 0.0

  ALLOCATE ( inv_prt( tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  inv_prt(:,:,:) = 0.0

  ALLOCATE ( disprate( tdims%i_start : tdims%i_end,                           &
                       tdims%j_start : tdims%j_end,                           &
                                   1 : tdims%k_end ) )
  disprate(:,:,:) = 0.0

  ALLOCATE ( inv_mt(  tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  inv_mt(:,:,:) = 0.0

  ALLOCATE ( si_avg(  tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  si_avg(:,:,:) = 0.0

  ALLOCATE ( dcfl_mp( tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  dcfl_mp(:,:,:) = 0.0

  ALLOCATE ( sigma2_s( tdims%i_start : tdims%i_end,                           &
                       tdims%j_start : tdims%j_end,                           &
                                   1 : tdims%k_end ) )
  sigma2_s(:,:,:) = 0.0

ELSE

  ! Failsafe option for when option is switched off

  ALLOCATE ( qcl_mpt( 1,1,1) )
  ALLOCATE ( tau_d(   1,1,1) )
  ALLOCATE ( inv_prt( 1,1,1) )
  ALLOCATE ( disprate(1,1,1) )
  ALLOCATE ( inv_mt(  1,1,1) )
  ALLOCATE ( si_avg(  1,1,1) )
  ALLOCATE ( dcfl_mp( 1,1,1) )
  ALLOCATE ( sigma2_s(1,1,1) )

END IF

IF (model_type /= mt_single_column) THEN
  IF ( l_subgrid_qcl_mp .AND. l_subgrid_cfl_mp_by_erosion  .AND.              &
      sf(302,4) ) THEN
    ! Need to allocate space for the diagnostic
    ALLOCATE ( cfl_incr_mp( tdims%i_start : tdims%i_end,                      &
                            tdims%j_start : tdims%j_end,                      &
                                        1 : tdims%k_end ) )

    ALLOCATE ( qcl_incr_mp_pc2( tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                  &
                                            1 : tdims%k_end ) )

    ! Initialise them to zero
    cfl_incr_mp(:,:,:)     = 0.0
    qcl_incr_mp_pc2(:,:,:) = 0.0
  ELSE
    ALLOCATE ( cfl_incr_mp(1,1,1)     )
    ALLOCATE ( qcl_incr_mp_pc2(1,1,1) )
  END IF
ELSE
  ALLOCATE ( cfl_incr_mp(1,1,1)     )
  ALLOCATE ( qcl_incr_mp_pc2(1,1,1) )
END IF

! Allocate reflectivity diagnostics where required
!-----------------------------------------------------------
! N.B. There is no need to allocate surface reflectivity or
! reflectivity at 1 km AGL as they are only calculated later

! Allocate all diagnostics if any is required. This avoids issues
! with arrays being required for total diagnostics but then not being
! allocated properly.
IF (l_ref_diag) THEN
  ALLOCATE (dbz_tot( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  dbz_tot(:,:,:) = ref_lim

  ALLOCATE (dbz_g( tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ) )
  dbz_g(:,:,:) = ref_lim

  ALLOCATE (dbz_i( tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ) )
  dbz_i(:,:,:) = ref_lim

  ALLOCATE (dbz_i2( tdims%i_start : tdims%i_end,                              &
                    tdims%j_start : tdims%j_end,                              &
                                1 : tdims%k_end ) )
  dbz_i2(:,:,:) = ref_lim

  ALLOCATE (dbz_r( tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ) )
  dbz_r(:,:,:) = ref_lim

  ALLOCATE (dbz_l( tdims%i_start : tdims%i_end,                               &
                   tdims%j_start : tdims%j_end,                               &
                               1 : tdims%k_end ) )
  dbz_l(:,:,:) = ref_lim

ELSE
  ALLOCATE ( dbz_tot(1,1,1) )
  ALLOCATE ( dbz_g(  1,1,1) )
  ALLOCATE ( dbz_i(  1,1,1) )
  ALLOCATE ( dbz_i2( 1,1,1) )
  ALLOCATE ( dbz_r(  1,1,1) )
  ALLOCATE ( dbz_l(  1,1,1) )
END IF

! Allocate required seeder feeder diagnostics
! Subgrid orographic water mixing ratio
IF (l_sfwater_diag) THEN
  ALLOCATE ( sfwater( tdims%i_start : tdims%i_end,                            &
                      tdims%j_start : tdims%j_end,                            &
                                  1 : tdims%k_end ) )
  sfwater(:,:,:) = 0.0
ELSE
  ALLOCATE ( sfwater(1,1,1) )
END IF

! Subgrid orographic rain production
IF (l_sfrain_diag) THEN
  ALLOCATE ( sfrain( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  sfrain(:,:,:) = 0.0
ELSE
  ALLOCATE ( sfrain(1,1,1) )
END IF

! Subgrid orographic snow production
IF (l_sfsnow_diag) THEN
  ALLOCATE ( sfsnow( tdims%i_start : tdims%i_end,                             &
                     tdims%j_start : tdims%j_end,                             &
                                 1 : tdims%k_end ) )
  sfsnow(:,:,:) = 0.0
ELSE
  ALLOCATE ( sfsnow(1,1,1) )
END IF


! ----------------------------------------------------------------------
! Section CLD.1.a Calculate diagnostic RHcrit or read as namelist param.
! ----------------------------------------------------------------------

! i_rhcpt_if1:
IF (i_rhcpt == rhcpt_horiz_var) THEN
  !       RHCRIT is 3D diagnosed variable

  CALL ls_calc_rhcrit( ext_p_layer_centres,                                   &
  !            Array dimensions
         rhc_row_length, rhc_rows,                                            &
  !            Prognostic Fields
         ext_tl, ext_ql, ext_qcf, ext_land_frac, ext_ice_frac,                &
  !            Output
         rhcpt)

ELSE ! i_rhcpt = rhcpt_off or rhcpt_tke_based
  !       RHCRIT is 1D Parameter read in from namelist
  DO k = 1, tdims%k_end
    rhcpt(1,1,k) = rhcrit(k)
  END DO
END IF  ! i_rhcpt_if1

IF (i_cld_vn /= i_cld_pc2 .AND. i_rhcpt /= rhcpt_tke_based .AND.              &
    l_inc_cld_prec) THEN

  ! ----------------------------------------------------------------------
  ! Section CLD.1.b Calculate (2A) large-scale cloud fraction/ condensate.
  ! This call is redundant as it should calculate the exact same cloud
  ! fraction as at the end of the previous time step.
  ! However, since the call is to ls_cld and not to ls_arcld, no area cloud
  ! fraction is taken into account here and the result will be different and
  ! inconsistent if i_cld_area != acf_off. The switch l_inc_cld_prec should 
  ! be set to .FALSE. if this call is not needed, but for now its default 
  ! value is set to .TRUE. to not affect UKV-type configurations. Longer term, 
  ! this call should be removed. 
  ! ----------------------------------------------------------------------
  ! DEPENDS ON: ls_cld
  CALL ls_cld( p_layer_centres(1,1,1), rhcpt,                                 &
               tdims%k_end, bl_levels,                                        &
               rhc_row_length, rhc_rows,                                      &
               ntml, cumulus, l_mr_physics, t_work(1,1,1),                    &
               cf_work(1,1,1), q_work(1,1,1), qcf_work(1,1,1),                &
               qcl_work(1,1,1), cfl_work(1,1,1), cff_work(1,1,1),             &
               error_code )

END IF  ! .not. i_cld_pc2

! Call timer for cloud code
IF (ltimer) CALL timer ('AP1M LS Cloud',6)

IF (error_code  ==  0 ) THEN

  ! ----------------------------------------------------------------------
  ! Section LSP.1 Call Large Scale precipitation scheme.
  ! ----------------------------------------------------------------------
  ! Call timer for LS rain code
  IF (ltimer) CALL timer ('AP1M LS Rain ',5)

  ! Set global land fraction

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)
!$OMP DO SCHEDULE(STATIC)
  DO j = tdims%j_start, tdims%j_end

    DO i = tdims%i_start, tdims%i_end

      land_frac(i,j) = 0.0

    END DO

  END DO
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, land_points

    j = (land_index(k)-1)/ ( tdims%i_len ) + 1

    i = land_index(k) - (j-1)*( tdims%i_len )

    land_frac(i,j) = fland(k)

  END DO
!$OMP END DO

!$OMP END PARALLEL


! Estimate sub-grid orographic blocking for Seeder Feeder

! Initialise 2D variables passed down to ls_ppn
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      hmteff(i,j)  = 0.0
      zb(i,j)      = 0.0
    END DO
  END DO


  IF (l_orograin) THEN

!   Initially set zb and hmteff as if no blocking
    DO l = 1, land_points
      zb_land(l) = 0.0
      hmteff_land(l) = nsigmasf * sd_orog_land(l)
    END DO


    IF (l_orograin_block) THEN

!     Gather variables required by lsp_froude_moist onto landpoints
      DO k=1,model_levels
        DO l=1,land_points

          j = (land_index(l)-1)/ ( tdims%i_len ) + 1
          i = land_index(l) - (j-1)*( tdims%i_len )

          u_land(l,k) = u_on_p(i,j,k)
          v_land(l,k) = v_on_p(i,j,k)
          rho_land(l,k) = rho_r2(i,j,k)
          theta_land(l,k) = theta_n(i,j,k)

          rho_levels_land(l,k) = r_rho_levels(i,j,k) -                        &
                               r_theta_levels(i,j,0)

          theta_levels_land(l,k) = r_theta_levels(i,j,k) -                    &
                                   r_theta_levels(i,j,0)

        END DO
      END DO

!     Estimate amount of low-level blocking and effective mountain height
      CALL lsp_froude_moist(  model_levels, land_points,                      &
                          u_land, v_land, rho_land, theta_land,               &
                          rho_levels_land, theta_levels_land,                 &
                          hmteff_land, zb_land )

    END IF  !l_orograin_block


!   Put zb and hmteff into full 2D arrays
    DO l = 1, land_points

      j = (land_index(l)-1)/ ( tdims%i_len ) + 1
      i = land_index(l) - (j-1)*( tdims%i_len )

      zb(i,j)     = zb_land(l)
      hmteff(i,j) = hmteff_land(l)

    END DO

  END IF  !l_orograin
! End of seeder-feeder sub-grid orographic blocking calculations


  IF (l_mcr_qcf2) THEN  ! Second cloud ice variable in use
      ! Store ice variables in qcf_work and qcf2_work for LS_PPN
    qcf_work(:,:,:)  = qcf_n(:,:,:)
    qcf2_work(:,:,:) = qcf2_n(:,:,:)
  END IF

  ! Copy the data of the input fields using tdims_s bounds to fields
  ! using tdims bounds as needed by ls_ppn.

!$OMP PARALLEL DEFAULT(NONE)                                          &
!$OMP SHARED( tdims,                                                  &
!$OMP         so4_accu, so4_accu_work, so4_diss, so4_diss_work,       &
!$OMP         aged_bmass_work, aged_bmass, cloud_bmass_work,          &
!$OMP         cloud_bmass, aged_ocff_work, aged_ocff,                 &
!$OMP         cloud_ocff_work, cloud_ocff, nitr_acc_work,             &
!$OMP         nitr_acc, nitr_diss_work, nitr_diss,aerosol_work,       &
!$OMP         aerosol, l_use_sulphate_autoconv)                       &
!$OMP PRIVATE ( i, j, k )

  IF( l_use_sulphate_autoconv ) THEN
!$OMP DO SCHEDULE(STATIC)
    DO k = 1,tdims%k_end
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          ! so4_aitken_work(i,j,k)  = so4_aitken(i,j,k) ! Not used in ls_ppn
          so4_accu_work(i,j,k)    = so4_accu(i,j,k)
          so4_diss_work(i,j,k)    = so4_diss(i,j,k)
          aged_bmass_work(i,j,k)  = aged_bmass(i,j,k)
          cloud_bmass_work(i,j,k) = cloud_bmass(i,j,k)
          aged_ocff_work(i,j,k)   = aged_ocff(i,j,k)
          cloud_ocff_work(i,j,k)  = cloud_ocff(i,j,k)
          nitr_acc_work(i,j,k)    = nitr_acc(i,j,k)
          nitr_diss_work(i,j,k)   = nitr_diss(i,j,k)
        END DO
      END DO
    END DO
!$OMP END DO NOWAIT
  END IF

!$OMP DO SCHEDULE(STATIC)
  DO k = 1,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        aerosol_work(i,j,k)     = aerosol(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT

!$OMP END PARALLEL


  ! theta_star holds current estimate of temperature at new time level.
  ! 3A
  CALL ls_ppn(                                                                &
              p_layer_boundaries, p_layer_centres(1,1,1),                     &
              land_sea_mask, deltaz,                                          &
  ! primary fields and
  ! cloud fractions
               cf_work, cfl_work, cff_work,                                   &
               rhcpt,                                                         &
               lspice_dim1,lspice_dim2,lspice_dim3,                           &
               rho_r2, q_work, qcf_work, qcl_work, t_work,                    &
               qcf2_work, qrain_work, qgraup_work,                            &
  ! Wind field and grid size
  ! for lateral displacement
  ! of falling ice by shear
               u_on_p, v_on_p,                                                &
  ! Aerosol variables
               sea_salt_film, sea_salt_jet,                                   &
               salt_dim1, salt_dim2, salt_dim3,                               &
               ukca_cdnc,                                                     &
               cdnc_dim1, cdnc_dim2, cdnc_dim3,                               &
               easyaerosol_cdnc,                                              &
               biogenic,                                                      &
               snow_depth, land_frac,                                         &
               so4_aitken_work, so4_accu_work,                                &
               so4_diss_work, aged_bmass_work, cloud_bmass_work,              &
               aged_ocff_work, cloud_ocff_work, nitr_acc_work,                &
               nitr_diss_work, aerosol_work,                                  &
               n_arcl_species, n_arcl_compnts, i_arcl_compnts, arcl,          &
   ! Other variables for mphys
               ls_rain, ls_snow, ls_graup,                                    &
               ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,                  &
               n_drop_pot, n_drop_3d,                                         &
               rhc_row_length, rhc_rows,                                      &
               rhodz_dry, rhodz_moist,                                        &
   ! variables for Jules rainfrac code
               ls_rainfrac, land_points, land_index,                          &
   ! COSP variables
               l_cosp_lsp,                                                    &
   ! For subgrid orographic rain enhancement
               hmteff, zb,                                                    &
   ! Variables for stochastic physics random parameters
               m_ci,                                                          &
               error_code )

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)


  ! Copy the data of the output fields using tdims bounds to fields
  ! using tdims_s bounds as returned by ls_ppn.

!$OMP DO SCHEDULE(STATIC)                   
  DO k = 1,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        aerosol(i,j,k) = aerosol_work(i,j,k)
      END DO
    END DO
  END DO
!$OMP END DO NOWAIT


  ! Calculate increment fields
!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end

    DO j = tdims%j_start, tdims%j_end

      DO i = tdims%i_start, tdims%i_end

        t_inc(i,j,k) = t_work(i,j,k) - t_n(i,j,k) + t_inc(i,j,k)
        q_inc(i,j,k) = q_work(i,j,k) - q_n(i,j,k) + q_inc(i,j,k)
        qcl_inc(i,j,k) = qcl_work(i,j,k) - qcl_n(i,j,k)                       &
                       + qcl_inc(i,j,k)
        qcf_inc(i,j,k) = qcf_work(i,j,k) - qcf_n(i,j,k)                       &
                       + qcf_inc(i,j,k)
        cf_inc(i,j,k) = cf_work(i,j,k)   - cf_n(i,j,k)                        &
                       + cf_inc(i,j,k)
        cfl_inc(i,j,k) = cfl_work(i,j,k)   - cfl_n(i,j,k)                     &
                       + cfl_inc(i,j,k)
        cff_inc(i,j,k) = cff_work(i,j,k)   - cff_n(i,j,k)                     &
                       + cff_inc(i,j,k)

      END DO

    END DO

  END DO
!$OMP END DO

!$OMP END PARALLEL

  ! Save microphys tendencies to physics_tendencies_mod module
  ! When l_fix_lsp_incs_to_spt is removed from code then this
  ! whole block can be removed.
  IF (l_retain_mic_tendencies .AND. .NOT. l_fix_lsp_incs_to_spt) THEN

    DO k = 1,tdims%k_end
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          dt_mic(i,j,k) = t_work(i,j,k) - t_n(i,j,k)
        END DO
      END DO
    END DO

    DO k = 1,tdims%k_end
      DO j = tdims%j_start,tdims%j_end
        DO i = tdims%i_start,tdims%i_end
          dq_mic(i,j,k) = q_work(i,j,k) - q_n(i,j,k)
        END DO
      END DO
    END DO
  END IF !end if over l_retain_mic_tendencies

  IF (l_subgrid_qcl_mp) THEN

    ! Generate mixed phase cloud by turbulent processes

    CALL mphys_turb_gen_mixed_phase(q_work, t_work, qcl_work, qcf_work,   &
                                    q_inc, qcl_inc, cfl_inc, cf_inc,      &
                                    t_inc, dqcl_mp,  bl_levels,           &
                                    bl_w_var,                             &
                                    cff_work, cfl_work, cf_work, q_n,     &
                                    cfl_n, cf_n, p_layer_centres,         &
                                    rhodz_dry, rhodz_moist, deltaz,       &
                                    qcl_mpt, tau_d, inv_prt, disprate,    &
                                    inv_mt, si_avg, dcfl_mp, sigma2_s )

  END IF ! l_subgrid_qcl_mp

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i, j, k)

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, bl_levels

    DO j = tdims%j_start, tdims%j_end

      DO i = tdims%i_start, tdims%i_end

            ! Calculate and save microphys increments to TL and QW,
            ! required in boundary layer.

        micro_tends(i,j,1,k) = t_inc(i,j,k)                                   &
                              - (      lc*qcl_inc(i,j,k) )/cp                 &
                              - ( (lc+lf)*qcf_inc(i,j,k) )/cp

        micro_tends(i,j,2,k) = q_inc(i,j,k) + qcl_inc(i,j,k)                  &
                                            + qcf_inc(i,j,k)

        micro_tends(i,j,1,k) = micro_tends(i,j,1,k)/timestep
        micro_tends(i,j,2,k) = micro_tends(i,j,2,k)/timestep

      END DO

    END DO

  END DO
!$OMP END DO

!$OMP END PARALLEL

  ! Call pc2 erosion to apply the mp qcl increments

  IF (l_subgrid_qcl_mp .AND. l_subgrid_cfl_mp_by_erosion) THEN

    IF( model_type /= mt_single_column) THEN
      IF ( sf(302,4) ) THEN
        ALLOCATE ( cfl_copy( tdims%i_start : tdims%i_end,                     &
                             tdims%j_start : tdims%j_end,                     &
                                       1 : tdims%k_end ) )

        ALLOCATE ( qcl_copy( tdims%i_start : tdims%i_end,                     &
                             tdims%j_start : tdims%j_end,                     &
                                         1 : tdims%k_end ) )

        DO k = 1, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              cfl_copy(i,j,k) = cfl_work(i,j,k)
              qcl_copy(i,j,k) = qcl_work(i,j,k)
            END DO
          END DO
        END DO
      END IF ! sf
    END IF ! not SCM

    ! DEPENDS ON: pc2_turbulence_ctl 
    CALL pc2_turbulence_ctl (                                                 &
    ! Primary fields passed in, unchanged on exit
         T_work, q_work, qcl_work, cf_work, cfl_work, cff_work,               &
         p_layer_centres(1,1,1),                                              &
         dqcl_mp,                                                             &
    ! diagnostic info
         stashwork4,                                                          &
    ! SCM diagnostics switches (dummy in full UM)
         nSCMDpkgs, L_SCMDiags,                                               &
    ! Increment fields passed in, updated on exit
         t_inc, q_inc, qcl_inc, cf_inc, cfl_inc, l_subgrid_qcl_mp )

    IF( model_type /= mt_single_column) THEN
      IF ( sf(302,4) ) THEN
        DO k = 1, tdims%k_end
          DO j = tdims%j_start, tdims%j_end
            DO i = tdims%i_start, tdims%i_end
              cfl_incr_mp(i,j,k)     = cfl_work(i,j,k) - cfl_copy(i,j,k)
              qcl_incr_mp_pc2(i,j,k) = qcl_work(i,j,k) - qcl_copy(i,j,k)
            END DO
          END DO
        END DO

        DEALLOCATE ( qcl_copy )
        DEALLOCATE ( cfl_copy )
      END IF ! sf
    END IF ! not SCM
  END IF  ! l_subgrid_qcl_mp etc

  IF (l_mcr_qcf2)                                                             &
                     ! Second cloud ice variable in use
    qcf2_inc(:,:,:) = qcf2_work(:,:,:) - qcf2_n(:,:,:)                        &
                    + qcf2_inc(:,:,:)

  IF (l_mcr_qrain) THEN
                      ! Prognostic rain in use
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                      &
!$OMP SHARED( tdims, qrain_inc, qrain_work, qrain_n )                 &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qrain_inc(i,j,k) = qrain_work(i,j,k) - qrain_n(i,j,k)               &
                           + qrain_inc(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF  ! l_mcr_qrain

  IF (l_mcr_qgraup) THEN
                       ! Prognostic graupel in use
!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, qgraup_inc, qgraup_work, qgraup_n )               &
!$OMP PRIVATE( i, j, k )
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          qgraup_inc(i,j,k) = qgraup_work(i,j,k) - qgraup_n(i,j,k)            &
                            + qgraup_inc(i,j,k)
        END DO
      END DO
    END DO
!$OMP END PARALLEL DO
  END IF

  ! Call timer for LS rain code
  IF (ltimer) CALL timer ('AP1M LS Rain ',6)
END IF ! on error code zero

!-----------------------------------------------------------------------
! Section on thunderstorm electrification
! - Only set to work when prognostic graupel is in use -
!-----------------------------------------------------------------------

IF (l_use_electric .AND. l_mcr_qgraup) THEN

  CALL electric_main( qcf_work, qcf2_work, qgraup_work, rhodz_dry,             &
                      rhodz_moist, t_n, w, at_extremity, stashwork21,          &
                      flash_pot(:, :, 1 : tdims%k_end ) )

END IF

! ----------------------------------------------------------------------
! Section LSP.1.1 Tracer scavenging
! ----------------------------------------------------------------------

! below cloud scavenging of dust

! Scavenging of mineral dust
! (at present this is simplistic scheme based on 4.5 code)

!Initialise tracer scavenging to zero. Moved from ls_ppn, where
! it was being initialised, but not used.

!$OMP PARALLEL DO DEFAULT(NONE) SCHEDULE(STATIC)                       &
!$OMP SHARED( tdims, lscav_so2, lscav_nh3, lscav_so4ait, lscav_so4acc, &
!$OMP         lscav_so4dis )                                           &
!$OMP PRIVATE( i, j )
DO j = tdims%j_start,tdims%j_end
  DO i = tdims%i_start,tdims%i_end
    lscav_so2(i,j)    = 0.0
    lscav_nh3(i,j)    = 0.0
    lscav_so4ait(i,j) = 0.0
    lscav_so4acc(i,j) = 0.0
    lscav_so4dis(i,j) = 0.0
  END DO
END DO
!$OMP END PARALLEL DO

IF (l_dust) THEN

  CALL mass_calc( row_length, rows,                                           &
                  r_rho_levels(1:row_length,1:rows,1:model_levels),           &
                  r_theta_levels(1:row_length,1:rows,0:model_levels),         &
                  timestep, rho_r2(1:row_length,1:rows,1:model_levels),       &
                  q_n, qcl_n, qcf_n, dm )


  !   Put dust arrays together for simplicity

!$OMP PARALLEL DEFAULT(NONE)                                           &
!$OMP SHARED( tdims, ls_rain3d, krain_Dust, dust_all, ndiv,            &
!$OMP         lscav_dust_all, timestep, ksnow_dust, ls_snow3d, dm,     &
!$OMP         dust_div1, dust_div2, dust_div3, dust_div4, dust_div5,   &
!$OMP         dust_div6, l_twobin_dust )                               &
!$OMP PRIVATE( i, j, k, rainrate, snowrate, rate, delta_dust, idiv, jj,&
!$OMP          omp_block)

! Previous initialisation of omp_block will only be valid if not compiled
! with OpenMP, otherwise can have a private version to avoid race condition
!$  omp_block=CEILING(REAL(((tdims%j_end-tdims%j_start)+1))/                   &
!$    omp_get_num_threads())

  DO idiv = 1, ndiv

!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end

      DO i = tdims%i_start, tdims%i_end

        lscav_dust_all( i, j, idiv ) = 0.0

      END DO

    END DO
!$OMP END DO NOWAIT

  END DO ! idiv

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end

    DO j = tdims%j_start, tdims%j_end

      DO i = tdims%i_start, tdims%i_end

        dust_all(i,j,k,1)=dust_div1(i,j,k)
        dust_all(i,j,k,2)=dust_div2(i,j,k)
        IF (.NOT. l_twobin_dust) THEN
          dust_all(i,j,k,3)=dust_div3(i,j,k)
          dust_all(i,j,k,4)=dust_div4(i,j,k)
          dust_all(i,j,k,5)=dust_div5(i,j,k)
          dust_all(i,j,k,6)=dust_div6(i,j,k)
        END IF

      END DO

    END DO

  END DO ! k
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
  DO jj = tdims%j_start, tdims%j_end, omp_block
    DO k = tdims%k_end, 1, -1

      DO j = jj, MIN(jj+omp_block-1,tdims%j_end)

        DO i = tdims%i_start, tdims%i_end

          !            deal with possible -ive precip

          rainrate = MAX( ls_rain3d(i,j,k), 0.0)
          snowrate = MAX( ls_snow3d(i,j,k), 0.0)

          !              calc proportion of dust mixing ratio scavenged
          !CDIR UNROLL=NDIV

          DO idiv = 1, ndiv

            rate = ( rainrate * krain_dust(idiv) +                            &
                snowrate * ksnow_dust(idiv) ) *                               &
                3600.0 * timestep

            delta_dust=dust_all(i,j,k,idiv)*(1.0-1.0/(1.0+ rate ))
            !                calc mass of dust removed
            lscav_dust_all(i,j,idiv)=lscav_dust_all(i,j,idiv) +               &
            delta_dust * dm (i,j,k)
            !                decrement mixing ratio
            dust_all(i,j,k,idiv)=dust_all(i,j,k,idiv) - delta_dust

          END DO ! idiv

        END DO ! i

      END DO ! j

    END DO ! k
  END DO ! jj
!$OMP END DO

  !       put newly calculated values into main arrays

!$OMP DO SCHEDULE(STATIC)
  DO k = 1, tdims%k_end

    DO j = tdims%j_start, tdims%j_end

      DO i = tdims%i_start, tdims%i_end

        dust_div1(i,j,k)=dust_all(i,j,k,1)
        dust_div2(i,j,k)=dust_all(i,j,k,2)
        IF (.NOT. l_twobin_dust) THEN
          dust_div3(i,j,k)=dust_all(i,j,k,3)
          dust_div4(i,j,k)=dust_all(i,j,k,4)
          dust_div5(i,j,k)=dust_all(i,j,k,5)
          dust_div6(i,j,k)=dust_all(i,j,k,6)
        END IF

      END DO

    END DO

  END DO
!$OMP END DO
!$OMP END PARALLEL

END IF ! l_dust


IF (l_sulpc_so2) THEN

  CALL rainout_intctl(row_length, rows,                                       &
               off_x, off_y,                                                  &
               halo_i, halo_j,                                                &
               rho_r2, q_n,                                                   &
               qcf_work, qcl_work,                                            &
               qcf_n, qcl_n,                                                  &
               ls_rain3d, ls_snow3d,                                          &
               timestep,                                                      &
               so4_diss(:,:,1:), so4_accu(:,:,1:),                            &
               rnout_tracer)

  CALL sl3dwash(row_length, rows,                                             &
                off_x, off_y,                                                 &
                halo_i, halo_j,                                               &
                timestep,                                                     &
                rho_r2,                                                       &
                q_n, qcl_n, qcf_n, so2(:,:,1:),                               &
                ls_rain3d,                                                    &
                lscav_tr)

END IF ! On test for sulphur cycle

IF ( (l_sulpc_so2 .OR. l_nitrate) .AND. l_sulpc_nh3) THEN

  CALL nh3dwash(row_length, rows,                                             &
                off_x, off_y,                                                 &
                halo_i, halo_j,                                               &
                timestep,                                                     &
                rho_r2,                                                       &
                q_n, qcl_n, qcf_n, nh3(:,:,1:),                               &
                ls_rain3d,                                                    &
                lscav_nh3)
END IF ! On test for (L_sulpc_SO2 .OR. L_nitrate) .AND. L_sulpc_NH3

IF (l_soot) THEN  ! IF soot modelling is included

  CALL rainout_intctl(row_length, rows,                                       &
               off_x, off_y,                                                  &
               halo_i, halo_j,                                                &
               rho_r2, q_n,                                                   &
               qcf_work, qcl_work,                                            &
               qcf_n, qcl_n,                                                  &
               ls_rain3d, ls_snow3d,                                          &
               timestep,                                                      &
               cloud_soot(:,:,1:), aged_soot(:,:,1:),                         &
               rnout_soot)


  !       LS washout of soot was neglected at 4.5, although code was
  !       allegedly written to perform the calculation. We can't use
  !       the same treatment for soot as we use for the S cycle.
  !       For now, we will again neglect this process, but note that
  !       we may wish to include it in the future.

  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_soot(i,j)=0.0
    END DO
  END DO

END IF  ! l_soot

IF (l_biomass) THEN  ! IF biomass modelling is included

  CALL rainout_intctl(row_length, rows,                                       &
               off_x, off_y,                                                  &
               halo_i, halo_j,                                                &
               rho_r2, q_n,                                                   &
               qcf_work, qcl_work,                                            &
               qcf_n, qcl_n,                                                  &
               ls_rain3d, ls_snow3d,                                          &
               timestep,                                                      &
               cloud_bmass(:,:,1:), aged_bmass(:,:,1:),                       &
               rnout_bmass)


  !       As with soot, LS washout of biomass smoke aerosol is
  !       currently neglected. Again, though, we may wish to include
  !       it in the future.

  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_bmass(i,j)=0.0
    END DO
  END DO

END IF  ! l_biomass

IF (l_ocff) THEN  ! IF fossil-fuel org carb modelling is included

  CALL rainout_intctl(row_length, rows,                                       &
               off_x, off_y,                                                  &
               halo_i, halo_j,                                                &
               rho_r2, q_n,                                                   &
               qcf_work, qcl_work,                                            &
               qcf_n, qcl_n,                                                  &
               ls_rain3d, ls_snow3d,                                          &
               timestep,                                                      &
               cloud_ocff(:,:,1:), aged_ocff(:,:,1:),                         &
               rnout_ocff)


  !       As with soot and biomass, LS washout of fossil-fuel organic
  !       carbon aerosol is currently neglected. Again, though, we may
  !       wish to include it in the future.

  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_ocff(i,j)=0.0
    END DO
  END DO

END IF  ! l_ocff

IF (l_nitrate) THEN  ! IF ammonium nitrate modelling is included

  CALL rainout_intctl(row_length, rows,                                       &
               off_x, off_y,                                                  &
               halo_i, halo_j,                                                &
               rho_r2, q_n,                                                   &
               qcf_work, qcl_work,                                            &
               qcf_n, qcl_n,                                                  &
               ls_rain3d, ls_snow3d,                                          &
               timestep,                                                      &
               nitr_diss(:,:,1:), nitr_acc(:,:,1:),                           &
               rnout_nitrate)

  !       As with soot, biomass, and OCFF, LS washout of ammonium
  !       nitrate is currently neglected. Again, though, we may
  !       wish to include it in the future.

  DO j = 1, lspice_dim2
    DO i = 1, lspice_dim1
      lscav_nitrate(i,j)=0.0
    END DO
  END DO

END IF  ! l_nitrate

! Save microphys tendencies to physics_tendencies_mod module
IF (l_retain_mic_tendencies .AND. l_fix_lsp_incs_to_spt) THEN

  DO k = 1,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        dt_mic(i,j,k) = t_work(i,j,k) - t_n(i,j,k)
      END DO
    END DO
  END DO

  DO k = 1,tdims%k_end
    DO j = tdims%j_start,tdims%j_end
      DO i = tdims%i_start,tdims%i_end
        dq_mic(i,j,k) = q_work(i,j,k) - q_n(i,j,k)
      END DO
    END DO
  END DO
END IF !end if over l_retain_mic_tendencies

! ----------------------------------------------------------------------
! Section LSP.2 Output Diagnostics
! ----------------------------------------------------------------------

! Save aggregate fraction for COSP
IF (l_cosp_lsp) cosp_frac_agg = frac_agg

! Check that microphysics diagnostics requested this timestep
SELECT CASE (model_type)

CASE DEFAULT
  IF (error_code  ==  0 .AND. sf(0,4)) THEN

    CALL diagnostics_lsrain(                                                  &
                     lspice_dim1,lspice_dim2,lspice_dim3,                     &
                     p_layer_centres, deltaz, rhodz_dry,                      &
                     t_work, q_work, qcl_work, qcf_work,                      &
                     qrain_work, qgraup_work, qcf2_work,                      &
                     cf_work, cfl_work, cff_work,                             &
                     t_n, q_n, qcl_n, qcf_n,                                  &
                     qrain_n, qgraup_n, qcf2_n,                               &
                     cf_n, cfl_n, cff_n,                                      &
                     ls_rain, ls_snow, ls_graup,                              &
                     ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,            &
                     rnout_tracer,lscav_dust_all,lscav_tr,                    &
                     lscav_nh3,                                               &
                     rnout_soot, lscav_soot,                                  &
                     rnout_bmass, lscav_bmass,                                &
                     rnout_ocff, lscav_ocff,                                  &
                     rnout_nitrate, lscav_nitrate,                            &
                     dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_r, dbz_l,             &
                     psdep,psaut,psacw,psacr,                                 &
                     psaci,psmlt,psmltevp,                                    &
                     praut,pracw,prevp,                                       &
                     pgaut,pgacw,pgacs,pgmlt,                                 &
                     pifrw,pifrr,piprm,pidep,piacw,                           &
                     piacr,pimlt,pimltevp,                                    &
                     pifall,psfall,prfall,pgfall,                             &
                     plset, plevpset, n_drop_pot, n_drop_3d,                  &
                     frac_agg, mphys_pts,                                     &
                     vm_cry, vm_agg, vtbranch_flag, vm_used,                  &
                     sfwater,sfrain,sfsnow,                                   &
                     qcl_mpt, tau_d, inv_prt, disprate,                       &
                     inv_mt, si_avg, dcfl_mp, sigma2_s,                       &
                     cfl_incr_mp, qcl_incr_mp_pc2,                            &
                     stashwork4                                               &
     )

  END IF   ! on error_code and sf(0,4)

  ! Copy 3D rain and snow into COSP input structure
  ! The mass of snow is already contained in the LS ice water variable
  ! (aggregates + crystals +[graupel]), so it does not need to be
  ! considered here.

  IF (L_cosp_lsp) THEN
    cosp_gbx%rain_ls(:,1:lspice_dim3) =                                       &
      RESHAPE(ls_rain3d, (/lspice_dim1*lspice_dim2, lspice_dim3/))
  END IF


CASE (mt_single_column)

  !-----------------------------------------------------------------
  !     SCM Large Scale Precip OR Increments Diagnostics Package
  !-----------------------------------------------------------------
  IF (L_SCMDiags(scmdiag_lsp) .OR. L_SCMDiags(scmdiag_incs)) THEN

    !       Stash 4,181
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = t_work(i,j,k) - t_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    CALL scmoutput(work_3d,'dt_lsr',                                          &
         'Temperature increment, large scale rain','K',                       &
         t_avg,d_all,default_streams,'',RoutineName)

    !       Stash 4,182
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = q_work(i,j,k) - q_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    CALL scmoutput(work_3d,'dq_lsr',                                          &
         'Specific humidity increment, large scale rain','kg/kg',             &
         t_avg,d_wet,default_streams,'',RoutineName)

    !       Stash 4,183
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = qcl_work(i,j,k) - qcl_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    CALL scmoutput(work_3d,'dqcl_lsr',                                        &
         'QCL increment, large scale rain','kg/kg',                           &
         t_avg,d_wet,default_streams,'',RoutineName)

    !       Stash 4,184
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = qcf_work(i,j,k) - qcf_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    CALL scmoutput(work_3d,'dqcf_lsr',                                        &
         'QCF increment, large scale rain','kg/kg',                           &
         t_avg,d_wet,default_streams,'',RoutineName)

    !       Stash 4,192
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = cf_work(i,j,k) - cf_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    CALL scmoutput(work_3d,'dbcf_lsr',                                        &
         'Bulk cloud frac increment, large scale rain','Fraction',            &
         t_avg,d_wet,default_streams,'',RoutineName)

    !       Stash 4,193
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = cfl_work(i,j,k) - cfl_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    CALL scmoutput(work_3d,'dcfl_lsr',                                        &
         'Liquid cloud frac increment, large scale rain','Fraction',          &
         t_avg,d_wet,default_streams,'',RoutineName)

    !       Stash 4,194
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          work_3d(i,j,k) = cff_work(i,j,k) - cff_n(i,j,k)
        END DO ! i
      END DO ! j
    END DO ! k

    CALL scmoutput(work_3d,'dcff_lsr',                                        &
         'Frozen cloud frac increment, large scale rain','Fraction',          &
         t_avg,d_wet,default_streams,'',RoutineName)

  END IF ! L_SCMDiags(SCMDiag_lsp) .OR. L_SCMDiags(SCMDiag_incs)

  !---------------------------------------------------------------------
  !     SCM Large Scale Precip Diagnostics Package
  !---------------------------------------------------------------------
  IF (L_SCMDiags(scmdiag_lsp)) THEN

    !       Stash 4,222
    CALL scmoutput(ls_rain3d,'ls_rain3d',                                     &
         'Rainfall rate out of model levels','kg/m2/s',                       &
         t_mult,d_wet,default_streams,'sec_day',RoutineName)

    !       Stash 4,223
    CALL scmoutput(ls_snow3d,'ls_snow3d',                                     &
         'Snowfall rate out of model levels','kg/m2/s',                       &
         t_mult,d_wet,default_streams,'sec_day',RoutineName)

    !       (copied from daglsran)
    !       Produce diagnostic 225 before 224 in order to save memory.
    !       Supercooled 3D rain content. It is equal to
    !       the 3D rainrate at T < 0 and equal to 0 at T > 0
    !       Alter the array LS_RAIN3D directly
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          IF (t_work(i,j,k) >= zerodegc) THEN
            ! Warm temperatures
            ls_rain3d(i,j,k) = 0.0
          END IF
        END DO ! i
      END DO ! j
    END DO ! k

    !       Stash 4,225
    CALL scmoutput(ls_rain3d,'ls_rain_scool',                                 &
         'Supercooled rain out of model levels','kg/m2/s',                    &
         t_avg,d_wet,default_streams,'',RoutineName)

    !       Supercooled liquid water content at TIME LEVEL N. It is equal
    !       to the liquid water content at T < 0 and equal to 0 at T > 0
    !       Use LS_RAIN3D as the array in order to save memory
    DO k = 1, model_levels
      DO j = 1, rows
        DO i = 1, row_length
          IF (t_work(i,j,k) < zerodegc) THEN
            ! Supercooled temperatures
            ! Use time level n fields in this diagnostic
            ls_rain3d(i,j,k) = qcl_n(i,j,k)
          ELSE
            ! Warm temperatures
            ls_rain3d(i,j,k) = 0.0
          END IF
        END DO ! i
      END DO ! j
    END DO ! k

    !       Stash 4,224
    CALL scmoutput(ls_rain3d,'ls_water_scool',                                &
         'Supercooled liquid water content','kg/kg',                          &
         t_avg,d_wet,default_streams,'',RoutineName)

    !       Stash 4,227
    CALL scmoutput(rainfrac3d,'rainfrac3d',                                   &
         'Rain fraction out of model levels','Fraction',                      &
         t_avg,d_wet,default_streams,'',RoutineName)

  END IF ! L_SCMDiags(SCMDiag_lsp)

END SELECT ! model_type

! REMINDER: DEALLOCATE FIELDS IN REVERSE ORDER OF ALLOCATION

DEALLOCATE ( sfsnow )        ! 56
DEALLOCATE ( sfrain )        ! 55
DEALLOCATE ( sfwater )       ! 54
DEALLOCATE ( dbz_l )         ! 53
DEALLOCATE ( dbz_r )         ! 52
DEALLOCATE ( dbz_i2 )        ! 51
DEALLOCATE ( dbz_i )         ! 50
DEALLOCATE ( dbz_g )         ! 49
DEALLOCATE ( dbz_tot )       ! 48
DEALLOCATE ( qcl_incr_mp_pc2)! 47
DEALLOCATE ( cfl_incr_mp )   ! 46
DEALLOCATE ( sigma2_s )      ! 45
DEALLOCATE ( dcfl_mp )       ! 44
DEALLOCATE ( si_avg )        ! 43
DEALLOCATE ( inv_mt )        ! 42
DEALLOCATE ( disprate )      ! 41
DEALLOCATE ( inv_prt )       ! 40
DEALLOCATE ( tau_d )         ! 39
DEALLOCATE ( qcl_mpt)        ! 38
DEALLOCATE ( vm_used )       ! 37
DEALLOCATE ( vm_agg )        ! 36
DEALLOCATE ( vm_cry )        ! 35
DEALLOCATE ( vtbranch_flag ) ! 34
DEALLOCATE ( plevpset )      ! 33
DEALLOCATE ( plset )         ! 32
DEALLOCATE ( pgfall )        ! 31
DEALLOCATE ( prfall )        ! 30
DEALLOCATE ( psfall )        ! 29
DEALLOCATE ( pifall )        ! 28
DEALLOCATE ( pgmlt )         ! 27
DEALLOCATE ( pgacs )         ! 26
DEALLOCATE ( pgacw )         ! 25
DEALLOCATE ( pgaut )         ! 24
DEALLOCATE ( prevp )         ! 23
DEALLOCATE ( pracw )         ! 22
DEALLOCATE ( praut )         ! 21
DEALLOCATE ( psaci )         ! 20
DEALLOCATE ( psaut )         ! 19
DEALLOCATE ( psmlt )         ! 18
DEALLOCATE ( pimlt )         ! 17
DEALLOCATE ( psmltevp )      ! 16
DEALLOCATE ( pimltevp )      ! 15
DEALLOCATE ( psacr )         ! 14
DEALLOCATE ( piacr )         ! 13
DEALLOCATE ( psacw )         ! 12
DEALLOCATE ( piacw )         ! 11
DEALLOCATE ( psdep )         ! 10
DEALLOCATE ( pidep )         ! 9
DEALLOCATE ( piprm )         ! 8
DEALLOCATE ( pifrr )         ! 7
DEALLOCATE ( pifrw )         ! 6
DEALLOCATE ( mphys_pts )     ! 5
DEALLOCATE ( frac_agg )      ! 4
DEALLOCATE ( qgraup_work )   ! 3
DEALLOCATE ( qrain_work )    ! 2
DEALLOCATE ( qcf2_work )     ! 1

! end of routine microphys_ctl
IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE microphys_ctl
