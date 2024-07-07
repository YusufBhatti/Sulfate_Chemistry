! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Description:
!  Call the UK Chemistry and Aerosols submodel components.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
!  Method:
! 1) Interface with atmosphere model:
!   Prognostic fields are selected from D1 using the item codes
!   defined in the routine UKCA_SetD1Defs, these include tracers
!   from section 34 (UKCA).
!   Diagnostic fields are selected from D1 using a tag of 98, these
!   need to be set in the STASH panel of the gui.
! 2) UKCA routines are called depending on the scheme selected:
!    - Woodward dust scheme
!    - Emissions routine
!    - photolysis routine
!    - Chemistry control routine
! 3) Updated tracer arrays are returned to the D1 array
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA
!
!  Code Description:
!    Language:  Fortran 2003
!
! ----------------------------------------------------------------------
MODULE ukca_main1_mod


 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

PRIVATE
PUBLIC :: ukca_main1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'UKCA_MAIN1_MOD'

CONTAINS
!
! Subroutine Interface:
! ----------------------------------------------------------------------
SUBROUTINE ukca_main1(tracer_ukca_um, q_um)
! ----------------------------------------------------------------------
USE ukca_option_mod,        ONLY:                                       &
    l_ukca_set_trace_gases, l_ukca_chem,        l_ukca_trop,            &
    l_ukca_aerchem,         l_ukca_raq,         l_ukca_intdd,           &
    l_ukca_advh2o,          l_ukca_mode,        l_ukca_strattrop,       &
    l_ukca_strat,           l_ukca_stratcfc,    l_ukca_trophet,         &
    l_ukca_het_psc,         l_ukca_sa_clim,     l_ukca_use_background_aerosol, &
    l_ukca_nr_aqchem,       l_ukca_qch4inter,   l_ukca_ageair,          &
    l_ukca_prescribech4,    l_ukca_arg_act,     jpctr,                  &
    jpspec,                 jppj,               jpdd,                   &
    jpdw,                   ukca_int_method,    l_ukca,                 &
    i_ukca_photol,          l_ukca_offline,     l_ukca_offline_be,      &
    l_ukca_chem_plev,       l_ukca_asad_plev,   chem_timestep,          &
    l_ukca_classic_hetchem, l_ukca_raqaero,     mode_parfrac,           &
    l_ukca_asad_columns,    l_ukca_limit_nat

USE ukca_um_interf_mod,     ONLY:                                       &
    getd1flds,              putd1flds,            ukca_pr_inputs,       &
    ukca_um_d1_initialise,  ukca_um_dealloc,      conv_cloud_base,      &
    conv_cloud_top,         kent,                 kent_dsc,             &
    land_sea_mask,          ntp_data,             all_ntp,              &
    emission1,              all_emissions,        dust_flux,            &
    theta,                  soil_moisture_layer1, q,                    &
    qcf,                    conv_cloud_lwp,       Tstar,                &
    zbl,                    Rough_length,         seaice_frac,          &
    um_ozone,               chloro_sea,           so4_aitken,           &
    so4_accum,              soot_fresh    ,       soot_aged,            &
    ocff_fresh    ,         ocff_aged     ,       dms_sea_conc,         &
    vertvel,                conv_cloud_amount,    Frac_types,           &
    laift_lp,               canhtft_lp,           canwctile_lp,         &
    Tstar_tile,             z0tile_lp,            Snow_tile,            &
    rho_r2,                 qcl,                  exner_rho_levels,     &
    area_cloud_fraction,    cloud_frac,           cloud_liq_frac,       &
    biogenic      ,         fland,                land_albedo,          &
    exner_theta_levels,     p_rho_levels,         p_theta_levels,       &
    pstar,                  net_surf_SW,          tot_surf_SW,          &
    sea_salt_film ,         sea_salt_jet  ,       sulphate_od,          &
    rhokh_mix,              dtrdz_charney_grid,   we_lim,               &
    t_frac,                 zrzi,                 we_lim_dsc,           &
    t_frac_dsc,             zrzi_dsc,             zhsc,                 &
    U_scalar_10m,           surf_hf,              stcon,                &
    u_s,                    bl_tke,               cloud_liq_water,      &
    ls_rain3d,              ls_snow3d,            ls_ppn_frac,          &
    ice_melt,               snow_melt,            rim_cry,              &
    rim_agg,                autoconv,             accretion,            &
    conv_rain3d,            conv_snow3d,          ch4_wetl_emiss,       &
    pv_on_theta_mlevs,      tropopause_height,    trmol_post_atmstep,   &
    mode_diags,             strat_fluxdiags,      totnodens,            &
    rel_humid_frac,         qsvp,                 qsmr

USE dust_parameters_mod,    ONLY: ndiv
USE tuning_segments_mod,    ONLY: ukca_mode_seg_size
USE conversions_mod,        ONLY: pi, isec_per_day

USE asad_mod,               ONLY:                                       &
                            asad_mod_final,   method,     nnaf,         &
                            interval,         advt,       fch4,         &
                            fco2,             fh2,        fn2,          &
                            fo2,              spt,        ntrkx,        &
                            jpspt

USE c_sulchm_mod,           ONLY: rad_ait, rad_acc, chi, sigma
USE ukca_tracer_vars
USE ukca_cspecies
USE ukca_constants,         ONLY: zboltz, c_co2, c_h2, c_n2, c_o2,      &
                                  c_ch4, rhosul, avc, imdi
USE ukca_tropopause,        ONLY:                                       &
    p_tropopause,           theta_trop,            pv_trop,             &
    tropopause_level,       l_troposphere,         ukca_calc_tropopause

USE ukca_mode_setup,        ONLY: nmodes, ncp, mode, component
USE ukca_mode_verbose_mod,  ONLY: glob_verbose
USE ukca_setup_indices
USE ukca_chem1_dat
USE ukca_chem_offline,      ONLY: o3_offline_diag, oh_offline_diag,     &
                                  no3_offline_diag, ho2_offline_diag
USE ukca_read_offline_oxidants_mod,                                     &
                            ONLY: ukca_read_offline_oxidants_ctl,       &
                                  ukca_offline_oxidants_diags
USE ukca_scavenging_diags_mod,                                          &
                            ONLY: item_no_plume,                        &
    ukca_plume_scav_diags_2d,     plume_scav_diag_2d,                   &
    plume_scav_diag_ping,         d_tracers,                            &
    msect38

USE ukca_d1_defs,           ONLY:                                       &
    n_chem_tracers,         n_aero_tracers,         em_chem_spec,       &
    L_ukca_stratflux,       n_strat_fluxdiags,      n_mode_diags,       &
    Nukca_D1items,          UkcaD1Codes,            istrat_first,       &
    n_boundary_vals,        lbc_mmr,                n_mode_tracers,     &
    n_use_emissions,        n_chem_emissions,       UKCA_sect,          &
    lbc_spec,               n_use_tracers,          nmax_mode_diags,    &
    imode_first,            item1_mode_diags,       nmax_strat_fluxdiags,&
    UKCA_diag_sect,         l_ukca_plume_diags

USE ukca_trace_gas_mixratio
USE ukca_transform_halogen_mod,                                         &
                            ONLY: ukca_transform_halogen
USE ukca_mode_tracer_maps_mod,                                          &
                            ONLY: ukca_aero_tracer_init

USE nlsizes_namelist_mod,   ONLY:                                       &
    global_row_length,      global_rows,          land_field,           &
    row_length,             rows,                 tr_ukca,              &
    len_tot,                model_levels,         n_obj_d1_max,         &
    bl_levels,              a_len_inthd,          a_len_cfi1,           &
    a_len_cfi2,             a_len_cfi3,           a_len_realhd,         &
    a_len1_levdepc,         a_len2_levdepc,       a_len1_coldepc,       &
    a_len2_coldepc,         a_len1_rowdepc,       a_len2_rowdepc,       &
    a_len1_flddepc,         a_len2_flddepc,       a_len_extcnst,        &
    a_len2_lookup,          st_levels,            sm_levels,            &
    tpps_ozone_levels,      theta_field_size,     tr_vars,              &
    tr_lbc_vars,            tr_lbc_ukca,          aocpl_row_length,     &
    aocpl_p_rows,           n_rows,               ozone_levels,         &
    len_fixhd,              len_dumphist,         len1_lookup,          &
    mpp_len1_lookup,        n_cca_lev,            ntiles,               &
    len1_lbc_comp_lookup,   tr_levels

USE asad_flux_dat,          ONLY: asad_load_default_fluxes

USE asad_chem_flux_diags,   ONLY:                                       &
    asad_allocate_chemdiag ,    asad_setstash_chemdiag,                 &
    asad_init_chemdiag,         asad_tendency_ste,                      &
    asad_mass_diagnostic,       asad_output_tracer,                     &
    asad_flux_put_stash,        calculate_STE,                          &
    calculate_tendency,         L_asad_use_chem_diags,                  &
    L_asad_use_STE,             L_asad_use_tendency,                    &
    L_asad_use_mass_diagnostic, L_asad_use_output_tracer,               &
    L_asad_use_trop_mask,       asad_tropospheric_mask

USE ukca_chem_schemes_mod,  ONLY: int_method_nr, int_method_be_explicit

USE ukca_diurnal_oxidant,   ONLY: ukca_int_cosz
USE ukca_emiss_ctl_mod,     ONLY: ukca_emiss_ctl
USE ukca_chemistry_ctl_be_mod, ONLY: ukca_chemistry_ctl_be

USE ukca_scenario_ctl_mod,  ONLY: ukca_scenario_ctl
USE ukca_set_array_bounds_mod, ONLY: ukca_set_array_bounds
USE ukca_chem_diags_allts_mod, ONLY: ukca_chem_diags_allts
USE ukca_chem_diags_mod,    ONLY: ukca_chem_diags
USE ukca_chemistry_ctl_mod, ONLY: ukca_chemistry_ctl
USE ukca_chemistry_ctl_col_mod, ONLY: ukca_chemistry_ctl_col
USE ukca_activate_mod,      ONLY: ukca_activate
USE ukca_aero_ctl_mod,      ONLY: ukca_aero_ctl
USE ukca_mode_diags_mod,    ONLY: ukca_mode_diags_alloc, ukca_mode_diags, &
                                  mdwat_diag, l_ukca_cmip6_diags,         &
                                  l_ukca_pm_diags
USE um_stashcode_mod,       ONLY: stashcode_so4_nuc_sol,          &
                                  stashcode_glomap_sec,           &
                                  stashcode_h2o_mmr,              &
                                  stashcode_pm_first,             &
                                  stashcode_pm_last      
USE ukca_age_air_mod,       ONLY: ukca_age_air
USE carbon_options_mod,     ONLY: l_co2_interactive
USE level_heights_mod,      ONLY: r_theta_levels, r_rho_levels,         &
                                eta_theta_levels
USE cderived_mod,           ONLY: delta_lambda, delta_phi, base_phi
USE planet_constants_mod,   ONLY: repsilon, g, two_omega, planet_radius
USE jules_surface_types_mod,ONLY: ntype, npft, soil
USE set_rad_steps_mod,      ONLY: l_rad_step_prog
USE turb_diff_mod,          ONLY: q_pos_method, qlimit,                 &
              l_qpos_diag_pr, qpos_diag_limit, q_pos_tracer_method
USE ereport_mod,            ONLY: ereport
USE Field_Types
USE UM_ParParams
USE lbc_mod
USE missing_data_mod,       ONLY: rmdi
USE dump_headers_mod,       ONLY: rh_z_top_theta, a_inthd, a_realhd
USE o3intp_mod,             ONLY: io3_3dspec,        io3_2dspec,        &
    io3_2dmasscon,                io3_trop_map,      io3_trop_map_masscon
USE cppxref_mod,            ONLY: ppxref_codelen,    ppxref_charlen
USE nlstgen_mod,            ONLY: steps_per_periodim, secs_per_periodim
USE submodel_mod,           ONLY:                                       &
    atmos_sm,               submodel_for_sm,          atmos_im
USE stash_array_mod,        ONLY:                                       &
    len_stlist, stash_maxlen,   stindex, stlist,      num_stash_levels, &
    stash_levels,               si,                   sf
USE stparam_mod,            ONLY: st_d1pos, st_macrotag
USE nlstcall_mod,           ONLY: ltimer
USE Control_Max_Sizes,      ONLY: max_req_thpv_levs

! required for clear-sky RH calculation
USE cloud_inputs_mod,       ONLY: rhcrit
USE lsp_subgrid_mod,        ONLY: lsp_qclear

USE model_time_mod,         ONLY:                                       &
    i_day,                  i_day_number,                  i_hour,      &
    i_minute,               i_month,                       i_year,      &
    previous_time,          secs_per_stepim,               stepim

USE land_soil_dimensions_mod, ONLY: land_points, land_index

USE atm_fields_bounds_mod,  ONLY: tdims, tdims_l, tdims_s, udims, wdims

USE umPrintMgr
USE UM_ParVars
USE UM_ParCore,   ONLY: mype, nproc

USE ukca_ntp_mod, ONLY: ntp_init, name2ntpindex

USE trignometric_mod,     ONLY: true_longitude,                   &
                                true_latitude,                    &
                                sin_theta_longitude,              &
                                sin_theta_latitude,               &
                                sin_v_latitude,                   &
                                cos_v_latitude,                   &
                                FV_cos_theta_latitude,            &
                                cos_theta_longitude,              &
                                tan_theta_latitude

USE rad_input_mod,        ONLY: a_sw_radstep_prog, l_use_biogenic,&
                                l_use_seasalt_direct,             &
                                l_use_seasalt_indirect,           &
                                l_use_arclsulp
USE run_aerosol_mod,       ONLY: l_sulpc_so2, l_soot, l_ocff,     &
                                 l_use_seasalt_sulpc
USE mphys_inputs_mod,     ONLY: l_use_seasalt_autoconv

USE dyn_coriolis_mod,       ONLY: f3_at_u
USE ukca_photo_scheme_mod,  ONLY: i_ukca_fastjx
USE spec_sw_lw,             ONLY: sw_spectrum
USE ukca_all_tracers_copy_mod, ONLY: ukca_all_tracers_copy_in,          &
                               ukca_all_tracers_copy_out,     all_tracers

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

USE ukca_iniasad_mod,       ONLY: ukca_iniasad
USE ukca_read_aerosol_mod,  ONLY: ukca_read_aerosol
USE ukca_setd1defs_mod,     ONLY: ukca_setd1defs
USE ukca_solang_mod,        ONLY: ukca_solang

USE tilepts_mod,            ONLY: tilepts
USE ukca_plev_diags_mod,    ONLY: ukca_plev_diags

USE solpos_mod,             ONLY: solpos

!Redirect routine names to avoid clash with existing qsat routines
USE qsat_mod, ONLY: qsat_new         => qsat,                           &
                    qsat_wat_mix_new => qsat_wat_mix,                   &
                    l_new_qsat_chem_aero !Currently defaults to FALSE

IMPLICIT NONE

! ----------------------------------------------------------------------
! Declarations
! ----------------------------------------------------------------------

! UKCA Tracers coming from the UM code via subroutine interface
REAL, INTENT(INOUT) :: tracer_ukca_um                           &
              (tdims_s%i_start:tdims_s%i_end,                   &
               tdims_s%j_start:tdims_s%j_end,                   &
               tdims_s%k_start:tdims_s%k_end,1:tr_ukca)
! Specific humidity coming from the UM code via subroutine interface
REAL, INTENT(INOUT) :: q_um                                     &
               (tdims_l%i_start:tdims_l%i_end,                  &
                tdims_l%j_start:tdims_l%j_end,                  &
                tdims_l%k_start:tdims_l%k_end)

! Local scalars
INTEGER    :: section           ! stash section
INTEGER    :: item              ! stash item
INTEGER    :: item_first        ! first stash item
INTEGER    :: item_last         ! last stash item
INTEGER    :: i,i1,i2           ! loop variables
INTEGER    :: ii                ! loop variables
INTEGER    :: j,j1,j2           ! loop variables
INTEGER    :: k,l,n             ! loop variables
INTEGER    :: kk                ! loop counter
INTEGER    :: imode             ! loop counter for modes
INTEGER    :: icp               ! loop counter for components
INTEGER    :: m                 ! loop counter
INTEGER    :: m_atm_modl        ! sub model
INTEGER    :: timestep_number   ! no. of atmos timesteps since bas
INTEGER    :: nd_o3             ! size of um ozone array
INTEGER    :: irad_prev         ! Indices for controlling
INTEGER    :: irad_next         ! albedo calculation
INTEGER    :: days_in_period    ! Real days covered by secs_per_period
INTEGER    :: asteps_per_day    ! Atmos steps per day
INTEGER    :: A_Steps_per_hr    ! Atmos steps per hour
INTEGER    :: Radsteps_per_day  ! Radiation TS per day
INTEGER    :: icnt              ! counter
INTEGER    :: icode=0           ! local error status
INTEGER    :: im_index          ! internal model index
INTEGER, SAVE :: n_day          ! test whether new day has been reached
INTEGER, SAVE :: wetox_in_aer   ! set for wet oxidation in MODE (=1)
                                ! or UKCA (=0)
INTEGER, SAVE :: uph2so4inaer   ! update H2SO4 tracer in MODE (=1)
                                ! or UKCA (=0) depending on chemistry
INTEGER, SAVE :: k_be_top       ! top level for offline oxidants (BE)

REAL       :: timestep                  ! atmosphere model TS
REAL       :: r_minute                  ! real equiv of i_minute
REAL       :: min_surf_albedo=0.1       ! fill-in value
REAL       :: fx                        ! interpolation fraction
REAL       :: eq_time                   ! Equation of time
REAL       :: Sindec                    ! Solar declination
REAL       :: sindec_obs, eqt_obs       ! dummy variables here
REAL       :: secondssincemidnight      ! day time
REAL       :: scs                       ! solar constant scaling factor
REAL       :: ssmn_incr                 ! sec. since midnight incr. counter

CHARACTER(LEN=10)  :: prods(2)          ! Products

LOGICAL       :: L_moses_ii             ! T if MOSES_II in use
LOGICAL       :: l_classic_aerosols     ! True if CLASSIC aerosols are modelled

! Local arrays
! Cache blocking feature : here sized to a maximum but within aero_ctl will be
! set to match the case in progress
INTEGER :: stride_seg   ! row_length*rows for striding through by column
INTEGER :: nseg         ! the number of segments on this MPI
INTEGER :: nbox_s(row_length*rows) ! the number of boxes on each segment
INTEGER :: ncol_s(row_length*rows) ! number of columns on each segment
INTEGER :: seg_rem                 ! used to establish if more segments needed
INTEGER :: lbase(row_length*rows)  ! start point of segment in the 3d array
INTEGER :: nbox_max                ! largest chunk has highest number of boxes
INTEGER :: ik,lb                   ! local loop iterators

! MPP-related Arrays
! Table of number of points on a row
INTEGER :: g_row_length(0:nproc-1)
! Table number of rows in theta field
INTEGER :: g_rows(0:nproc-1)
! Tile index etc
INTEGER :: tile_index(land_points,ntype) !
INTEGER :: tile_pts(ntype)   ! No of tile p

REAL, ALLOCATABLE :: qsatmr(:,:,:)        ! qsat mixing ratio
REAL, ALLOCATABLE :: qsatmr_wat(:,:,:)  ! qsatmr w.r.t. wat
REAL, ALLOCATABLE :: rel_humid_frac_clr(:,:,:)
REAL, ALLOCATABLE :: rel_humid_frac_clr1D(:)
REAL, ALLOCATABLE :: q_clr1D(:)
REAL, ALLOCATABLE :: q1d(:)
REAL, ALLOCATABLE :: qsatmr1d(:)
REAL, ALLOCATABLE :: qsatmr_wat1d(:)
REAL, ALLOCATABLE :: qcf1d(:)
REAL, ALLOCATABLE :: cloud_liq_frac1d(:)
REAL, ALLOCATABLE :: cloud_frac1d(:)
REAL, ALLOCATABLE :: rhcrit1d(:)
REAL, ALLOCATABLE :: ls_ppn3d(:,:,:)
REAL, ALLOCATABLE :: conv_ppn3d(:,:,:)
REAL, ALLOCATABLE :: delso2_wet_h2o2(:,:,:)
REAL, ALLOCATABLE :: delso2_wet_o3(:,:,:)
REAL, ALLOCATABLE :: delh2so4_chem(:,:,:)  ! Rate of net chem Pdn H2SO4
REAL, ALLOCATABLE :: delso2_drydep(:,:,:)
REAL, ALLOCATABLE :: delso2_wetdep(:,:,:)
REAL, ALLOCATABLE :: um_ozone3d(:,:,:) ! after O3_to_3D
REAL, ALLOCATABLE :: um_ozone1d(:) ! ip to O3_to_3D
REAL, ALLOCATABLE   :: Tile_Frac(:,:)
REAL, ALLOCATABLE, SAVE :: land_albedo_all(:,:,:)  ! all r
REAL, ALLOCATABLE   :: surf_albedo(:,:)  ! interpolated
REAL, ALLOCATABLE   :: climoz2d(:,:) !climatological O3
REAL, ALLOCATABLE   :: theta_latitude(:,:) ! gridbox lat
REAL, ALLOCATABLE   :: v_latitude(:,:)     ! boundary lat
REAL, ALLOCATABLE   :: sinv_latitude(:,:)  ! sin(boundary lat)
REAL, ALLOCATABLE   :: cos_zenith_angle(:,:)
REAL, ALLOCATABLE, SAVE   :: int_zenith_angle(:,:) ! integral over sza
REAL, ALLOCATABLE, SAVE             :: daylength(:,:)         ! day length
REAL, ALLOCATABLE   :: rb_dust_div1(:,:)

! solid angle of grid cell saved for mass calculation
REAL, SAVE, ALLOCATABLE :: mass (:,:,:)
REAL, SAVE, ALLOCATABLE :: solid_angle(:,:)   ! solid angle
REAL, SAVE, ALLOCATABLE :: z_half(:,:,:)
REAL, SAVE, ALLOCATABLE :: z_half_alllevs(:,:,:)
REAL, SAVE, ALLOCATABLE :: volume(:,:,:)   ! gridbox volume on theta points
REAL, SAVE, ALLOCATABLE :: area(:,:,:)     ! gridbox area
REAL, SAVE, ALLOCATABLE :: delta_r(:,:,:)  ! delta radius
REAL, SAVE, ALLOCATABLE :: interf_z(:,:,:) ! altitude of interfaces (in
                                           ! metres) above ground level
REAL, SAVE, ALLOCATABLE :: f3u_by_omega(:,:)  ! Variable representing
! f3_at_u/two_omega. Needs to be calculated seperately for ND and EG as
! indexing is different between the two grids.

! diagnostics from chemistry_ctl
! z dimension is always 1:model_levels - this is consistent
! with STASH metadata

! Nitric acid trihydrate (kg(nat)/kg(air))
REAL :: nat_psc(row_length,rows,model_levels)

! Trop CH4 burden in moles
REAL :: trop_ch4_mol(row_length,rows,model_levels)

! Trop O3 burden in moles
REAL :: trop_o3_mol(row_length,rows,model_levels)

! Trop OH burden in moles
REAL :: trop_oh_mol(row_length,rows,model_levels)

! Strat CH4 burden in moles
REAL :: strat_ch4_mol(row_length,rows,model_levels)

! Stratospheric CH4 loss (Moles/s)
REAL :: strat_ch4loss(row_length,rows,model_levels)

! Sin(latitude)
REAL :: sinlat(udims%i_start:udims%i_end, udims%j_start:udims%j_end)

! Atmospheric Burden of CH4 in moles
REAL :: atm_ch4_mol(row_length,rows,model_levels)

! Atmospheric Burden of CO in moles
REAL :: atm_co_mol(row_length,rows,model_levels)

! Atmospheric Burden of Nitrous Oxide (N2O) in moles
REAL :: atm_n2o_mol(row_length,rows,model_levels)

! Atmospheric Burden of CFC-12 in moles
REAL :: atm_cf2cl2_mol(row_length,rows,model_levels)

! Atmospheric Burden of CFC-11 in moles
REAL :: atm_cfcl3_mol(row_length,rows,model_levels)

! Atmospheric Burden of CH3Br in moles
REAL :: atm_mebr_mol(row_length,rows,model_levels)

! Atmospheric Burden of H2 in moles 
REAL :: atm_h2_mol(row_length,rows,model_levels)

! Local dynamic arrays
REAL, ALLOCATABLE :: STASHwork34(:)
REAL, ALLOCATABLE :: STASHwork38(:)
REAL, ALLOCATABLE :: STASHwork50(:)

! Derived data arrays, allocatable if necessary
REAL :: land_fraction(row_length,rows)

! Temperature on theta levels
REAL, ALLOCATABLE    :: t_theta_levels(:,:,:)

! thickness of BL layers in metres
REAL, ALLOCATABLE :: thick_bl_levels(:,:,:)

REAL, ALLOCATABLE :: p_layer_boundaries(:,:,:)
REAL, ALLOCATABLE :: dj(:,:,:,:)   ! photolysis rates
REAL, ALLOCATABLE :: so4_sa(:,:,:) ! aerosol surface area
REAL, ALLOCATABLE :: t_chem(:,:,:) ! Temperature for chemistry
REAL, ALLOCATABLE :: q_chem(:,:,:) ! Specific humidity for chemistry
! working arrays for heterogeneous rate calculations
REAL :: trop_het_n2o5(row_length,rows,model_levels)
REAL :: trop_het_ho2 (row_length,rows,model_levels)


REAL :: delta_latitude
REAL :: delta_longitude

REAL :: z_top_of_model  ! top of model

! weighted cdnc = total cdnc * cldflg [m-3]
REAL :: cdncflag(row_length,rows,model_levels)
! weighted <cdnc^-1/3> = <cdnc^-1/3> * cldflg [m-3]
REAL :: cdnc3flag(row_length,rows,model_levels)
! weighted cdnc = total cdnc * liq_cloud_frac [m-3]
REAL :: cdncwt(row_length,rows,model_levels)

! Variables for COS initialisation
REAL :: y
REAL :: ymin
REAL :: ymax
REAL :: pmin
REAL :: pmax
REAL :: lp
REAL :: lpmin
REAL :: lpmax

REAL, SAVE :: lambda_aitken, lambda_accum ! parameters for computation
                                          ! of surface area density
LOGICAL :: do_chemistry
LOGICAL, SAVE :: firstchem = .TRUE.       ! Logical for any operations
                   ! specific to first chemical timestep, which is usually
                   ! different to the first ukca_main call. This should
                   ! only be used within 'do_chemistry' sections

! required for ASAD Flux Diagnostics
INTEGER :: ierr, stashsize
REAL, ALLOCATABLE :: fluxdiag_all_tracers(:,:,:,:)

! Min value to assign SO4 used in FastJX if no CLASSIC
REAL, PARAMETER :: min_SO4_val = 1.0e-25

REAL, PARAMETER :: max_z_for_offline_chem = 20000.0
! For backward-Euler Offline oxidants: only integrate chemistry below
!  this height. Note that in order to retain processor bit comparability
!  this must be set above the level corresponding to the
!  first_constant_r_rho_level

INTEGER :: first_constant_r_rho_level

LOGICAL, SAVE :: first       = .TRUE.   ! true only on 1st call

LOGICAL :: l_planet_obs = .FALSE. !  Calculate angles towards distant observer
LOGICAL :: l_store_value           ! T to store value, F to make difference

LOGICAL, PARAMETER :: unlump_species = .TRUE.
LOGICAL, PARAMETER :: lump_species = .FALSE.

! Mask to limit formation of Nat below specified height
LOGICAL, SAVE, ALLOCATABLE  :: have_nat3d(:,:,:)
LOGICAL                     :: level_above_ht
REAL,PARAMETER :: nat_limit_ht = 1000.  ! 1km
REAL           :: h_atmos   ! depth of atmosphere

CHARACTER (LEN=11) :: emiss_input  ! 'Ancillary  ' or  'NetCDF file'

! ErrorStatus
INTEGER                    :: errcode=0     ! Error flag (0 = OK)
CHARACTER(LEN=errormessagelength)   :: cmessage      ! Error return message

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='UKCA_MAIN1'

!- End of header
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 1. Initial set up
! ----------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

z_top_of_model = A_realhd(rh_z_top_theta)

m_atm_modl = submodel_for_sm(atmos_im)   ! submodel code
timestep = secs_per_stepim(atmos_im) ! timestep in seconds
timestep_number = stepim(atmos_im)   ! no. of steps since basis
A_Steps_per_hr = 3600*STEPS_PER_PERIODim(atmos_im)/                  &
                      SECS_PER_PERIODim(atmos_im)
radsteps_per_day=A_Steps_per_hr*24/a_sw_radstep_prog

! Determine number of model steps in a day
days_in_period = secs_per_periodim(atmos_im)/isec_per_day 
asteps_per_day = steps_per_periodim(atmos_im)/days_in_period

r_minute  = REAL(i_minute)     ! Real copy of i_minute
! Set Verbosity level for GLOMAP routines
glob_verbose = 0
IF ( PrintStatus > PrStatus_Oper ) glob_verbose = 2

! Set mpp arrays from parvars.h information

DO i= 0,nproc-1
  g_row_length(i) = g_lasize(1,fld_type_p,halo_type_no_halo,i)
  g_rows      (i) = g_lasize(2,fld_type_p,halo_type_no_halo,i)
END DO ! i processors

! Allocate diagnostic space for STASH
ALLOCATE (STASHwork34(STASH_maxlen(34,atmos_im)))
ALLOCATE (STASHwork38(STASH_maxlen(38,atmos_im)))
ALLOCATE (STASHwork50(STASH_maxlen(50,atmos_im)))
STASHwork34=0.0
STASHwork38=0.0
STASHwork50=0.0

! ----------------------------------------------------------------------
! Initialise ASAD constants - take UM values if requested
!  (values calculated in atmos_physics1 called from atm_step)
! code initially from asad_cinit (called from ukca_iniasad),
! but moved here since may need to be updated every timestep
! ----------------------------------------------------------------------
IF (L_ukca_set_trace_gases) THEN
  IF (l_co2_interactive) THEN
    ! Set to rmdi as fco2 should not be used
    fco2 = rmdi
  ELSE
    fco2 = um_co2_for_ukca/c_co2
  END IF
  fh2  = um_h2_for_ukca/c_h2
  fn2  = um_n2_for_ukca/c_n2
  fo2  = um_o2_for_ukca/c_o2
  fch4 = um_ch4_for_ukca/c_ch4
ELSE ! present day values
  fco2 = 350.0e-6
  fh2  = 5.0e-7
  fn2  = 0.78084
  fo2  = 0.20945
  fch4 = 1.76e-6
END IF

IF (first) THEN
  IF (PrintStatus >= PrStatus_Oper) THEN
    WRITE(umMessage,'(A)') 'First call to ukca_main'
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,E14.2)') 'timestep = ',timestep
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    WRITE(umMessage,'(A,I8)') 'timestep number = ',timestep_number
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END IF

  ! Initialise chemical definition arrays - needs to be here in case we are
  ! using ASAD
  CALL ukca_chem1_init()


  IF (L_ukca_chem) THEN

     IF (L_ukca_asad_columns) THEN
        CALL ukca_iniasad(model_levels)
     ELSE
        CALL ukca_iniasad(theta_field_size)
     END IF

    ! Initialise chemical diagnostics for N-R chemistry; set L_asad_use_chem_diags
    IF (.NOT. (L_ukca_trop .OR. L_ukca_aerchem .OR. L_ukca_raq .OR.   &
        l_ukca_raqaero)) THEN
      CALL asad_load_default_fluxes()
      ! check to see if chemical diagnostics are required, and
      ! set up stash for them.
      CALL asad_setstash_chemdiag(row_length,rows,model_levels)
      ! initialise chemical diagnostics
      IF (L_asad_use_chem_diags) &
           CALL asad_init_chemdiag(icode)
    END IF

! END IF  ! L_ukca_chem

  ! Check consistency of dry deposition scheme and number of tiles
  IF ( (ntiles /= ntype) .AND. l_ukca_intdd ) THEN
    WRITE(umMessage,'(A,A,I4,A,I4)') 'Cannot use interactive dry dep.',&
    ' You have ',ntiles,' surface tiles but should have ',ntype
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    errcode = 1
    cmessage='UKCA: Cannot use interactive dry dep'

    CALL ereport('UKCA_MAIN',errcode,cmessage)
  END IF

  ! Check that some CLASSIC aerosol types/processes are modelled if
  ! heterogeneous reactions on CLASSIC aerosols are used in UKCA.
  l_classic_aerosols = l_sulpc_so2            .OR. l_soot              .OR. &
                       l_ocff                 .OR. l_use_biogenic      .OR. &
                       l_use_seasalt_autoconv .OR. l_use_seasalt_sulpc .OR. &
                       l_use_seasalt_indirect .OR. l_use_seasalt_direct
  !
  IF (l_ukca_classic_hetchem .AND. .NOT. l_classic_aerosols) THEN
    cmessage = 'Cannot use heterogeneous reactions on CLASSIC aerosols ' // &
               'if these are not being modelled'
    errcode  = 5
    CALL ereport ('UKCA_MAIN', errcode, cmessage)
  END IF

  ! Heterogeneous reactions on CLASSIC aerosols are only allowed for the
  ! RAQ chemistry scheme, although this can be exteded in the future.
  IF (l_ukca_classic_hetchem .AND. .NOT. l_ukca_raq) THEN
    cmessage = 'Heterogeneous chemistry on CLASSIC aerosols is a '       // &
               'feature only available for the RAQ chemistry scheme'
    errcode  = 6
    CALL ereport ('UKCA_MAIN', errcode, cmessage)
  END IF

  END IF  ! L_ukca_chem
  ! Set up structure holding non-transported prognostics
  ! Needs to be called after ukca_chem1_init and before
  ! ukca_setd1defs
  CALL ntp_init(all_ntp)

  ! Define D1 definitions, emissions and diagnostics required
  CALL ukca_setd1defs(row_length, rows, n_rows, bl_levels,              &
           tr_levels, land_points, sm_levels, ntiles, tr_ukca, all_ntp)

  ! UKCA_CALC_CSPECIES needs to be called after SETD1DEFS
  IF (L_ukca_chem) THEN
    ALLOCATE(c_species(n_chem_tracers+n_aero_tracers))
    ALLOCATE(c_na_species(jpspec-jpctr))
    CALL ukca_calc_cspecies()

    errcode = 0

    DO i=1,jpctr        ! to allow for case when jpctr < n_chem_tracers...
      IF (c_species(i) < 0.001) THEN
        ! check that ratio of molec. wts. array has been initialised.
        errcode = errcode +1
        icode = -1*errcode
        cmessage=' c_species array is zero for '//advt(i)
        CALL ereport('UKCA_MAIN',icode,cmessage)
      END IF
    END DO

    IF (errcode > 0) THEN
      cmessage=' c_species array has zero values'//               &
                   ', check UKCA_CSPECIES routine'
      CALL ereport('UKCA_MAIN',errcode,cmessage)
    END IF

    DO i=1,nnaf
      IF (c_na_species(i) < 0.001) THEN
        ! check that ratio of molec. wts. array has been initialised.
        cmessage=' c_na_species array contains zero value'
        errcode = i
        CALL ereport('UKCA_MAIN',errcode,cmessage)
      END IF
    END DO
    ! L_ukca_advh2o *MUST* be .TRUE. if using H2O as 'TR'
    IF ((n_h2o > 0) .AND. (.NOT. L_ukca_advh2o)) THEN
      cmessage=                                                  &
           ' H2O is defined as a tracer, but L_UKCA_ADVH2O=F'
      errcode = n_h2o
      CALL ereport('UKCA_MAIN',errcode,cmessage)
    END IF
  END IF  ! L_ukca_chem

! Ensure that particulate fraction of SO2 emissions is between 0 and 5%.
  IF ( ( ANY(em_chem_spec == 'SO2_high  ') .OR.                  &
         ANY(em_chem_spec == 'SO2_nat   ') .OR.                  &
         ANY(em_chem_spec == 'SO2_low   ') ) .AND. (             &
         mode_parfrac < 0.0 .OR. mode_parfrac > 5.0) ) THEN
    WRITE(cmessage,'(A,F16.8)') 'SO2 emissions are on but ' //   &
          'mode_parfrac is invalid (should be 0.0-5.0): ', mode_parfrac
    errcode = 1
    CALL ereport('UKCA_MAIN',errcode, cmessage)
  END IF

  ! Setup a logical array for permitting NitricAcidTrihydrate
  ! formation, depending on height of each grid box
  IF ( .NOT. ALLOCATED(have_nat3d))                             &
       ALLOCATE(have_nat3d(row_length,rows,model_levels))
  DO k = 1, model_levels
   DO j = 1, rows
    DO i = 1, row_length
      ! True (allow NAT everywhere) if l_ukca_limit_nat=False,
      !      which is the default value for logical.
      ! False (supress everywhere) initially if logical is True
      have_nat3d(i,j,k) = (.NOT. l_ukca_limit_nat)
    END DO
   END DO
  END DO
  IF ( l_ukca_limit_nat ) THEN
    h_atmos = r_theta_levels(1,1,model_levels) - planet_radius
    DO k = 1, model_levels
     level_above_ht = ( (h_atmos * eta_theta_levels(k)) > nat_limit_ht ) 
     DO j = 1, rows
      DO i = 1, row_length
        have_nat3d(i,j,k) = level_above_ht
      END DO
     END DO
    END DO
  END IF

END IF ! first

! Set logical for CMIP6 MODE diagnostics if any item in range is requested
! Note that this is done at every timestep, as these diagnostics are only
!  required when the diagnostics is requested
item_first = stashcode_so4_nuc_sol - 1000*stashcode_glomap_sec
item_last = stashcode_h2o_mmr - 1000*stashcode_glomap_sec
IF (l_ukca_mode .AND. ANY(sf(item_first:item_last,stashcode_glomap_sec)))   &
  l_ukca_cmip6_diags = .TRUE.

! Set logical for PM10 and PM2.5 MODE diagnostics if any item in range is
! requested
item_first = stashcode_pm_first - 1000*stashcode_glomap_sec
item_last = stashcode_pm_last - 1000*stashcode_glomap_sec
IF (l_ukca_mode .AND. ANY(sf(item_first:item_last,stashcode_glomap_sec)))   &
  l_ukca_pm_diags = .TRUE.

! ----------------------------------------------------------------------
! Loop through all the objects in D1 to find address and characteristics
! of the requested items
! ----------------------------------------------------------------------
CALL ukca_um_d1_initialise(&
m_atm_modl, tdims%k_start, first)

IF (first) THEN

  IF ( l_ukca_chem )  THEN
  ! set constants needed in sulphur chemistry
  IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
     (L_ukca_stratcfc)) THEN
    lambda_aitken = 3.0/rhosul / chi /                              &
                   rad_ait * EXP(-2.5 * (LOG(sigma)**2))            &
                   * 0.01 ! revert from m^-1 to cm^2/cm^3
    lambda_accum  = 3.0/rhosul / chi /                              &
                   rad_acc * EXP(-2.5 * (LOG(sigma)**2))            &
                   * 0.01 ! revert from m^-2 to cm^2/cm^3
  END IF

! Check that solver interval has been set correctly in ukca_init
! This is not relevant if running without chem i.e. Age-of-air mode

  icode = 0
  IF (ukca_int_method == int_method_nr .OR.                        &
      ukca_int_method == int_method_be_explicit) THEN
    ! interval can be greater than 1 for BE explicit and N-R solvers
    IF (interval == imdi) icode = 101
  ELSE
    ! only allowed to use interval = 1 for other solvers
    IF (interval == imdi .OR. interval /= 1) icode = 102
  END IF
  IF (icode > 0) THEN
    cmessage = 'Error in UKCA_MAIN for solver interval'
    WRITE(umMessage,'(A40,A12,I6)') cmessage,' Interval: ',interval
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    CALL ereport('UKCA_MAIN', icode, cmessage)
  END IF
  END IF      ! l_ukca_chem

END IF  ! first

! allocate ASAD chemical diagnostics and load the requested diagnostics
IF (L_asad_use_chem_diags) THEN
  CALL asad_allocate_chemdiag(row_length,rows)
END IF


! decide whether to do chemistry
IF ( L_ukca_chem ) THEN
  do_chemistry = (MOD(timestep_number, interval) == 0)
ELSE
  do_chemistry = .FALSE.
END IF

! call subroutine which populates the all_tracers array using
! data from tracer_ukca_um and q_um
CALL ukca_all_tracers_copy_in(tracer_ukca_um, q_um)

! set up nmr and mmr indices on first call, only if using 
! GLOMAP-mode
IF (l_ukca_mode .AND. first) CALL ukca_aero_tracer_init()

! ----------------------------------------------------------------------
! 3. Calculate any derived variables
! ----------------------------------------------------------------------

! Derived variables for chemistry
IF (L_ukca_chem) THEN
  i1 = LBOUND(theta,dim=1)
  i2 = UBOUND(theta,dim=1)
  j1 = LBOUND(theta,dim=2)
  j2 = UBOUND(theta,dim=2)
  ALLOCATE(t_theta_levels(i1:i2,j1:j2,model_levels))
  ALLOCATE(rel_humid_frac(row_length,rows,model_levels))
  ALLOCATE(qsvp(1:row_length,1:rows,model_levels))
  ALLOCATE(qsmr(i1:i2,j1:j2,model_levels))
  ALLOCATE(Thick_bl_levels(1:row_length,1:rows,bl_levels))
  ALLOCATE(tile_frac(land_points,ntiles))              ! see below
  ALLOCATE(p_layer_boundaries(row_length,rows,0:model_levels))
  ALLOCATE(surf_albedo(row_length,rows))
  ALLOCATE(totnodens(row_length,rows,model_levels))
  ALLOCATE(ls_ppn3d(row_length,rows,model_levels))
  ALLOCATE(conv_ppn3d(row_length,rows,model_levels))
  ALLOCATE(so4_sa(row_length, rows, model_levels)) ! sulphate area density
! SO2 fluxes are required in chemistry_ctl
  ALLOCATE(delSO2_wet_h2o2(row_length,rows,model_levels))
  ALLOCATE(delSO2_wet_o3(row_length,rows,model_levels))
  ALLOCATE(delh2so4_chem(row_length,rows,model_levels))
  ALLOCATE(delSO2_drydep(row_length,rows,model_levels))
  ALLOCATE(delSO2_wetdep(row_length,rows,model_levels))
  ALLOCATE(qsatmr(i1:i2,j1:j2,model_levels))
  ALLOCATE(qsatmr_wat(i1:i2,j1:j2,model_levels))
  ALLOCATE(rel_humid_frac_clr(1:row_length,1:rows,model_levels))
  ALLOCATE(q1d(row_length*rows))
  ALLOCATE(q_clr1d(row_length*rows))
  ALLOCATE(qsatmr1d(row_length*rows))
  ALLOCATE(qsatmr_wat1d(row_length*rows))
  ALLOCATE(qcf1d(row_length*rows))
  ALLOCATE(cloud_liq_frac1d(row_length*rows))
  ALLOCATE(cloud_frac1d(row_length*rows))
  ALLOCATE(rhcrit1d(row_length*rows))
  ALLOCATE(rel_humid_frac_clr1d(row_length*rows))

  ! Initialise fluxes to zero as may not be filled by chemistry
  delSO2_wet_h2o2(:,:,:)=0.0
  delSO2_wet_o3(:,:,:)=0.0
  delh2so4_chem(:,:,:)=0.0
  delSO2_drydep(:,:,:)=0.0
  delSO2_wetdep(:,:,:)=0.0
! END IF ukca_chem

! Derived variables for stratopheric flux diagnostics
IF (.NOT. ALLOCATED(strat_fluxdiags)) THEN
  IF (.NOT. L_ukca_stratflux) n_strat_fluxdiags=0
  IF (n_strat_fluxdiags > 0) THEN
    ALLOCATE(strat_fluxdiags(row_length,rows,model_levels,          &
         n_strat_fluxdiags))
    strat_fluxdiags=0.0
  END IF
END IF

! Required in call to UKCA_MODE, but may be unallocated
IF (.NOT. ALLOCATED(mode_diags)) THEN
  IF (.NOT. L_ukca_mode) n_mode_diags=1
  ALLOCATE(mode_diags(row_length,rows,model_levels,               &
           n_mode_diags))
  mode_diags=0.0
END IF

! Required in call to UKCA_EMISS_CTL, but may be unallocated
IF (.NOT. ALLOCATED(dms_sea_conc)) THEN
  ALLOCATE(dms_sea_conc(row_length,rows))
  dms_sea_conc(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(chloro_sea)) THEN
  ALLOCATE(chloro_sea(row_length,rows))
  chloro_sea(:,:) = 0.0
END IF
IF (.NOT. ALLOCATED(so4_aitken)) THEN
  ALLOCATE(so4_aitken(row_length,rows,model_levels))
  so4_aitken(:,:,:) = min_SO4_val
END IF
IF (.NOT. ALLOCATED(so4_accum)) THEN
  ALLOCATE(so4_accum(row_length,rows,model_levels))
  so4_accum(:,:,:) = min_SO4_val
END IF
IF ( .NOT. ALLOCATED(dust_flux)) THEN
! Not allocated if L_ukca_dust=False, but passed to ukca_emiss_ctl
  ALLOCATE(dust_flux(row_length,rows,ndiv))
  dust_flux(:,:,:) = 0.0
END IF

! Required in call to UKCA_AERO_CTL, but may be unallocated
IF (.NOT. ALLOCATED(rim_agg)) THEN
  ALLOCATE(rim_agg(row_length,rows,model_levels))
  rim_agg=0.0
END IF
IF (.NOT. ALLOCATED(rim_cry)) THEN
  ALLOCATE(rim_cry(row_length,rows,model_levels))
  rim_cry=0.0
END IF

! initialise budget diagnostics at top of routine
nat_psc = 0.0
trop_ch4_mol = 0.0
trop_o3_mol = 0.0
trop_oh_mol = 0.0
strat_ch4_mol = 0.0
strat_ch4loss = 0.0
atm_ch4_mol = 0.0
atm_co_mol = 0.0
atm_n2o_mol = 0.0
atm_cf2cl2_mol = 0.0
atm_cfcl3_mol = 0.0
atm_mebr_mol = 0.0
atm_h2_mol = 0.0

  IF ( first ) THEN

    ALLOCATE(theta_latitude(row_length, rows))
    ALLOCATE(volume(row_length,rows,model_levels))
    ALLOCATE(delta_r(row_length,rows,model_levels))
    ALLOCATE(interf_z(row_length,rows,0:model_levels))
    ALLOCATE(area(row_length,rows,model_levels))
    ALLOCATE(mass(row_length, rows, model_levels))
    ALLOCATE(p_tropopause(row_length,rows))
    ALLOCATE(tropopause_level(row_length,rows))
    ALLOCATE(theta_trop(row_length,rows))
    ALLOCATE(pv_trop(row_length,rows))
    ALLOCATE(L_troposphere(row_length,rows,model_levels))
    ALLOCATE(solid_angle(row_length, rows))
    IF (.NOT. ALLOCATED(f3u_by_omega))                                &
       ALLOCATE(f3u_by_omega(row_length,rows))


    ! For V-At-Poles, V latitude already bounds the grid-cells on both sides

    ALLOCATE(v_latitude(row_length, rows+1))
    ALLOCATE(sinv_latitude(row_length, rows+1))
    DO j=1,rows+1
      v_latitude(:,j) = base_phi + (datastart(2)+j-2) * delta_phi
    END DO
    sinv_latitude = SIN(v_latitude)
    solid_angle(:,1:rows) = delta_lambda * (sinv_latitude(:,2:rows+1) -  &
                                              sinv_latitude(:,1:rows))

    DO k=1,model_levels
      area(1:row_length,1:rows,k) = 2*pi*                     &
        r_theta_levels(1:row_length,1:rows,k)**2 *                       &
        (sinv_latitude(1:row_length,2:rows+1) -                          &
         sinv_latitude(1:row_length,1:rows))/                            &
         global_row_length
    END DO

   ! Calculate f3_at_u/two_omega once, since used multiple times.
   ! For EG f3_at_u(0:row_length-1,1:rows), hence i+1,j
    DO j= udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        f3u_by_omega(i+1,j) = f3_at_u(i,j)/two_omega
      END DO
    END DO

   ! Compute the grid-box volume for the theta grid.
   ! This is not valid for use with quantities held on rho grid 
   ! e.g. density and relying on global conservation
    DO k=2,model_levels-1
      delta_r(:,:,k) = (r_rho_levels(1:row_length,1:rows,k+1) -     &
        r_rho_levels(1:row_length,1:rows,k))
    END DO
    delta_r(:,:,1) = (r_rho_levels(1:row_length,1:rows,2) -         &
      r_theta_levels(1:row_length,1:rows,0))
    delta_r(:,:,model_levels) =                                     &
      (r_theta_levels(1:row_length,1:rows,model_levels) -           &
       r_rho_levels(1:row_length,1:rows,model_levels))

    DO k=1,model_levels
      volume(:,:,k) = r_theta_levels(1:row_length,1:rows,k) *       &
        r_theta_levels(1:row_length,1:rows,k) *                     &
        delta_r(:,:,k) * delta_phi * delta_lambda *                 &
        FV_cos_theta_latitude(1:row_length,1:rows)
    END DO

    ! Calculate altitude above ground for the interfaces (in metres).
    ! Needed by some vertical profiles in the new emission system.
    interf_z         (1:row_length, 1:rows, 0) = 0.0

    DO k = 1, model_levels-1
      interf_z          (1:row_length, 1:rows, k) =                  &
         r_rho_levels   (1:row_length, 1:rows, k+1) -                &
         r_theta_levels (1:row_length, 1:rows, 0)

    END DO

    interf_z         (1:row_length, 1:rows, model_levels) =         &
      r_theta_levels (1:row_length, 1:rows, model_levels) -         &
      r_theta_levels (1:row_length, 1:rows, 0)

    IF (ALLOCATED(theta_latitude)) DEALLOCATE(theta_latitude)
    IF (ALLOCATED(v_latitude)) DEALLOCATE(v_latitude)
    IF (ALLOCATED(sinv_latitude)) DEALLOCATE(sinv_latitude)

    ALLOCATE(z_half(row_length,rows,bl_levels))
    ALLOCATE(z_half_alllevs(row_length,rows,model_levels))

    DO k = 1,model_levels
      DO j = 1, rows
        DO i= 1, row_length
          z_half_alllevs(i,j,k)=                                    &
                       r_rho_levels(i,j,k)-r_theta_levels(i,j,0)
        END DO
      END DO
    END DO
    DO k = 1,bl_levels
      z_half(:,:,k) = z_half_alllevs(:,:,k)
    END DO

! Calculate last model level for integration of offline oxidants (BE)
!  and check that it is higher than the level of first_constant_r_rho_level
    k_be_top = model_levels        ! default
    IF (l_ukca_offline_be) THEN
      DO k = 2,model_levels
        IF ((MINVAL(r_theta_levels(:,:,k)) - planet_radius) >            &
             max_z_for_offline_chem) THEN
          k_be_top = k
          EXIT
        END IF
      END DO

      IF (PrintStatus >=  PrStatus_Oper) THEN
        WRITE(umMessage,'(A50,I5)') 'Top level for offline oxidants'//  &
                                    ' integration: ',k_be_top
        CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      END IF

      first_constant_r_rho_level = a_inthd(24)
      IF (max_z_for_offline_chem <                                      &
          (MAXVAL(r_rho_levels(:,:,first_constant_r_rho_level))         &
                                               - planet_radius)) THEN
        cmessage = ' max_z_for_offline_chem is set too low for processor'// &
                   ' bit comparability'
        errcode = 1
        WRITE(umMessage,'(A74,2E12.4)') cmessage,max_z_for_offline_chem, &
             MINVAL(r_rho_levels(:,:,first_constant_r_rho_level))-planet_radius
        CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        CALL ereport('UKCA_MAIN1',errcode,cmessage)
      END IF
    END IF    ! l_ukca_offline_be

  END IF     ! first

  ! Read and interpolate aerosol surface area density.
  ! Depending on the choice of logicals, the array can contain
  ! surface area densities based on:
  ! - Climatological aerosol in whole domain
  ! - Climatology (above 12km)+ CLASSIC/AeroClim aerosol (below 12km)
  ! - Climatology (above 12km)+ GLOMAP aerosol  (below 12km)
  ! - GLOMAP aerosol in the whole domain
  so4_sa(:,:,:) = 0.0
  cmessage = ' '
  IF (do_chemistry .AND. L_ukca_het_psc) THEN
    IF (L_ukca_sa_clim) THEN
      cmessage =                                                  &
        'UKCA_Het_PSC: Using only climatological aerosol surface area'
      CALL ukca_read_aerosol(i_year, i_month, row_length, rows,   &
           model_levels,                                          &
           sin_theta_latitude(1:row_length, 1:rows),              &
           eta_theta_levels(1:model_levels) * z_top_of_model,     &
           L_ukca_use_background_aerosol, so4_sa)
      ! Below uses SO4 tracer from classic scheme for troposphere only
      IF ((l_sulpc_so2 .OR. l_use_arclsulp) .AND.                 &
          (ALLOCATED(so4_aitken) .AND. ALLOCATED(so4_accum))) THEN
        cmessage =                                                &
         'UKCA_Het_PSC: Using Climatology + CLASSIC/climatology'//&
         ' aerosol surface area'
        ! Calculate surface area density fields from aerosol tracers
        !   and copy into cells below 12 km
        DO k=1,model_levels
          IF (eta_theta_levels(k) * z_top_of_model < 12000.0) THEN
            so4_sa(:,:,k) = (                                     &
             lambda_aitken * so4_aitken(1:row_length,1:rows,k)    &
            +lambda_accum  * so4_accum (1:row_length,1:rows,k))   &
            *rho_r2(1:row_length,1:rows,k)/                       &
             (r_theta_levels(1:row_length,1:rows,k)**2)
          END IF
        END DO
      ELSE IF ( L_ukca_mode )  THEN
        ! Use aerosol surface area from GLOMAP for cells below 12km
        cmessage =                                                &
         'UKCA_Het_PSC: Using Climatology + GLOMAP aerosol surface area'
        i = name2ntpindex(all_ntp,'surfarea  ')
        DO k=1,model_levels
          IF (eta_theta_levels(k) * z_top_of_model < 12000.0)     &
              so4_sa(1:row_length,1:rows,k) =                     &
                   all_ntp(i)%data_3d(1:row_length,1:rows,k)
        END DO
      END IF    ! l_sulpc_so2/ Mode

    ELSE        ! Not ukca_sa_clim, take aerosol fully from GLOMAP

      IF (L_ukca_mode ) THEN
        cmessage='UKCA_Het_PSC: Using only GLOMAP aerosol surface area '
        i = name2ntpindex(all_ntp,'surfarea  ')
        so4_sa(1:row_length,1:rows,1:model_levels) =               &
             all_ntp(i)%data_3d(1:row_length,1:rows,1:model_levels)
      ELSE
        cmessage='UKCA HET_PSC: Surface area undefined. At least '// &
          'one from l_ukca_sa_clim and l_ukca_glomap has to be selected.'
        errcode = 1
        CALL ereport('UKCA_MAIN1',errcode,cmessage)
      END IF     ! l_ukca_mode
    END IF       ! L_ukca_sa_clim

    WHERE(so4_sa(:,:,:) < 0.0)    ! Remove any negative values
      so4_sa(:,:,:) = 0.0
    END WHERE
    ! Inform users of the choice active for this run
    IF ( firstchem .AND. PrintStatus > PrStatus_Min )             &
         CALL umPrint(cmessage,src='ukca_main1-ukca_main1')

  END IF         ! do_chemistry and l_ukca_het_psc

  t_theta_levels(:,:,1:model_levels)=                              &
       exner_theta_levels(:,:,1:model_levels) *                    &
       theta(:,:,1:model_levels)
  Thick_bl_levels(1:row_length,1:rows,1) = 2.0*                    &
        (r_theta_levels(1:row_length,1:rows,1) -                   &
         r_theta_levels(1:row_length,1:rows,0))
  DO k=2,bl_levels
    Thick_bl_levels(1:row_length,1:rows,k) =                       &
          r_rho_levels(1:row_length,1:rows,k+1) -                  &
          r_rho_levels(1:row_length,1:rows,k)
  END DO
  p_layer_boundaries(1:row_length,1:rows,0)                =       &
          pstar(1:row_length,1:rows)
  p_layer_boundaries(1:row_length,1:rows,model_levels)     =       &
          p_theta_levels(1:row_length,1:rows,model_levels)
  p_layer_boundaries(1:row_length,1:rows,1:model_levels-1) =       &
          p_rho_levels(1:row_length,1:rows,2:model_levels)


  ! Mass:
  ! We make the hydrostatic assumption in the diagnostic mass calculation
  ! here so that mass = -b*(solid_angle/(3*g))*(r_top^3 - r_bottom^3)
  ! where   b =                                                           &
  !    &    (p_layer_boundaries(:,:,k) - p_layer_boundaries(:,:,k-1))/    &
  !    &    (r_theta_levels(:,:,k)     - r_theta_levels(:,:,k-1))         &


  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        mass(i,j,k) = (-solid_angle(i,j) / (3.0*g)) *               &
      ((p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1))/   &
       (r_theta_levels(i,j,k) - r_theta_levels(i,j,k-1)) ) *        &
       (r_theta_levels(i,j,k)**3.0 - r_theta_levels(i,j,k-1)**3.0)
      END DO
    END DO
  END DO


  DO k=1,model_levels
    DO j=1,rows
      DO i=1,row_length
        totnodens(i,j,k) = p_theta_levels(i,j,k)                   &
                    /(zboltz*t_theta_levels(i,j,k))  ! molec/m3
        conv_ppn3d(i,j,k) = conv_rain3d(i,j,k) + conv_snow3d(i,j,k)
        ls_ppn3d(i,j,k)   = ls_rain3d  (i,j,k) + ls_snow3d  (i,j,k)

        conv_ppn3d(i,j,k) = MAX(conv_ppn3d(i,j,k), 0.0)
        ls_ppn3d  (i,j,k) = MAX(ls_ppn3d  (i,j,k), 0.0)
      END DO
    END DO
  END DO

  i1=LBOUND(qsmr,1)  ! lower bound of 1st dimension of array
  j1=LBOUND(qsmr,2)  ! lower bound of 2nd dimension of array
  DO k = 1, model_levels

    IF (l_new_qsat_chem_aero) THEN
      CALL qsat_wat_mix_new(qsmr(:,:,k),t_theta_levels(:,:,k),        &
              p_theta_levels(:,:,k),                                  &
              SIZE(t_theta_levels(:,1,k)),SIZE(t_theta_levels(1,:,k)))
    ELSE
      ! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsmr(i1,j1,k),t_theta_levels(i1,j1,k),        &
              p_theta_levels(i1,j1,k),SIZE(t_theta_levels(:,:,k)),    &
              .TRUE.)
    END IF

    ! Derive saturated vapour pressure from saturated mixing ratio
    qsvp(:,:,k) = qsmr(1:row_length,1:rows,k)*                    &
                  p_theta_levels(1:row_length,1:rows,k)/          &
                  (repsilon + qsmr(1:row_length,1:rows,k)*        &
                  (1.0 - repsilon))

    ! Convert to relative humidity fraction
    rel_humid_frac(1:row_length,1:rows,k) =                       &
            q(1:row_length,1:rows,k)/                             &
            qsmr(1:row_length,1:rows,k)
  END DO
  IF (ALLOCATED(qsmr)) DEALLOCATE(qsmr)

  WHERE (rel_humid_frac < 0.0)
    rel_humid_frac=0.0           ! remove negatives
  END WHERE
  WHERE (rel_humid_frac >= 1.0)
    rel_humid_frac=0.999         ! remove values >= 1.0
  END WHERE

  !-----------------------------------------------------------------------
  ! Clear sky relative humidity calculation
  ! DEPENDS ON: qsat_wat_mix
  DO k = 1, model_levels

    !         qsatmr     is w.r.t  ice   when T < 0C
    !         qsatmr_wat is w.r.t  water when T < 0C

    IF (l_new_qsat_chem_aero) THEN
      CALL qsat_new(qsatmr(:,:,k),t_theta_levels(:,:,k),              &
              p_theta_levels(:,:,k),                                  &
              SIZE(t_theta_levels(:,1,k)),SIZE(t_theta_levels(1,:,k)))
      CALL qsat_wat_mix_new(qsatmr_wat(:,:,k),t_theta_levels(:,:,k),  &
              p_theta_levels(:,:,k),                                  &
              SIZE(t_theta_levels(:,1,k)),SIZE(t_theta_levels(1,:,k)))
    ELSE
      ! DEPENDS ON: qsat
      CALL qsat(qsatmr(i1,j1,k),t_theta_levels(i1,j1,k),              &
              p_theta_levels(i1,j1,k),SIZE(t_theta_levels(:,:,k)))
      ! DEPENDS ON: qsat_wat_mix
      CALL qsat_wat_mix(qsatmr_wat(i1,j1,k),t_theta_levels(i1,j1,k),  &
              p_theta_levels(i1,j1,k),SIZE(t_theta_levels(:,:,k)),    &
              .TRUE.)
    END IF

    ! Prepare 1D arrays for input to lsp_qclear
    q1d              = RESHAPE(q(1:row_length,1:rows,k),              &
         (/SIZE(q(1:row_length,1:rows,k))/))
    qsatmr1d           = RESHAPE(qsatmr(1:row_length,1:rows,k),       &
         (/SIZE(qsatmr(1:row_length,1:rows,k))/))
    qsatmr_wat1d       = RESHAPE(qsatmr_wat(1:row_length,1:rows,k),   &
         (/SIZE(qsatmr_wat(1:row_length,1:rows,k))/))
    qcf1d            = RESHAPE(qcf(1:row_length,1:rows,k),            &
         (/SIZE(qcf(1:row_length,1:rows,k))/))
    cloud_liq_frac1d = RESHAPE(cloud_liq_frac(1:row_length,1:rows,k), &
         (/SIZE(cloud_liq_frac(1:row_length,1:rows,k))/))
    cloud_frac1d     = RESHAPE(cloud_frac(1:row_length,1:rows,k),     &
         (/SIZE(cloud_frac(1:row_length,1:rows,k))/))
    rhcrit1d(:)      = rhcrit(k)

    CALL lsp_qclear(q1d, qsatmr1d, qsatmr_wat1d, &
         qcf1d, cloud_liq_frac1d, cloud_frac1d,    &
         rhcrit1d, q_clr1d, SIZE(q1d))

    rel_humid_frac_clr1d(:) = q_clr1d(:) / qsatmr_wat1d(:)

    rel_humid_frac_clr(1:row_length,1:rows,k) =                   &
         RESHAPE(rel_humid_frac_clr1d, (/row_length, rows/))

  END DO

  ! Cap RH_clr values: 0.0 <=  RH_CLR <= 0.999
  WHERE (rel_humid_frac_clr < 0.0)
    rel_humid_frac_clr=0.0           ! remove negatives
  END WHERE
  WHERE (rel_humid_frac_clr >= 1.0)
    rel_humid_frac_clr=0.999           ! RH >= 1.0
  END WHERE

  IF (ALLOCATED(qsatmr)) DEALLOCATE(qsatmr)
  IF (ALLOCATED(qsatmr_wat)) DEALLOCATE(qsatmr_wat)
  IF (ALLOCATED(q1D)) DEALLOCATE(q1D)
  IF (ALLOCATED(q_clr1D)) DEALLOCATE(q_clr1D)
  IF (ALLOCATED(qsatmr1D)) DEALLOCATE(qsatmr1D)
  IF (ALLOCATED(qsatmr_wat1D)) DEALLOCATE(qsatmr_wat1D)
  IF (ALLOCATED(qcf1D)) DEALLOCATE(qcf1D)
  IF (ALLOCATED(cloud_liq_frac1D)) DEALLOCATE(cloud_liq_frac1D)
  IF (ALLOCATED(cloud_frac1D)) DEALLOCATE(cloud_frac1D)
  IF (ALLOCATED(rhcrit1D)) DEALLOCATE(rhcrit1D)
  IF (ALLOCATED(rel_humid_frac_clr1D)) DEALLOCATE(rel_humid_frac_clr1D)

  !-----------------------------------------------------------------------

! If land fraction is not already allocated set to 1.
! This only happens if coastal tiling is off.
! When coastal tiling is off all land points have a land fraction of 1
IF (.NOT. ALLOCATED(fland) ) THEN
  ALLOCATE(fland(land_points))
  fland(1:land_points)=1.0
END IF
END IF   ! If ukca_chem

IF (L_ukca_chem .OR. L_ukca_mode) THEN
  ! Set up land fraction (required for MODE and DRY DEPN).
  ! Input data fland is a 1D packed array on land points only.
  ! We require a 2D data valid at all points.
  land_fraction(:,:)=0.0
  DO l = 1, land_points
    j = (land_index(l)-1)/row_length + 1
    i = land_index(l) - (j-1)*row_length
    land_fraction(i,j) = fland(l)
  END DO

! ----------------------------------------------------------------------
! Print out values for debugging
! ----------------------------------------------------------------------

  CALL ukca_pr_inputs(first)

! Print values derived locally here
  IF (first .AND. PrintStatus >= PrStatus_Oper ) THEN
    WRITE(umMessage,'(A,I8, 2E12.4)') 'Thick_bl_levels (level 1): ',mype,      &
      MAXVAL(Thick_bl_levels(:,:,1)),                                          &
      MINVAL(Thick_bl_levels(:,:,1))
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    DO k=1,model_levels
      WRITE(umMessage,'(A,I6)') 'LEVEL: ',k
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 't_theta_levels: ',mype,k,         &
        MAXVAL(t_theta_levels(1:row_length,1:rows,k)),                         &
        MINVAL(t_theta_levels(1:row_length,1:rows,k))
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    END DO
    DO k=1,model_levels
      WRITE(umMessage,'(A,I6)') 'WET LEVEL: ',k
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'rel_humid_frac: ',mype,k,         &
        MAXVAL(rel_humid_frac(1:row_length,1:rows,k)),                         &
        MINVAL(rel_humid_frac(1:row_length,1:rows,k))
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    END DO
    DO k=1,bl_levels
      WRITE(umMessage,'(A,I6)') 'BL LEVEL: ',k
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      WRITE(umMessage,'(A,I8, I6, 2E12.4)') 'z_half:     ',mype,k,             &
        MAXVAL(z_half(1:row_length,1:rows,k)),                                 &
        MINVAL(z_half(1:row_length,1:rows,k))
      CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
    END DO
    WRITE(umMessage,'(A16,I5,2E12.4)') 'land_fraction: ', mype,                &
      MAXVAL(land_fraction), MINVAL(land_fraction)
    IF (l_sulpc_so2 .AND. i_ukca_photol == i_ukca_fastjx) THEN
      DO k=1,model_levels
        WRITE(umMessage,*) 'so4_aitken: ', mype,k,                       &
                 MAXVAL(so4_aitken(1:row_length,1:rows,k)),              &
                 MINVAL(so4_aitken(1:row_length,1:rows,k))
        CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        WRITE(umMessage,*) 'so4_accum: ', mype,k,                        &
                MAXVAL(so4_accum(1:row_length,1:rows,k)),                &
                MINVAL(so4_accum(1:row_length,1:rows,k))
        CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      END DO
      IF (ALLOCATED(sulphate_od)) THEN
        DO k=1,sw_spectrum(1)%basic%n_band
          WRITE(umMessage,*) 'sulphate_od: ', mype,k,                    &
                  MAXVAL(sulphate_od(1:row_length,1:rows,k)),            &
                  MINVAL(sulphate_od(1:row_length,1:rows,k))
          CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        END DO
      END IF
    END IF        ! sulpc_so2 and Fast-JX
    WRITE(umMessage,*) 'kent:     ', mype,MAXVAL(kent),    MINVAL(kent)
    CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
  END IF

  !       Warning statements about possible mismatch between interactive
  !       methane emissions and emissions ancillary
  IF (first) THEN
    emiss_input = 'NetCDF file'

    IF (L_ukca_qch4inter .AND. (.NOT. L_ukca_prescribech4)) THEN
      cmessage = 'CH4 WETLANDS EMS ARE ON - '  // TRIM(emiss_input) // &
                 ' SHOULD NOT contain wetland ems'
      errcode=-1

      CALL ereport('UKCA_MAIN1',errcode,cmessage)
    ELSE IF (.NOT. L_ukca_prescribech4 .AND. ANY(advt == 'CH4       ')) THEN
      cmessage = 'CH4 WETLANDS EMS ARE OFF - ' // TRIM(emiss_input) // &
                 ' SHOULD contain wetland ems'
      errcode=-2
      CALL ereport('UKCA_MAIN1',errcode,cmessage)
    END IF
  END IF  ! End of IF (first)

END IF    ! End of IF (L_ukca_chem .OR. L_ukca_mode)

!     Set up tile info
IF (L_ukca_chem) THEN
  CALL tilepts(land_points,frac_types,tile_pts,tile_index)

  !     Set tile fractions to 1 if aggregate tiles are used (NTILES=1).
  !     Otherwise, set tile fractions to surface type fractions.

  tile_frac=0.0
  IF (ntiles == 1) THEN
    DO l=1,land_points
      tile_frac(l,1) = 1.0
    END DO
  ELSE
    DO n=1,ntiles
      DO j=1,tile_pts(n)
        l = tile_index(j,n)
        tile_frac(l,n) = frac_types(l,n)
      END DO
    END DO
  END IF

END IF    ! L_ukca_chem


! ----------------------------------------------------------------------
! 4. Calls to science subroutines
! ----------------------------------------------------------------------
! 4.1 Age-of-air
! ----------------------------------------------------------------------

IF ( L_ukca_ageair ) THEN
  CALL ukca_age_air(timestep, z_top_of_model, all_tracers(:,:,1:,:))
END IF    ! l_ukca_ageair

! ----------------------------------------------------------------------
! 4.2 Chemistry related pre-processing
! ----------------------------------------------------------------------

IF (L_ukca_chem) THEN
  IF (ltimer) CALL timer('UKCA CHEMISTRY MODEL',5)

  ! Transform tracers to ensure elemental conservation
  IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR.                   &
     (L_ukca_stratcfc)) THEN
    CALL ukca_transform_halogen(tr_ukca,rows,row_length,            &
                   model_levels, tdims_s%halo_i,tdims_s%halo_j,     &
                   all_tracers(:,:,1:model_levels,:),               &
                   tdims_s%halo_i,tdims_s%halo_j,                   &
                   q(:,:,1:model_levels), unlump_species,timestep_number)
  END IF


  !       Calculate tropopause pressure using a combined
  !       theta and PV surface

  CALL ukca_calc_tropopause(row_length, rows, model_levels,      &
       sin_theta_latitude(1:row_length,1:rows),                  &
       theta(1:row_length,1:rows,1:model_levels),                &
       pv_on_theta_mlevs(1:row_length,1:rows,1:model_levels),    &
       p_layer_boundaries(1:row_length,1:rows,0:model_levels),   &
       p_theta_levels(1:row_length,1:rows,1:model_levels),       &
       p_tropopause(1:row_length,1:rows),                        &
       tropopause_level(1:row_length,1:rows))

  IF (L_asad_use_chem_diags .AND. L_asad_use_trop_mask) &
       CALL asad_tropospheric_mask( &
       rows,row_length,model_levels,ierr)



  !       Calculate the difference in tracers per timestep in the
  !       troposphere (moles/s) due to transport using the
  !       trmol_post_atmstep array from this timestep and the
  !       trmol_post_chem array from the previous timestep.
  !       Do this only for those stratospheric flux diagnostics
  !       which are switched on.

  IF (first .AND. n_strat_fluxdiags > 0) THEN

    ALLOCATE(trmol_post_chem(row_length,rows,model_levels,       &
                               n_chem_tracers))  ! moles
    trmol_post_chem = 0.0

    DO l=1,n_strat_fluxdiags
      DO k=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            strat_fluxdiags(i,j,k,l) = 0.0
          END DO
        END DO
      END DO
    END DO

  ELSE IF ((.NOT. first) .AND. n_strat_fluxdiags > 0) THEN

    ALLOCATE(trmol_post_atmstep(row_length,rows,model_levels,     &
                          n_chem_tracers))  ! moles
    trmol_post_atmstep = 0.0

    DO l=1,n_chem_tracers
      DO k=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            trmol_post_atmstep(i,j,k,l) = all_tracers(i,j,k,l)    &
                *totnodens(i,j,k)*volume(i,j,k)/(c_species(l) *avc)   ! moles
          END DO
        END DO
      END DO
    END DO

    icnt = 0
    DO l=1,n_chem_tracers
      IF (UkcaD1codes(istrat_first+l-1)%item /= imdi) THEN
        icnt = icnt + 1

        DO k=1,model_levels
          DO j=1,rows
            DO i=1,row_length
              IF (L_troposphere(i,j,k)) THEN              ! troposphere
                strat_fluxdiags(i,j,k,icnt) =                     &
                 (trmol_post_atmstep(i,j,k,l)-                    &
                  trmol_post_chem(i,j,k,l))/timestep      ! moles/sec
              ELSE                                        ! stratosphere
                strat_fluxdiags(i,j,k,icnt) = 0.0
              END IF
            END DO
          END DO
        END DO
      END IF
    END DO     ! n_chem_tracers

    IF (ALLOCATED(trmol_post_atmstep)) DEALLOCATE(trmol_post_atmstep)  ! moles

  END IF      ! End of IF first and n_strat_fluxdiags statement

  IF (.NOT. first) THEN ! DO NOT CALL ON FIRST TIMESTEP
     ! ASAD Flux Diags STE
    l_store_value = .FALSE.
    IF (L_asad_use_chem_diags .AND. L_asad_use_STE) &
         CALL asad_tendency_STE( &
         row_length,rows,model_levels, &
         n_chem_tracers,all_tracers, &
         totnodens,volume,mass,timestep, &
         calculate_STE,l_store_value,ierr)
  END IF

  l_store_value = .TRUE.
  IF (do_chemistry .AND. (L_asad_use_chem_diags &
       .AND. L_asad_use_tendency)) &
       CALL asad_tendency_STE( &
       row_length,rows,model_levels, &
       n_chem_tracers,all_tracers, &
       totnodens,volume,mass,timestep, &
       calculate_tendency,l_store_value,ierr)

  ! force area cloud fraction to be zero or positive to prevent
  ! errors in Fastj-X
  WHERE(area_cloud_fraction < 0.)
      area_cloud_fraction(:,:,:) = 0.
  END WHERE

! ----------------------------------------------------------------------
! 4.3 Photolysis
! ----------------------------------------------------------------------

  ! Do pre-processing for photolysis scheme

  IF ( first ) THEN
    ALLOCATE(land_albedo_all(row_length,rows,Radsteps_per_day))
    land_albedo_all=min_surf_albedo
  END IF

  ! Update land_albedo_all
  !  irad_prev and irad_next are indices to the slices of the albedo array
  !  from 'previous' and the 'next' Radiation timestep 
  !  which we use for the interpolation in time
  !  irad_prev should be 1 for timesteps 1 to a_sw_radstep_prog
  !  irad_prev 2 should be 2 for timesteps a_sw_radstep_prog+1 etc

  irad_prev = 1+MOD((timestep_number-1)/a_sw_radstep_prog,Radsteps_per_day)

  ! irad_next is irad_prev+1 unless irad_prev = Radsteps_per_day when it is 1
  irad_next = 1+MOD(irad_prev,Radsteps_per_day)

  ! Update land_albedo_all
  IF (l_rad_step_prog) THEN
    IF ( .NOT. ALLOCATED(land_albedo) ) THEN
      ! If coastal tiling is off then we do not have the land_albedo 
      ! array allocated
      i1=LBOUND(tot_surf_sw,1)
      i2=UBOUND(tot_surf_sw,1)
      j1=LBOUND(tot_surf_sw,2)
      j2=UBOUND(tot_surf_sw,2)
      ALLOCATE(land_albedo(i1:i2,j1:j2))
    END IF
    WHERE (tot_surf_sw > 0.0 .AND. tot_surf_sw < 1.5e3)
      land_albedo=(tot_surf_sw-net_surf_sw)/tot_surf_sw
    ELSEWHERE
      land_albedo=min_surf_albedo
    END WHERE
    land_albedo_all(:,:,irad_prev)=land_albedo(:,:)
  END IF

  IF (ALLOCATED(net_surf_sw)) DEALLOCATE(net_surf_sw)
  IF (ALLOCATED(tot_surf_sw)) DEALLOCATE(tot_surf_sw)

  ! The following (Fast-JX calculation) only needs to be done when
  ! chemistry is called:

  IF (do_chemistry) THEN

  ! Interpolate surf_albedo to timestep value - deactivated for all schemes 
  ! except RAQ (for now) as this is not NRUN-CRUN bit-reproducible
  
    IF ( l_ukca_raq .OR. l_ukca_raqaero ) THEN
      ! fx is the interpolation weight for irad_prev
      ! fx=1 for timesteps 1, 1+a_sw_radstep_prog etc
      ! fx=1/12 for timesteps a_sw_radstep_prog-1, 2*a_sw_radstep_prog-1 etc

      fx=1.0-REAL(MOD((timestep_number-1),a_sw_radstep_prog))/          &
                                      REAL(a_sw_radstep_prog)
      surf_albedo(:,:)=fx*land_albedo_all(:,:,irad_prev)+               &
                   (1.0-fx)*land_albedo_all(:,:,irad_next)
    ELSE
      ! Do not interpolate, use value from last Radiation call
      surf_albedo(:,:)=land_albedo(:,:)
    END IF
   
    ALLOCATE(um_ozone3d(row_length,rows,ozone_levels))
    IF (SIZE(um_ozone,dim=1) == 1) THEN   ! Use zonal average
      DO i=1,row_length
        um_ozone3d(i,:,:)=um_ozone(1,:,1:ozone_levels)
      END DO
    ELSE                                  ! 3-D ozone specifed
      um_ozone3d(:,:,:)= um_ozone(:,:,1:ozone_levels)
    END IF

    !       Convert 2D or 3D ozone field from UM to a 1D field

    nd_o3=SIZE(um_ozone)
    ALLOCATE(um_ozone1D(nd_o3))
    um_ozone1D = RESHAPE(um_ozone,(/nd_o3/))


    ALLOCATE(dj(row_length,rows,model_levels,jppj))
    dj=0.0

    IF ((i_ukca_photol == i_ukca_fastjx)) THEN
      IF (ltimer) CALL timer('UKCA FASTJX MODEL  ',5)

      ! DEPENDS ON: ukca_fastjx
      CALL ukca_fastjx(                                          &
      dj(1:row_length,1:rows,1:model_levels,1:jppj),             &
      p_layer_boundaries(1:row_length,1:rows,0:model_levels),    &
      pstar(1:row_length,1:rows),                                &
      p_theta_levels(1:row_length,1:rows,1:model_levels),        &
      t_theta_levels(1:row_length,1:rows,1:model_levels),        &
      f3u_by_omega(1:row_length,1:rows),                         &
      rho_r2(1:row_length,1:rows,1:model_levels),                &
      z_top_of_model,                                            &
      so4_aitken(1:row_length,1:rows,1:model_levels),            &
      so4_accum(1:row_length,1:rows,1:model_levels),             &
      qcl(1:row_length,1:rows,1:model_levels),                     &
      qcf(1:row_length,1:rows,1:model_levels),                     &
      rel_humid_frac(1:row_length,1:rows,:),                     &
      area_cloud_fraction(1:row_length,1:rows,1:model_levels),     &
      conv_cloud_lwp(1:row_length,1:rows),                       &
      conv_cloud_top(1:row_length,1:rows),                       &
      conv_cloud_base(1:row_length,1:rows),                      &
      conv_cloud_amount(1:row_length,1:rows,1:model_levels),       &
      surf_albedo(1:row_length,1:rows),                          &
      tropopause_level(1:row_length, 1:rows),                    &
      all_tracers(1:row_length,1:rows,1:model_levels,n_o3),      &
      land_fraction)

      IF (ltimer) CALL timer('UKCA FASTJX MODEL  ',6)
    END IF  ! End of IF (i_ukca_fastjx) statement

    !       Clear workspace

    IF (ALLOCATED(um_ozone1D)) DEALLOCATE(um_ozone1D)
    IF (ALLOCATED(surf_albedo)) DEALLOCATE(surf_albedo)

  END IF   ! do_chemistry

! ----------------------------------------------------------------------
! 4.4 Emissions
! ----------------------------------------------------------------------
  !       Allocate array for wetland ch4 emissions when
  !       L_ukca_qch4inter is false. Set emissions to zero.
  IF (.NOT. (L_ukca_qch4inter)) THEN
    ALLOCATE(ch4_wetl_emiss(1:row_length,1:rows))
    ch4_wetl_emiss(1:row_length,1:rows) = 0.0
  END IF

  ! Emission part. Before calling emissions,
  ! set up array of lower-boundary mixing ratios for stratospheric species.
  ! Note that some mixing ratios (CH4, N2O, F11, F12, F113, H22) are defined
  ! in the RUN_UKCA namelist held in UKCA module ukca_option_mod through the
  ! gui, others are not. If for some species
  ! actual emissions (not prescribed lower boundary conditions) are required,
  ! remove these species from lbc_spec, e.g, for CH4.
  ! Mass mixing ratios are best estimates, taken from the WMO Scientific
  ! assessment of ozone depletion (2002), for 2000 conditions.
  ! H2 corresponds to 0.5 ppmv.
  ! Select lower boundary MMRs. Either from gui or predefined constants
  ! Note that for bromine compounds (MeBr, H1211, H1301) are increased
  ! by 24% to account for non-included species and HBr is set to 0.
  !
  IF (L_ukca_strat .OR. L_ukca_stratcfc .OR. L_ukca_strattrop) THEN

    ! LBC code contained within UKCA_SCENARIO_CTL
    CALL ukca_scenario_ctl(n_boundary_vals, lbc_spec,            &
                           i_year, i_day_number, lbc_mmr,        &
                           .NOT. (L_ukca_stratcfc))

  END IF    ! L_ukca_strat etc

  !       Calculate solar zenith angle

  ALLOCATE(cos_zenith_angle(row_length, rows))
  IF (.NOT. ALLOCATED(int_zenith_angle))                        &
     ALLOCATE(int_zenith_angle(row_length, rows))
  IF (.NOT. ALLOCATED(daylength))                               &
     ALLOCATE(daylength(row_length, rows))

  Eq_Time = 0.0
  Sindec  = 0.0
  scs     = 0.0

  secondssincemidnight = REAL(previous_time(4)*3600             &
                       +      previous_time(5)*60               &
                       +      previous_time(6))

  CALL solpos (previous_time(7), previous_time(1),              &
       secondssincemidnight, timestep, l_planet_obs ,           &
       Eq_Time, Sindec, scs, sindec_obs, eqt_obs)

  ! Compute COS(sza) integral only when day changes
  ! i.e. on first timestep of 'next' day, determined by checking
  ! if previous timestep was fully divisible by asteps_per_day.
  IF ( first .OR.                                                &
       MOD(timestep_number-1, asteps_per_day) == 0 ) THEN
    int_zenith_angle(:,:) = 0.0
    DO kk = 1, (A_Steps_per_hr*24)
      ssmn_incr = secondssincemidnight + timestep*(kk-1)
      CALL ukca_solang(Sindec, ssmn_incr,                   &
                       timestep, Eq_Time,                   &
                       sin_theta_latitude, true_longitude,  &
                       theta_field_size, cos_zenith_angle )

      WHERE (cos_zenith_angle > 0.0)
        int_zenith_angle(:,:) = int_zenith_angle(:,:) +     &
           cos_zenith_angle(:,:) * timestep
      END WHERE
    END DO
  END IF

  ! This call needs to be after above loop to prevent errors with value in
  !  cos_zenith_angle
  CALL ukca_solang(Sindec, secondssincemidnight,                &
              timestep, Eq_Time,                                &
              sin_theta_latitude, true_longitude,               &
              theta_field_size, cos_zenith_angle )

  IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
    ! Compute COS(sza) integral and daylength, also stores cos(zenith) for
    !  Offline chemistry diurnal cycle
    CALL ukca_int_cosz(row_length, rows, true_latitude, sindec, &
                       cos_zenith_angle,int_zenith_angle,       &
                       daylength)

    ! Read in netCDF offline oxidant fields
    CALL ukca_read_offline_oxidants_ctl(row_length, rows,       &
                        model_levels, global_row_length,        &
                        global_rows, timestep_number,           &
                        timestep, first)

  END IF

  ! Emissions system (based on NetCDF emission input in the case of
  ! non-interactive emissions)
  
  ! Read NetCDF emission fields and calculate online emissions, fill in
  ! the emissions structure, inject emissions, do tracer mixing and
  ! finally output emission diagnostics.
  CALL ukca_emiss_ctl (                                                 &
      row_length,   rows,                                               &
      bl_levels, global_row_length, global_rows,                        &
      n_chem_tracers + n_aero_tracers + n_mode_tracers,                 &
      i_year, i_month, i_day, i_hour, timestep_number, timestep, first, &
      land_points, land_index,                                          &
      conv_cloud_base       (1:row_length, 1:rows),                     &
      conv_cloud_top        (1:row_length, 1:rows),                     &
      delta_lambda, delta_phi,                                          &
      r_theta_levels        (1:row_length, 1:rows, 0:model_levels),     &
      FV_cos_theta_latitude (1:row_length, 1:rows),                     &
      cos_zenith_angle      (1:row_length, 1:rows),                     &
      int_zenith_angle      (1:row_length, 1:rows),                     &
      sin_theta_latitude    (1:row_length, 1:rows),                     &
      tan_theta_latitude    (1:row_length, 1:rows),                     &
      true_longitude        (1:row_length, 1:rows),                     &
      tropopause_height     (1:row_length, 1:rows),                     &
      f3u_by_omega          (1:row_length, 1:rows),                     &
      land_fraction         (1:row_length, 1:rows),                     &
      r_rho_levels          (1:row_length, 1:rows, 1:model_levels),     &
      p_theta_levels        (1:row_length, 1:rows, 1:model_levels),     &
      t_theta_levels        (1:row_length, 1:rows, 1:model_levels),     &
      p_layer_boundaries    (1:row_length, 1:rows, 0:model_levels),     &
      totnodens             (1:row_length, 1:rows, 1:model_levels),     &
      volume                (1:row_length, 1:rows, 1:model_levels),     &
      mass                  (1:row_length, 1:rows, 1:model_levels),     &
      interf_z              (1:row_length, 1:rows, 0:model_levels),     &
      land_sea_mask         (1:row_length, 1:rows),                     &
      theta                 (1:row_length, 1:rows, 1:model_levels),     &
      q                     (1:row_length, 1:rows, 1:model_levels),     &
      qcl                   (1:row_length, 1:rows, 1:model_levels),     &
      qcf                   (1:row_length, 1:rows, 1:model_levels),     &
      exner_rho_levels      (1:row_length, 1:rows, 1:model_levels+1),   &
      rho_r2                (1:row_length, 1:rows, 1:model_levels),     &
      kent                  (1:row_length, 1:rows),                     &
      kent_dsc              (1:row_length, 1:rows),                     &
      rhokh_mix             (1:row_length, 1:rows, 1:bl_levels),        &
      dtrdz_charney_grid    (1:row_length, 1:rows, 1:bl_levels),        &
      we_lim                (1:row_length, 1:rows, 1:3),                &
      t_frac                (1:row_length, 1:rows, 1:3),                &
      zrzi                  (1:row_length, 1:rows, 1:3),                &
      we_lim_dsc            (1:row_length, 1:rows, 1:3),                &
      t_frac_dsc            (1:row_length, 1:rows, 1:3),                &
      zrzi_dsc              (1:row_length, 1:rows, 1:3),                &
      zbl                   (1:row_length, 1:rows),                     &
      zhsc                  (1:row_length, 1:rows),                     &
      z_half                (1:row_length, 1:rows, 1:bl_levels),        &
      ch4_wetl_emiss        (1:row_length, 1:rows),                     &
      seaice_frac           (1:row_length, 1:rows),                     &
      area                  (1:row_length, 1:rows, 1:model_levels),     &
      dust_flux             (:,:,:),                                    &
      u_scalar_10m          (:,:),                                      &
      tstar                 (1:row_length, 1:rows),                     &
      dms_sea_conc          (1:row_length, 1:rows),                     &
      chloro_sea            (1:row_length, 1:rows),                     &
      all_tracers           (1:row_length, 1:rows, 1:model_levels,      &
                             1:n_chem_tracers+n_aero_tracers+           &
                               n_mode_tracers),  &
     STASH_maxlen(38,atmos_im), STASHwork38,                            &
     STASH_maxlen(50,atmos_im), STASHwork50)

  IF (L_asad_use_chem_diags .AND. L_asad_use_mass_diagnostic) &
       CALL asad_mass_diagnostic( &
       row_length, &
       rows, &
       model_levels, &
       mass, &
       ierr)

! ----------------------------------------------------------------------
! 4.5 Call chemistry routines
! ----------------------------------------------------------------------
  ! Do chemistry calculation here (at chemistry timesteps)
  IF (do_chemistry) THEN

    IF (.NOT. ALLOCATED(t_chem)) ALLOCATE(t_chem(row_length,rows,model_levels))
    IF (.NOT. ALLOCATED(q_chem)) ALLOCATE(q_chem(row_length,rows,model_levels))
    t_chem(:,:,:) = t_theta_levels(1:row_length,1:rows,1:model_levels)
    q_chem(:,:,:) = q(1:row_length,1:rows,1:model_levels)

    !======================================
    ! PASSIVE OZONE - CheST/StratChem ONLY
    !  n_passive is set in ukca_calc_cspecies
             ! Do passive ozone
    IF (n_passive > 0) THEN
       ! assumes 360 day calendar
      IF ((i_day_number == 91) .OR. (i_day_number == 301)) THEN
        IF (i_hour == 0) THEN
          all_tracers(:,:,1:model_levels,n_passive) =          &
                all_tracers(:,:,1:model_levels,n_o3)
        END IF
      END IF
    END IF
    !======================================
    l_store_value = .TRUE.
    IF (L_asad_use_chem_diags .AND. L_asad_use_tendency) &
         CALL asad_tendency_STE( &
         row_length,rows,model_levels, &
         n_chem_tracers,all_tracers, &
         totnodens,volume,mass,timestep, &
         calculate_tendency,l_store_value,ierr)


    IF (interval > 1) THEN
       ! MAY NEED TO CALL A SECOND TIME TO HAVE CORRECT MEANING PERIOD
      CALL ukca_solang(Sindec, secondssincemidnight,              &
           REAL(chem_timestep), Eq_Time,                          &
           sin_theta_latitude, true_longitude,                    &
           theta_field_size, cos_zenith_angle )
    END IF

    IF (l_ukca_mode) THEN
      IF (L_ukca_aerchem) THEN
        ! SO2 wet oxidation and H2SO4 updating done in UKCA
        wetox_in_aer = 0
        uph2so4inaer = 0
      ELSE IF (l_ukca_nr_aqchem .OR. l_ukca_offline_be .OR.       &
               l_ukca_raqaero) THEN
        ! SO2 wet oxidation in UKCA, and H2SO4 updating done in MODE
        wetox_in_aer = 0
        uph2so4inaer = 1
      ELSE                        ! No SO2 oxidation in UKCA chemistry
        wetox_in_aer = 1          ! Sulphur emissions are still needed
        uph2so4inaer = 1
      END IF
      IF (printstatus >= prstatus_oper) THEN
        WRITE(umMessage,'(A,I5)') 'uph2so4inaer is set to: ',      &
              uph2so4inaer
        CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
        WRITE(umMessage,'(A,I5)') 'wetox_in_aer is set to: ',      &
              wetox_in_aer
        CALL umPrint(umMessage,src='ukca_main1-ukca_main1')
      END IF
    ELSE
      wetox_in_aer = 0
      uph2so4inaer = 0
    END IF                 ! l_ukca_mode

    IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
      ! Allocate diagnostic arrays for offline chemistry
      IF (.NOT. ALLOCATED(o3_offline_diag))                         &
        ALLOCATE(o3_offline_diag(theta_field_size,1:model_levels))
      IF (.NOT. ALLOCATED(oh_offline_diag))                         &
        ALLOCATE(oh_offline_diag(theta_field_size,1:model_levels))
      IF (.NOT. ALLOCATED(no3_offline_diag))                        &
        ALLOCATE(no3_offline_diag(theta_field_size,1:model_levels))
      IF (.NOT. ALLOCATED(ho2_offline_diag))                        &
        ALLOCATE(ho2_offline_diag(theta_field_size,1:model_levels))
      ! Initialise to zero
      o3_offline_diag(:,:) = 0.0
      oh_offline_diag(:,:) = 0.0
      no3_offline_diag(:,:) = 0.0
      ho2_offline_diag(:,:) = 0.0
    END IF

    IF (l_ukca_offline_be) THEN
! Offline chemistry with explicit backward-Euler solver
       CALL ukca_chemistry_ctl_be(i_month, i_day_number, i_hour,      &
            r_minute - timestep/60.,                                  &
            chem_timestep,                                            &
            k_be_top,                                                 &
            n_chem_tracers+n_aero_tracers,                            &
            f3u_by_omega(1:row_length,1:rows),                        &
            fv_cos_theta_latitude(1:row_length,1:rows),               &
            true_longitude,                                           &
            p_theta_levels(1:row_length,1:rows,1:model_levels),       &
            t_chem(1:row_length,1:rows,1:model_levels),               &
            q_chem(1:row_length,1:rows,1:model_levels),               &
            qcl(1:row_length, 1:rows,1:model_levels),                 &
            rel_humid_frac(1:row_length,1:rows,1:model_levels),       &
            p_layer_boundaries(1:row_length,1:rows,0:model_levels),   &
            r_theta_levels(1:row_length, 1:rows, :),                  &
            all_tracers(1:row_length,1:rows,1:model_levels,           &
            1:n_chem_tracers+n_aero_tracers),                         &
            tstar,                                                    &
            thick_bl_levels,                                          &
            rough_length,                                             &
            u_s,                                                      &
            ls_ppn3d, conv_ppn3d,                                     &
            cloud_frac(1:row_length, 1:rows, :),                      &
            volume(:,:,:),                                            &
            mass(:,:,:),                                              &
            ! Extra variables for new dry dep scheme
            land_points, land_index,                                  &
            tile_pts, tile_index, tile_frac,                          &
            zbl, surf_hf, seaice_frac, stcon,                         &
            soil_moisture_layer1, fland,                              &
            laift_lp, canhtft_lp,                                     &
            z0tile_lp, tstar_tile, canwctile_lp,                      &
            pv_on_theta_mlevs,                                        &
            theta(1:row_length,1:rows,1:model_levels),                &
            uph2so4inaer,                                             &
            delso2_wet_h2o2,                                          &
            delso2_wet_o3,                                            &
            delh2so4_chem,                                            &
            delso2_drydep,                                            &
            delso2_wetdep                                             &
            )

    ELSE
       IF (l_ukca_asad_columns .AND. (.NOT. (L_ukca_trop .OR.         &
            L_ukca_aerchem .OR. L_ukca_raq .OR. l_ukca_raqaero))) THEN  
          CALL ukca_chemistry_ctl_col(i_month, i_day_number, i_hour,  &
               r_minute - timestep/60.0,                              &
               REAL(chem_timestep),                                   &
               n_chem_tracers+n_aero_tracers,                         &
               f3u_by_omega(1:row_length,1:rows),                     &
               FV_cos_theta_latitude(1:row_length,1:rows),            &
               true_longitude,                                        &
               p_theta_levels(1:row_length,1:rows,1:model_levels),    &
               T_chem(1:row_length,1:rows,1:model_levels),            &
               Q_chem(1:row_length,1:rows,1:model_levels),            &
               qcf(1:row_length,1:rows,1:model_levels),               &
               qcl(1:row_length, 1:rows,1:model_levels),              &
               rel_humid_frac(1:row_length,1:rows,1:model_levels),    &
               p_layer_boundaries(1:row_length,1:rows,0:model_levels),&
               r_theta_levels(1:row_length, 1:rows, :),               &
               z_top_of_model,                                        &
               cos_zenith_angle,                                      &
               all_tracers(1:row_length,1:rows,1:model_levels,        &
               1:n_chem_tracers+n_aero_tracers),             &
               all_ntp,                                               &
               Tstar,                                                 &
               Thick_bl_levels,                                       &
               Rough_length,                                          &
               u_s,                                                   &
               ls_ppn3d, conv_ppn3d,                                  &
               cloud_frac(1:row_length, 1:rows, :),                   &
               dj(:,:,:,:),                                           &
               volume(:,:,:),                                         &
               mass(:,:,:),                                           &
               ! Extra variables for new dry dep scheme
               land_points, land_index,                               &
               tile_pts, tile_index, tile_frac,                       &
               zbl, surf_hf, seaice_frac, stcon,                      &
               soil_moisture_layer1, fland,                           &
               laift_lp, canhtft_lp,                                  &
               z0tile_lp, tstar_tile, canwctile_lp,                   &
               have_nat3d,                                            &
               pv_on_theta_mlevs,                                     &
               theta(1:row_length,1:rows,1:model_levels),             &
               um_ozone3d,                                            &
               uph2so4inaer,                                          &
               delso2_wet_h2o2,                                       &
               delso2_wet_o3,                                         &
               delh2so4_chem,                                         &
               delso2_drydep,                                         &
               delso2_wetdep,                                         &
               so4_sa,                                                &
               ! Diagnostics
               nat_psc(1:row_length,1:rows,1:model_levels),           &
               atm_ch4_mol(1:row_length,1:rows,1:model_levels),       &
               atm_co_mol(1:row_length,1:rows,1:model_levels),        &
               atm_n2o_mol(1:row_length,1:rows,1:model_levels),       &
               atm_cf2cl2_mol(1:row_length,1:rows,1:model_levels),    &
               atm_cfcl3_mol(1:row_length,1:rows,1:model_levels),     &
               atm_mebr_mol(1:row_length,1:rows,1:model_levels),      &
               atm_h2_mol(1:row_length,1:rows,1:model_levels)         &
               )

       ELSE
          ! Aerosol mmr / numbers (from CLASSIC aerosol scheme) are only
          ! used by the RAQ chemistry scheme if heterogeneous chemistry
          ! is ON. Even if that functionality is OFF they need to be
          ! allocated before the call to UKCA_CHEMISTRY_CTL. Note that
          ! so4_aitken and so4_accum (also used for heterogeneous
          ! chemistry) are allocated somewhere else in this routine.
          IF (.NOT. ALLOCATED(soot_fresh))                            &
              ALLOCATE (soot_fresh (row_length, rows, model_levels))

          IF (.NOT. ALLOCATED(soot_aged))                             &
              ALLOCATE (soot_aged  (row_length, rows, model_levels))

          IF (.NOT. ALLOCATED(ocff_fresh))                            &
              ALLOCATE (ocff_fresh (row_length, rows, model_levels))

          IF (.NOT. ALLOCATED(ocff_aged))                             &
              ALLOCATE (ocff_aged  (row_length, rows, model_levels))

          IF (.NOT. ALLOCATED(biogenic))                              &
              ALLOCATE (biogenic   (row_length, rows, model_levels))

          IF (.NOT. ALLOCATED(sea_salt_film))                         &
              ALLOCATE (sea_salt_film (row_length, rows, model_levels))

          IF (.NOT. ALLOCATED(sea_salt_jet))                          &
              ALLOCATE (sea_salt_jet  (row_length, rows, model_levels))

          CALL ukca_chemistry_ctl(i_month, i_day_number, i_hour,      &
               r_minute - timestep/60.0,                              &
               REAL(chem_timestep),                                   &
               n_chem_tracers+n_aero_tracers,                         &
               f3u_by_omega(1:row_length,1:rows),                     &
               FV_cos_theta_latitude(1:row_length,1:rows),            &
               true_longitude,                                        &
               p_theta_levels(1:row_length,1:rows,1:model_levels),    &
               T_chem(1:row_length,1:rows,1:model_levels),            &
               Q_chem(1:row_length,1:rows,1:model_levels),            &
               qcf(1:row_length,1:rows,1:model_levels),               &
               qcl(1:row_length, 1:rows,1:model_levels),              &
               rel_humid_frac(1:row_length,1:rows,1:model_levels),    &
               p_layer_boundaries(1:row_length,1:rows,0:model_levels),&
               r_theta_levels(1:row_length, 1:rows, :),               &
               z_top_of_model,                                        &
               cos_zenith_angle,                                      &
               all_tracers(1:row_length,1:rows,1:model_levels,        &
               1:n_chem_tracers+n_aero_tracers),                      &
               all_ntp,                                               &
               Tstar,                                                 &
               Thick_bl_levels,                                       &
               Rough_length,                                          &
               u_s,                                                   &
               ls_ppn3d, conv_ppn3d,                                  &
               cloud_frac(1:row_length, 1:rows, :),                   &
               dj(:,:,:,:),                                           &
               volume(:,:,:),                                         &
               mass(:,:,:),                                           &
               so4_aitken    (1:row_length, 1:rows, 1:model_levels),  &
               so4_accum     (1:row_length, 1:rows, 1:model_levels),  &
               soot_fresh    (1:row_length, 1:rows, 1:model_levels),  &
               soot_aged     (1:row_length, 1:rows, 1:model_levels),  &
               ocff_fresh    (1:row_length, 1:rows, 1:model_levels),  &
               ocff_aged     (1:row_length, 1:rows, 1:model_levels),  &
               biogenic      (1:row_length, 1:rows, 1:model_levels),  &
               sea_salt_film (1:row_length, 1:rows, 1:model_levels),  &
               sea_salt_jet  (1:row_length, 1:rows, 1:model_levels),  &
               ! Extra variables for new dry dep scheme
               land_points, land_index,                               &
               tile_pts, tile_index, tile_frac,                       &
               zbl, surf_hf, seaice_frac, stcon,                      &
               soil_moisture_layer1, fland,                           &
               laift_lp, canhtft_lp,                                  &
               z0tile_lp, tstar_tile, canwctile_lp,                   &
               have_nat3d,                                            &
               pv_on_theta_mlevs,                                     &
               theta(1:row_length,1:rows,1:model_levels),             &
               um_ozone3d,                                            &
               uph2so4inaer,                                          &
               delso2_wet_h2o2,                                       &
               delso2_wet_o3,                                         &
               delh2so4_chem,                                         &
               delso2_drydep,                                         &
               delso2_wetdep,                                         &
               so4_sa,                                                &
               ! Diagnostics
               nat_psc(1:row_length,1:rows,1:model_levels),           &
               trop_ch4_mol(1:row_length,1:rows,1:model_levels),      &
               trop_o3_mol(1:row_length,1:rows,1:model_levels),       &
               trop_oh_mol(1:row_length,1:rows,1:model_levels),       &
               strat_ch4_mol(1:row_length,1:rows,1:model_levels),     &
               strat_ch4loss(1:row_length,1:rows,1:model_levels),     &
               atm_ch4_mol(1:row_length,1:rows,1:model_levels),       &
               atm_co_mol(1:row_length,1:rows,1:model_levels),        &
               atm_n2o_mol(1:row_length,1:rows,1:model_levels),       &
               atm_cf2cl2_mol(1:row_length,1:rows,1:model_levels),    &
               atm_cfcl3_mol(1:row_length,1:rows,1:model_levels),     &
               atm_mebr_mol(1:row_length,1:rows,1:model_levels),      &
               atm_h2_mol(1:row_length,1:rows,1:model_levels),        &
               STASH_maxlen(50,atmos_im),                             &
               STASHwork50                                            &
               )
       END IF
    END IF

    IF (ALLOCATED(cos_zenith_angle)) DEALLOCATE(cos_zenith_angle)
    IF (ALLOCATED(T_chem)) DEALLOCATE(T_chem)
    IF (ALLOCATED(Q_chem)) DEALLOCATE(Q_chem)

    ! Add fields to section 50 diagnostics if requested
    IF (l_ukca_offline .OR. l_ukca_offline_be) THEN
      CALL ukca_offline_oxidants_diags(row_length, rows, model_levels,   &
                                       totnodens,                        &
                                    STASH_maxlen(50,atmos_im),STASHwork50)
    END IF

  END IF ! do_chemistry

  IF (n_strat_fluxdiags > 0) THEN
    DO l=1,n_chem_tracers
      DO k=1,model_levels
        DO j=1,rows
          DO i=1,row_length
            trmol_post_chem(i,j,k,l) = all_tracers(i,j,k,l)        &
                                *totnodens(i,j,k)*volume(i,j,k)    &
                                /(c_species(l)*avc)    ! moles
          END DO
        END DO
      END DO
    END DO
  END IF   ! End of IF n_strat_fluxdiags > 0 statement

  !  Save all_tracers array for STE calculation in ASAD Flux Diagnostics
  IF (.NOT. ALLOCATED(fluxdiag_all_tracers) .AND. &
       (L_asad_use_chem_diags .AND. L_asad_use_STE)) &
       THEN
    ALLOCATE(fluxdiag_all_tracers(&
         1-offx:row_length+offx,1-offy:rows+offy, &
         1:model_levels,1:n_chem_tracers))
    fluxdiag_all_tracers(:,:,1:model_levels,:) =                  &
         all_tracers(:,:,1:model_levels,1:n_chem_tracers)
  END IF

  l_store_value=.FALSE.
  IF (do_chemistry .AND. (L_asad_use_chem_diags                   &
       .AND. L_asad_use_tendency))                                &
    CALL asad_tendency_ste(row_length,rows,model_levels,          &
                           n_chem_tracers,all_tracers,            &
                           totnodens,volume,mass,timestep,        &
                           calculate_tendency,l_store_value,ierr)

  ! Transform halogen/nitrogen/hydrogen species back
  IF ((L_ukca_strat) .OR. (L_ukca_strattrop) .OR. (L_ukca_stratcfc)) THEN

    ! first, copy unlumped NO2, BrO and HCl fields into all_ntp structure

    ! NO2
    i = name2ntpindex(all_ntp,'NO2       ')
    all_ntp(i)%data_3d(1:row_length,1:rows,1:model_levels) =             &
         all_tracers(1:row_length,1:rows,1:model_levels,n_no2)

    ! BrO
    i = name2ntpindex(all_ntp,'BrO       ')
    all_ntp(i)%data_3d(1:row_length,1:rows,1:model_levels) =             &
         all_tracers(1:row_length,1:rows,1:model_levels,n_bro)

    ! HCl
    i = name2ntpindex(all_ntp,'HCl       ')
    all_ntp(i)%data_3d(1:row_length,1:rows,1:model_levels) =             &
         all_tracers(1:row_length,1:rows,1:model_levels,n_hcl)

    ! Now call ukca_transform halogen to do the lumping. After this,
    ! the all_tracers array contains the lumped species not NO2, BrO and HCl
    CALL ukca_transform_halogen(tr_ukca,rows,row_length,            &
                    model_levels, tdims_s%halo_i, tdims_s%halo_j,   &
                    all_tracers(:,:,1:model_levels,:),              &
                    tdims_s%halo_i, tdims_s%halo_j,                 &
                    q(:,:,1:model_levels), lump_species, 0)
  END IF

  IF (ltimer) CALL timer('UKCA CHEMISTRY MODEL',6)
END IF    ! End of IF (L_ukca_chem) statement

! ----------------------------------------------------------------------
! 4.6 GLOMAP-mode aerosol scheme
! ----------------------------------------------------------------------
IF (l_ukca_mode .AND. do_chemistry) THEN

  ! Allocate space for copies of fields required for later diagnostic
  ! calculations
  IF (l_ukca_cmip6_diags .OR. l_ukca_pm_diags)                       &
    CALL ukca_mode_diags_alloc(theta_field_size*model_levels) 

  ! Cache blocking scheme
  stride_seg = row_length*rows    ! i.e. the number of columns on MPI task
  nbox_max = ukca_mode_seg_size*model_levels   ! requested cols by levels
  IF ( (ukca_mode_seg_size < 1).OR.(ukca_mode_seg_size > stride_seg) ) THEN
    ! unrealistic value for ukca_mode_seg_size more than columns on MPI task
    WRITE(cmessage,'(A,i8,A,i8)') 'Select ukca_mode_seg_size as a factor of'// &
                 ' stride_seg. Current values are: ukca_mode_seg_size=',    &
                        ukca_mode_seg_size, ' and stride_seg=',stride_seg
    errcode = 6
    CALL ereport(routinename,errcode,cmessage)
    ! code will abort (no way to test any earlier)
  ELSE
    ! divide number of total columns by the ncpc
    ! this allows segments to cross rows.
    nseg = stride_seg/ukca_mode_seg_size             ! number of segments 
                                                     ! (check remainder)
    nbox_s(1:nseg) = ukca_mode_seg_size*model_levels ! all uniform initially
    nbox_s(nseg+1:stride_seg) = 0
    ncol_s(1:nseg) = ukca_mode_seg_size 
    ncol_s(nseg+1:stride_seg) = 0

    seg_rem = MOD(stride_seg,ukca_mode_seg_size)     ! Only allow whole 
                                                     ! factors of cols
    IF (seg_rem > 0) THEN 
      ! account for final additional segment
      nseg = nseg+1
      ncol_s(nseg) = seg_rem    ! a smaller number of columns
      nbox_s(nseg) = seg_rem*model_levels  ! fewer boxes on final segment
      IF ( first ) THEN
        WRITE(cmessage,'(A,i8,A,i8,A)')                                        &
       'UKCA_MAIN1: ukca_mode_seg_size is not a factor of columns on this PE ', &
          stride_seg,'. Segments will be non-uniform with final segment of ',  &
          seg_rem,' columns.'
        errcode = -6
        CALL ereport(routinename,errcode,cmessage)
      END IF  ! first
    END IF
    
    ! determine start of each segment in terms of column base
    lb = 1
    DO ik = 1, nseg 
      lbase(ik) = lb 
      lb = lb + ncol_s(ik)
    END DO

  END IF

  IF (ltimer) CALL timer('UKCA AEROSOL MODEL  ',5)

  CALL ukca_aero_ctl(i_month, i_day_number, i_hour,               &
       INT(r_minute - timestep/60.0),                             &
       REAL(chem_timestep),                                       &
       rows, row_length,                                          &
       global_row_length,global_rows,                             &
       n_chem_tracers+n_aero_tracers,                             &
       n_mode_tracers,                                            &
       area,                                                      &
       f3u_by_omega(1:row_length,1:rows),                         &
       FV_cos_theta_latitude(1:row_length,1:rows),                &
       true_longitude,                                            &
       p_theta_levels(1:row_length,1:rows,1:model_levels),        &
       t_theta_levels(1:row_length,1:rows,:),                     &
       q(1:row_length,1:rows,1:model_levels),                     &
       rel_humid_frac(1:row_length,1:rows,:),                     &
       rel_humid_frac_clr(1:row_length,1:rows,:),                 &
       p_layer_boundaries(1:row_length,1:rows,0:model_levels),    &
       all_tracers(1:row_length,1:rows,1:model_levels,            &
            1:n_chem_tracers+n_aero_tracers),                     &
       all_tracers(1:row_length,1:rows,1:model_levels,            &
            n_chem_tracers+n_aero_tracers+1:                      &
            n_chem_tracers+n_aero_tracers+                        &
            n_mode_tracers),                                      &
       Tstar(1:row_length,1:rows),                                &
       seaice_frac(1:row_length,1:rows),                          &
       Rough_length(1:row_length,1:rows),                         &
       u_s,                                                       &
       ls_rain3d(1:row_length,1:rows,1:model_levels),             &
       conv_rain3d(1:row_length,1:rows,1:model_levels),           &
       ls_snow3d(1:row_length,1:rows,1:model_levels),             &
       conv_snow3d(1:row_length,1:rows,1:model_levels),           &
       autoconv(1:row_length,1:rows,1:model_levels),              &
       accretion(1:row_length,1:rows,1:model_levels),             &
       rim_agg(1:row_length,1:rows,1:model_levels),               &
       rim_cry(1:row_length,1:rows,1:model_levels),               &
       land_fraction,                                             &
       nbox_max,                                                  &
       delso2_wet_h2o2,                                           &
       delso2_wet_o3,                                             &
       delh2so4_chem,                                             &
       mode_diags,                                                &
       cloud_frac(1:row_length,1:rows,1:model_levels),            &
       cloud_liq_frac(1:row_length,1:rows,1:model_levels),        &
       qcl(1:row_length,1:rows,1:model_levels),                   &
       offx, offy,                                                &
       z_half_alllevs, delta_r,                                   &
       volume,mass,zbl,                                           &
       uph2so4inaer,                                              &
       wetox_in_aer,                                              &
       all_ntp,                                                   &
       nseg, nbox_s, ncol_s, lbase, stride_seg                    &
       )

  IF (ltimer) CALL timer('UKCA AEROSOL MODEL  ',6)

  ! Call activation scheme if switched on in namelist
  IF (L_ukca_arg_act) THEN

    CALL ukca_activate(                                           &
        row_length, rows,                                         &
        bl_levels,                                                &
        theta_field_size,                                         &
        n_mode_tracers,                                           &
        n_mode_diags,                                             &
        p_theta_levels(1:row_length,1:rows,1:model_levels),       &
        t_theta_levels(1:row_length,1:rows,1:model_levels),       &
        q(1:row_length,1:rows,1:model_levels),                    &
        qsvp(1:row_length,1:rows,1:model_levels),                 &
        bl_tke(1:row_length,1:rows,1:bl_levels),                  &
        vertvel(1:row_length,1:rows,0:model_levels-1),            &
        cloud_liq_frac(1:row_length,1:rows,1:model_levels),       &
        qcl(1:row_length,1:rows,1:model_levels),                  &
        cdncflag,                                                 &
        cdnc3flag,                                                &
        cdncwt,                                                   &
        all_tracers(1:row_length,1:rows,1:model_levels,           &
                    n_chem_tracers+n_aero_tracers+1:              &
                    n_chem_tracers+n_aero_tracers+n_mode_tracers),&
        mode_diags                                                &
        )
    
    ! write CDNC and <CDNC^-1/3> to all_ntp structure
    
    i = name2ntpindex(all_ntp,'cdnc      ')
    all_ntp(i)%data_3d(:,:,:)=cdncflag(:,:,:)
    i = name2ntpindex(all_ntp,'cdnc3     ')
    all_ntp(i)%data_3d(:,:,:)=cdnc3flag(:,:,:)
    
  END IF   ! L_ukca_arg_act

END IF   ! L_ukca_mode

! ----------------------------------------------------------------------
! 5. Output prognostics to D1 and diagnostics to stash
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 5.1 Put prognostics into D1
! ----------------------------------------------------------------------
! Take a copy of tracer fields if requested for diagnostic purposes
IF (L_asad_use_chem_diags .AND. L_asad_use_output_tracer)         &
     CALL asad_output_tracer(row_length, rows, model_levels,      &
                        n_chem_tracers,                           &
                        all_tracers(:,:,1:model_levels,:), ierr)


! For ENDGAME Copy level1 to level0 before QPOS/copy to D1
all_tracers(:,:,0,:) = all_tracers(:,:,1,:)
IF ( L_ukca_chem ) q(:,:,0) =  q(:,:,1)

!     Return fields to D1
CALL putd1flds()

! Update the tracer_ukca_um and q_um arrays from the data used in 
! chemistry and aerosol modules. These are then passed back 
! to the UM.
CALL ukca_all_tracers_copy_out(tracer_ukca_um, q_um)

! ----------------------------------------------------------------------
! 5.2 Add output fields to STASHwork arrays.
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 5.2.1 MODE diagnostics  (items 38,201 - 38,212 now done in ukca_mode_ems_um)
! ----------------------------------------------------------------------

IF ( l_ukca_chem )  THEN
im_index = 1

    icnt=0
    DO l=1,nmax_mode_diags
      IF (UkcaD1codes(imode_first+l-1)%item /= imdi) THEN
        icnt = icnt + 1
        item = UkcaD1codes(imode_first+l-1)%item
        section = stashcode_glomap_sec
        IF (sf(item,section) .AND. item > item1_mode_diags+12 .AND.   &
            item < stashcode_so4_nuc_sol) THEN
          ! DEPENDS ON: copydiag_3d
          CALL copydiag_3d (stashwork38(si(item,section,im_index)),   &
            mode_diags(:,:,:,icnt),                                   &
            row_length,rows,model_levels,0,0,0,0, at_extremity,       &
            stlist(1,stindex(1,item,section,im_index)),len_stlist,    &
            stash_levels,num_stash_levels+1,                          &
            atmos_im,section,item,icode,cmessage)

          IF (icode >  0) THEN
            icode = section*1000+item
            CALL ereport('UKCA_MAIN1',icode,cmessage)
          END IF
        END IF
      END IF
    END DO       ! 1,nmax_mode_diags

    ! Copy CMIP6 diagnostics and/or PM diagnostics into STASHwork array
    IF (l_ukca_mode .AND. do_chemistry .AND.                                   &
        (l_ukca_cmip6_diags .OR. l_ukca_pm_diags)) THEN
      CALL ukca_mode_diags(theta_field_size*model_levels,                      &
                           n_mode_tracers,                                     &
                           p_theta_levels(1:row_length,1:rows,1:model_levels), &
                           t_theta_levels(1:row_length,1:rows,:),              &
                           all_tracers(1:row_length,1:rows,1:model_levels,     &
                                       n_chem_tracers+n_aero_tracers+1:        &
                                       n_chem_tracers+n_aero_tracers+          &
                                       n_mode_tracers),                        &
                           interf_z,                                           &
                           STASH_maxlen(38,atmos_im),                          &
                           STASHwork38)
    END IF

! ------------------------------
! Plume Scavenging Diagnostics
! ------------------------------

    ! Copy plume scavenging diagnostics into STASHwork array
    IF (l_ukca_mode .AND. l_ukca_plume_diags) THEN
      section = stashcode_glomap_sec
      DO imode=1,nmodes
        IF (mode(imode)) THEN
          DO icp = 1,ncp
            IF (component(imode,icp)) THEN
              item = item_no_plume(imode,icp) - msect38   ! find stash item no.
              IF (item < 0) THEN
                icode = 1
                cmessage = ' No valid item number identified for this mode'// &
                           ' and component'
                WRITE(umMessage,'(A70,2I6)') cmessage,imode,icp
                CALL umPrint(umMessage,src='ukca_main1')
                CALL ereport('UKCA_MAIN1',icode,cmessage)
              END IF

              IF (sf(item,section)) THEN
                ! Set plume_scav_diags for this mode and component
                CALL ukca_plume_scav_diags_2d(imode, icp,row_length, rows, &
                                        model_levels, tr_ukca)

                ! DEPENDS ON: copydiag
                CALL copydiag(stashwork38(si(item,section,im_index)),        &
                     plume_scav_diag_2d(:,:),                                &
                     row_length,rows,0,0,0,0, at_extremity,                  &
                     atmos_im,section,item,                                  &
                     icode,cmessage)

                IF (icode >  0) THEN
                  icode = section*1000+item
                  CALL ereport('UKCA_MAIN1',icode,cmessage)
                END IF
              END IF       ! sf
            END IF        ! component
          END DO         ! icp
        END IF          ! mode(imode)
      END DO           ! imode

    ! Set convective indicator to false ready for next step.
      plume_scav_diag_ping = .FALSE.
      IF (ALLOCATED(d_tracers)) DEALLOCATE(d_tracers)
      IF (ALLOCATED(plume_scav_diag_2d)) DEALLOCATE(plume_scav_diag_2d)
    END IF            ! l_ukca_mode ...

! ----------------------------------------------------------------------
! 5.2.2 Stratospheric flux diagnostics
! ----------------------------------------------------------------------
    icnt = 0
    DO l=1,nmax_strat_fluxdiags
      IF (UkcaD1codes(istrat_first+l-1)%item /= imdi) THEN
        icnt = icnt + 1
        item = UkcaD1codes(istrat_first+l-1)%item
        section = stashcode_glomap_sec
        IF (sf(item,section)) THEN
          ! DEPENDS ON: copydiag_3d
          CALL copydiag_3d (stashwork38(si(item,section,im_index)),   &
            strat_fluxdiags(:,:,:,icnt),                              &
            row_length,rows,model_levels,0,0,0,0, at_extremity,       &
            stlist(1,stindex(1,item,section,im_index)),len_stlist,    &
            stash_levels,num_stash_levels+1,                          &
            atmos_im,section,item,icode,cmessage)

          IF (icode >  0) THEN
            icode = section*1000+item
            CALL ereport('UKCA_MAIN1',icode,cmessage)
          END IF
        ELSE
          cmessage=' Strat flux item not found by STASH flag array'
          icode = section*1000+item
          CALL ereport('UKCA_MAIN1',icode,cmessage)
        END IF
      END IF
    END DO       ! 1,nmax_strat_fluxdiags

! ----------------------------------------------------------------------
! 5.2.3 diagnostics calculated from arrays passed back from
! chemistry_ctl
! ----------------------------------------------------------------------

    CALL ukca_chem_diags_allts(                                        &
      row_length, rows, model_levels, n_use_tracers,                   &
      p_tropopause(1:row_length,1:rows),                               &
      all_tracers(1:row_length,1:rows,1:model_levels,1:n_use_tracers), &
      p_theta_levels(1:row_length,1:rows,1:model_levels),              &
      p_layer_boundaries(1:row_length,1:rows,0:model_levels),          &
      z_top_of_model,                                                  &
      STASH_maxlen(50,atmos_im), STASHwork50)

    IF (do_chemistry) &
      CALL ukca_chem_diags(                                &
        row_length, rows, model_levels,                    &
        nat_psc(1:row_length,1:rows,1:model_levels),       &
        trop_ch4_mol(1:row_length,1:rows,1:model_levels),  &
        trop_o3_mol(1:row_length,1:rows,1:model_levels),   &
        trop_oh_mol(1:row_length,1:rows,1:model_levels),   &
        strat_ch4_mol(1:row_length,1:rows,1:model_levels), &
        strat_ch4loss(1:row_length,1:rows,1:model_levels), &
        atm_ch4_mol(1:row_length,1:rows,1:model_levels),   &
        atm_co_mol(1:row_length,1:rows,1:model_levels),    &
        atm_n2o_mol(1:row_length,1:rows,1:model_levels),   &
        atm_cf2cl2_mol(1:row_length,1:rows,1:model_levels),&
        atm_cfcl3_mol(1:row_length,1:rows,1:model_levels), &
        atm_mebr_mol(1:row_length,1:rows,1:model_levels),  &
        atm_h2_mol(1:row_length,1:rows,1:model_levels),    &
        so4_sa(1:row_length,1:rows,1:model_levels),        &
        dj (1:row_length,1:rows,1:model_levels,1:jppj),    &
        STASH_maxlen(50,atmos_im), STASHwork50)

! ----------------------------------------------------------------------
! 5.2.4 ASAD flux diagnostics: fill STASHwork - needs to be done
!         BEFORE ste calc below
! ----------------------------------------------------------------------
    IF (L_asad_use_chem_diags) THEN
      stashsize = SIZE(STASHwork50)
      CALL asad_flux_put_stash(                                      &
           stashsize,STASHwork50,                                    &
           row_length,rows,model_levels,do_chemistry)
    END IF


! ----------------------------------------------------------------------
! 5.2.5 ASAD Flux Diags Strat-Trop Exchange - done AFTER asad_flux_put_stash
! ----------------------------------------------------------------------
    l_store_value = .TRUE.
    IF (L_asad_use_chem_diags .AND. L_asad_use_STE) THEN
      CALL asad_tendency_STE(                                         &
                   row_length,rows,model_levels,                      &
                   n_chem_tracers,fluxdiag_all_tracers,               &
                   totnodens,volume,mass,timestep,                    &
                   calculate_STE,l_store_value,ierr)
    END IF
END IF ! l_ukca_chem

! ----------------------------------------------------------------------
! 5.2.6 Grid-cell volume. Needs to be outside of IF(ukca_chem) section
! ----------------------------------------------------------------------
section = ukca_diag_sect
item = 255
IF (sf(item,ukca_diag_sect) .AND. ALLOCATED(volume)) THEN
  ! DEPENDS ON: copydiag_3d
  CALL copydiag_3d (stashwork50(si(item,section,im_index)),       &
        volume(:,:,:),                                            &
        row_length,rows,model_levels,0,0,0,0, at_extremity,       &
        stlist(1,stindex(1,item,section,im_index)),len_stlist,    &
        stash_levels,num_stash_levels+1,                          &
        atmos_im,section,item,icode,cmessage)

  IF (icode >  0) THEN
    icode = section*1000+item
    CALL ereport('UKCA_MAIN1',icode,cmessage)
  END IF

END IF  ! SF()

! ----------------------------------------------------------------------
! 5.3 Processing Diagnostics on pressure levels: section 34 to section 51
! and section 50 to section 52
! ----------------------------------------------------------------------
IF (l_ukca_chem_plev .OR. (l_ukca_chem .AND. l_ukca_asad_plev)) THEN
    CALL ukca_plev_diags(l_ukca_chem, l_ukca_chem_plev, l_ukca_asad_plev, &
         row_length, rows, model_levels,                                  &
         all_ntp, exner_theta_levels, tracer_ukca_um, STASHwork50)
END IF      ! l_ukca_asad_plev

! ----------------------------------------------------------------------
! 5.4 Call stash routine for sections 34, 38 and 50.
! ----------------------------------------------------------------------
! DEPENDS ON: stash
CALL stash(atmos_sm,atmos_im,34,STASHwork34,                     &
  icode,cmessage)

IF (icode >  0) THEN
  cmessage=" Error in STASH"//cmessage
  errcode=34
  CALL ereport('UKCA_MAIN1',errcode,cmessage)
END IF

! DEPENDS ON: stash
CALL stash(atmos_sm,atmos_im,38,STASHwork38,                     &
    icode,cmessage)

IF (icode >  0) THEN
  cmessage=" Error in STASH"//cmessage
  errcode=38
  CALL ereport('UKCA_MAIN1',errcode,cmessage)
END IF

! DEPENDS ON: stash
CALL stash(atmos_sm,atmos_im,50,STASHwork50,                     &
    icode,cmessage)

IF (icode >  0) THEN
  cmessage=" Error in STASH"//cmessage
  errcode=50
  CALL ereport('UKCA_MAIN1',errcode,cmessage)
END IF

! ----------------------------------------------------------------------
! 6. Finally deallocate arrays
! ----------------------------------------------------------------------
! These are currently only allocated on first call - see ukca_iniasad
! CALL ASAD_MOD_FINAL

! Deallocate variables in ukca_uk_interf_mod
CALL ukca_um_dealloc()

! Deallocate variables local to ukca_main
IF (ALLOCATED(theta_latitude)) DEALLOCATE(theta_latitude)
IF (ALLOCATED(v_latitude)) DEALLOCATE(v_latitude)
IF (ALLOCATED(sinv_latitude)) DEALLOCATE(sinv_latitude)
IF (ALLOCATED(um_ozone1D)) DEALLOCATE(um_ozone1D)
IF (ALLOCATED(surf_albedo)) DEALLOCATE(surf_albedo)
IF (ALLOCATED(cos_zenith_angle)) DEALLOCATE(cos_zenith_angle)
IF (ALLOCATED(T_chem)) DEALLOCATE(T_chem)
IF (ALLOCATED(Q_chem)) DEALLOCATE(Q_chem)
IF (ALLOCATED(fluxdiag_all_tracers)) DEALLOCATE(fluxdiag_all_tracers)
IF (ALLOCATED(um_ozone3d)) DEALLOCATE(um_ozone3d)
IF (ALLOCATED(ls_ppn3d)) DEALLOCATE(ls_ppn3d)
IF (ALLOCATED(conv_ppn3d)) DEALLOCATE(conv_ppn3d)
IF (ALLOCATED(dj)) DEALLOCATE(dj)
IF (ALLOCATED(rel_humid_frac_clr)) DEALLOCATE(rel_humid_frac_clr)
IF (ALLOCATED(Tile_Frac)) DEALLOCATE(Tile_Frac)
IF (ALLOCATED(Thick_bl_levels)) DEALLOCATE(Thick_bl_levels)
IF (ALLOCATED(STASHwork34)) DEALLOCATE(STASHwork34)
IF (ALLOCATED(STASHwork38)) DEALLOCATE(STASHwork38)
IF (ALLOCATED(STASHwork50)) DEALLOCATE(STASHwork50)
IF (ALLOCATED(p_layer_boundaries)) DEALLOCATE(p_layer_boundaries)
IF (ALLOCATED(t_theta_levels)) DEALLOCATE(t_theta_levels)
IF (ALLOCATED(delSO2_wet_h2o2)) DEALLOCATE(delSO2_wet_h2o2)
IF (ALLOCATED(delSO2_wet_o3)) DEALLOCATE(delSO2_wet_o3)
IF (ALLOCATED(delh2so4_chem)) DEALLOCATE(delh2so4_chem)
IF (ALLOCATED(delSO2_drydep)) DEALLOCATE(delSO2_drydep)
IF (ALLOCATED(delSO2_wetdep)) DEALLOCATE(delSO2_wetdep)
IF (ALLOCATED(so4_sa)) DEALLOCATE(so4_sa)
IF (ALLOCATED(all_tracers)) DEALLOCATE(all_tracers)

first=.FALSE.
IF ( do_chemistry ) firstchem = .FALSE.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

END SUBROUTINE ukca_main1

END MODULE ukca_main1_mod
