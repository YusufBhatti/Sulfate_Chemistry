! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine SETCONA (ENDGAME VERSION) ------------------------------
!
!   Purpose : to calculate additional constants derived from
!        information in the input dump and pass out into the model via
!        modules.
!
!   Method:
!      0. Initialise vertical coordinates from eta values and surface
!         orography.
!      0.1 Initialise data arrays held in secondary dump space.
!      1. Trigonometric functions:
!          Convert from degrees to radians.
!          Calculate trig functions for this grid.
!          Update halos for trig functions.
!      2. Set up semi-lagrangian advection options.
!         Initialise land/soil points in index arrays.
!      3. Determine model level boundaries for L/M/H cloud diagnostics.
!      4. Add locally derived parameters.
!
!   Language: FORTRAN 77 + common extensions also in Fortran 90.
!   Programming standard; Unified Model Documentation Paper No. 3
!
!   Documentation : Unified Model Documentation Paper No P0
!
!  -------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Top Level

SUBROUTINE Setcona_4A                                                   &
! Input data (including some logical fields)
           (eta_theta_in, eta_rho_in, z_theta_in, z_rho_in,             &
            Smvcst, Land_sea_mask, Orog,                                &
            grad_x, grad_y, orog_unfilt,                                &
            rho,exner_rho_levels,                                       &
            orog_lbc,exner_lbc,                                         &
! Input size and control variables
            lenrima,lbc_sizea,lbc_starta,                               &
            rimwidtha, rimweightsa,                                     &
            global_row_length, global_rows,                             &
            rows,n_rows,row_length,                                     &
            land_field,boundary_layer_levels,                           &
            first_constant_r_rho_level,cloud_levels,z_top_of_model,     &
            ozone_levels,bl_levels,tr_levels, tr_vars, tr_ukca,         &
            height_gen_method,                                          &
! other constants
            Delta_lambda_in,Delta_phi_in,Base_phi_in,Base_lambda_in,    &
            lat_rot_NP_in,long_rot_NP_in,                               &
! VarRes grid info in degrees
            Lambda_p_in, Lambda_u_in, Phi_p_in, Phi_v_in,               &
! FixHd values:
            CoordSystem, RowDepCStart,                                  &
! Initialise variables in secondary D1:
            exner_theta_levels,p_theta_levels,p,p_star,                 &
! Output data
           a_len2_coldepc, a_len2_rowdepc,                              &
           icode,Cmessage)

USE swap_bounds_mv_mod, ONLY: swap_bounds_mv
USE conversions_mod, ONLY: pi_over_180
USE rimtypes
USE solinc_data, ONLY: L_orog, slope_angle, horiz_limit, l_skyview,     &
    slope_aspect, horiz_ang, horiz_aspect
USE level_heights_Mod, ONLY: eta_theta_levels, eta_rho_levels,          &
                             z_ref_theta,      z_ref_rho,               &
                             r_theta_levels,   r_rho_levels,            &
                             r_at_u, r_at_v, r_at_u_w, r_at_v_w,        &
                             z_top_theta, d_layer,                      &
                             r_layer_centres, r_layer_boundaries
USE trignometric_Mod
USE diag_ctl_mod, ONLY: init_print_diag
USE atm_step_const_mod
USE atm_step_timestep_mod, ONLY: atm_step_time_init
USE dyn_coriolis_Mod
USE dyn_var_res_Mod
USE dynamics_testing_mod, ONLY: l_idealised_data, problem_number
USE problem_mod,          ONLY: standard
USE diff_coeff_Mod
USE turb_diff_mod
USE rad_mask_trop_Mod
USE rad_input_mod, ONLY:                                                &
    L_use_orog_corr, L_use_grad_corr,                                   &
    l_use_skyview,L_rad_deg,l_orog_unfilt, l_use_cariolle, l_extra_top
USE ancil_info, ONLY: l_soil_point, l_lice_point
USE rot_coeff_Mod
USE swapable_field_mod, ONLY: swapable_field_pointer_type
USE atm_fields_bounds_Mod
USE planet_constants_mod, ONLY: planet_radius, two_omega, &
                                recip_kappa, p_zero
USE dust_parameters_mod, ONLY: l_twobin_dust

USE global_2d_sums_mod, ONLY: global_sum_method, global_sum_reprod_dd,  &
                              global_sum_reprod_orig, global_sum_fast

USE bl_option_mod, ONLY: nl_bl_levels, z_nl_bl_levels,                  &
                         i_bl_vn, i_bl_vn_9b, i_bl_vn_9c, off, alpha_cd,&
                         alpha_cd_in

USE mym_option_mod, ONLY: bdy_tke, mymodel3
USE stochastic_physics_run_mod, ONLY: zmin_pert_theta, minlev_pert_theta,&
                                      zmax_pert_theta, maxlev_pert_theta,&
                                      i_pert_theta
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE Control_Max_Sizes
USE cderived_mod, ONLY: num_cloud_types, h_split, low_bot_level,  &
                        low_top_level, med_bot_level, med_top_level,  &
                        high_bot_level, high_top_level, delta_lambda, &
                        delta_phi, base_lambda, base_phi, lat_rot_np, &
                        long_rot_np
USE UM_ParVars
USE UM_ParCore,   ONLY: mype
USE UM_ParParams, ONLY: Ndim_max, NHalo_max,                   &
                        halo_type_extended, halo_type_no_halo, &
                        peast, pwest, pnorth, psouth
USE Field_Types,  ONLY: Nfld_max, fld_type_p, fld_type_u, fld_type_v
USE ereport_mod, ONLY: ereport
USE umPrintMgr
USE set_horiz_grid_mod
USE set_var_horiz_grid_mod
USE horiz_grid_mod
USE alloc_grid_mod
USE calc_global_grid_spacing_mod

USE halo_exchange, ONLY: swap_bounds
USE mpp_conf_mod,  ONLY: swap_field_is_scalar

USE lam_lbc_weights_mod, ONLY: init_lbc_wghts

USE trophgt1_mod, ONLY: z_min_trop, z_max_trop
USE sl_input_mod, ONLY:                                                 &
          thmono_levels, L_mono,L_high,high_order_scheme,thmono_height



USE dynamics_input_mod, ONLY:   L_LBC_balance,L_lbc_new,L_transparent,  &
                                GCR_diagnostics,GCR_its_avg_step,       &
                                L_lbc_old

USE lam_config_inputs_mod, ONLY: n_rims_to_do

USE eg_parameters_mod, ONLY: pole_consts
USE idealise_run_mod, ONLY: f_plane, ff_plane,                          &
                l_shallow, l_const_grav, zero_orography

USE gcr_input_mod, ONLY:                                                &
          GCR_its_switch, GCR_max_its,GCR_min_its,GCR_max_time,         &
          GCR_min_time,GCR_sum_its

USE mphys_inputs_mod, ONLY: l_mcr_qcf2, l_mcr_qgraup, l_mcr_qrain

USE cloud_inputs_mod, ONLY: i_cld_vn
USE pc2_constants_mod, ONLY: i_cld_pc2
USE cv_run_mod, ONLY: l_conv_prog_group_1, l_conv_prog_group_2,         &
                      l_conv_prog_group_3, l_conv_prog_precip

USE murk_inputs_mod, ONLY: l_murk_advect

USE dust_parameters_mod, ONLY: l_dust

USE nlsizes_namelist_mod,  ONLY: super_array_size, moisture_array_size, &
     turb_array_size
USE carbon_options_mod, ONLY: l_co2_interactive
USE lbc_read_data_mod, ONLY: l_int_uvw_lbc
USE run_aerosol_mod, ONLY:   &
     l_sulpc_so2, l_sulpc_dms, l_sulpc_nh3, l_soot, l_ocff, l_biomass,  &
     l_nitrate
USE ac_control_mod, ONLY: long_e, long_w, lat_n, lat_s,                 &
                          long_w_model, long_e_model

USE missing_data_mod, ONLY: imdi, rmdi
USE highos_mod,       ONLY: quinticLagrange

USE coriolis_mod
USE gravity_mod
USE helmholtz_const_matrix_mod
USE eg_helmholtz_mod
USE departure_pts_mod
USE ref_pro_mod
USE eg_star_mod
USE eg_lam_domain_kind_mod
USE eg_sisl_setcon_mod

USE errormessagelength_mod, ONLY: errormessagelength

USE model_domain_mod, ONLY: l_regular, model_type, mt_global, mt_lam, &
                            mt_bi_cyclic_lam
USE dump_headers_mod, ONLY: fh_CoordSystem_Cartesian

USE nlsizes_namelist_mod, ONLY: model_levels

USE land_soil_dimensions_mod, ONLY: land_points, land_ice_points, &
     soil_points, land_index, land_ice_index, soil_index

USE aspang_mod, ONLY: aspang
USE aspang_ancil_mod, ONLY: aspang_ancil
USE coeffs_degrade_mod, ONLY: coeffs_degrade
USE pole_bearing_mod, ONLY: pole_bearing
USE rad_degrade_mask_mod, ONLY: rad_degrade_mask
USE skyview_mod, ONLY: skyview
IMPLICIT NONE

INTEGER ::                                                        &
       lenrima(Nfld_max,NHalo_max,Nrima_max),                     &
                                               ! IN LBC data len
       lbc_sizea(4,Nfld_max,NHalo_max,Nrima_max),                 &
                                                 ! IN LBC size
       lbc_starta(4,Nfld_max,NHalo_max,Nrima_max),                &
                                                   ! IN LBC start
       rimwidtha(Nrima_max),                                      &
                              ! IN RIM width
       global_row_length,                                         &
                           ! IN total number of point in a row
       global_rows,                                               &
                     ! IN total number of rows in model
       rows,n_rows,row_length,icode,                              &
       land_field,                                                &
                     ! IN Number of land points in model
       first_constant_r_rho_level,                                &
                                   !(IN) 1st constant height level
       boundary_layer_levels,                                     &
                              ! (IN) Num. of boundary layer levels
       height_gen_method,                                         &
                            ! (IN) Method for generating heights
       cloud_levels                                               &
                     ! IN No of cloudy levels in the model
      ,tr_levels, tr_vars, tr_ukca

! FixHd values:
INTEGER, INTENT(IN) :: CoordSystem    ! Chosen grid coordinate system
INTEGER :: RowDepCStart
                  ! IN Start of Row dependent constant
CHARACTER(LEN=errormessagelength) ::                              &
       cmessage              ! Error message if ICODE >0
CHARACTER(LEN=*) :: RoutineName
PARAMETER (   RoutineName='SETCONA_4A')

! The following are used in the decalaration of Lambda_u_in, Phi_v_in
!  to avoid out-of-bounds errors
INTEGER :: a_len2_coldepc, a_len2_rowdepc

REAL ::                                                           &
      ! Input arguments (grid constants)
       Delta_lambda_in                                            &
                       ! EW (x) grid spacing in degrees
      ,Delta_phi_in                                               &
                       ! NS (y) grid spacing in degrees
      ,Base_phi_in                                                &
                       ! Latitude of first theta point in degrees
      ,Base_lambda_in                                             &
                       ! Longitude of first theta point in degs
      ,lat_rot_NP_in                                              &
                     ! Real latitude of 'pseudo' N pole in degs
      ,long_rot_NP_in                                             &
                     ! Real longitude of 'pseudo' N pole in degs
      ,z_top_of_model                                             &
                      ! (IN) Height of top of model in metres
,      rimweightsa(rimwidtha(rima_type_norm))  ! IN RIM weights

! In the regular grid case A_LEN2_COLDEPC, A_LEN2_ROWDEPC are zero, but
! Lambda_u_in, Phi_v_in must have size 1 to avoid out-of-bounds errors
REAL ::                                                           &
      ! Input VarRes grid info in degrees
   Lambda_p_in(global_row_length)                                 &
                                         !IN EW and NS VarRes grid
  ,Lambda_u_in(MIN(global_row_length, &
                   global_row_length *a_len2_coldepc +1))         &
                                         !IN EW and NS u,v grid
  ,Phi_p_in(global_rows)                                          &
                                         !location in degrees
  ,Phi_v_in(MIN(global_rows+1,                                    &
                global_rows *a_len2_rowdepc +1)) !location in degrees

REAL ::                                                           &
      ! Array arguments with intent(IN):
     Smvcst(land_field)                                           &
                                ! IN Volumetric saturation point
    ,eta_theta_in(0:model_levels)                                 &
                                  ! IN Eta values for theta levs
    ,eta_rho_in(model_levels)                                     &
                                  ! IN Eta values for rho levels
    ,z_theta_in(0:model_levels)                                   &
                                  ! IN height above sea-level
    ,z_rho_in(model_levels)                                       &
                                  ! IN height above sea-level
    ,Orog(row_length, rows)                                       &
                                  ! IN Orography (on all points)
    ,orog_unfilt(row_length, rows)                                &
                                  ! IN Unfiltered orography
    ,grad_x(row_length*rows)                                      &
                                  ! IN Orographic X-gradient
    ,grad_y(row_length*rows)                                      &
                                  ! IN Orographic Y-gradient
    ,rho(pdims_s%i_start:pdims_s%i_end,                           &
         pdims_s%j_start:pdims_s%j_end,                           &
         pdims_s%k_start:pdims_s%k_end)                           &
                     ! density*(r**2): used in call to Calc_P_star
    ,exner_rho_levels(pdims_s%i_start:pdims_s%i_end,              &
                      pdims_s%j_start:pdims_s%j_end,              &
                      pdims_s%k_start:pdims_s%k_end+1)            &
                                               ! Exner at p levels
,      orog_lbc(lenrima(fld_type_p,halo_type_extended,            &
                        rima_type_orog))                          &
                                           ! Orography LBC
,      exner_lbc(lenrima(fld_type_p,halo_type_extended,           &
                                                        !Exner LBC
                         rima_type_norm),model_levels+1)

REAL ::                                                           &
      ! Output args (secondary arrays calculated from dump)
  p_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                 tdims_s%j_start:tdims_s%j_end,                   &
                 tdims_s%k_start:tdims_s%k_end)                   &
                                         ! press on theta levs Pa
, p(pdims_s%i_start:pdims_s%i_end,                                &
    pdims_s%j_start:pdims_s%j_end,                                &
    pdims_s%k_start:pdims_s%k_end+1)                              &
                                                        ! in Pa
, p_star(row_length, rows)                                        &
                           ! surface pressure (Pa)
, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,               &
                     tdims_s%j_start:tdims_s%j_end,               &
                     tdims_s%k_start:tdims_s%k_end)
                                             ! Exner at theta levs

LOGICAL ::                                                        &
       Land_sea_mask(row_length,rows)


! Local variables

REAL :: orog_halos(1-halo_i:row_length+halo_i,                    &
                   1-halo_j:rows+halo_j)

REAL :: orog_horiz(1-horiz_limit:row_length+horiz_limit,          &
                   1-horiz_limit:rows+horiz_limit)

REAL :: orog_global(global_row_length,global_rows)
REAL :: max_orog  ! Maximum orographic height

REAL :: phi_horiz(1-horiz_limit:rows+horiz_limit)
REAL :: dphi_horiz(1-horiz_limit:rows+horiz_limit)
REAL :: dlambda_horiz(1-horiz_limit:row_length+horiz_limit)

REAL ::                                                           &
   Delta_lambda_wk                                                &
                       ! EW (x) grid spacing in degrees
  ,Delta_phi_wk                                                   &
                       ! NS (y) grid spacing in degrees
  ,Base_phi_wk                                                    &
                       ! Latitude of first theta point in degrees
  ,Base_lambda_wk                                                 &
                       ! Longitude of first theta point in degrees
  ,latitude_rot(row_length, rows)                                 &
                                   ! rot. latit. in degrees (LAM)
  ,longitude_rot(row_length, rows)                                &
                                   ! rot. longit. in degrees (LAM)
  ,true_latitude_b(row_length, rows)                              &
                                     ! true latitude on B-grid
  ,true_longitude_b(row_length,rows)                              &
                                     ! true longitude on B-grid
  ,bear_rot_NP(row_length, rows)                                  &
                                   ! Bearing of 'pseudo' N pole
  ,temp1                                                          &
             ! Used in calculation of Coriolis terms (LAM only)
  ,temp2                                                          &
             ! Used in calculation of Coriolis terms (LAM only)
   , f_plane_rad                                                  &
                  ! f_plane latitude in radians
   , ff_plane_rad                                                 &
                   ! f_plane latitude in radians
   , f1_temp                                                      &
              ! f1 term for f-plane or ff-plane
   , f2_temp                                                      &
              ! f2 term for f-plane or ff-plane
   , f3_temp                                                      &
              ! f3 term for f-plane or ff-plane
  , SCALE                                                         &
  ,cloud_bound(num_cloud_types+1)                                 &
                                  ! boundaries of cloud types
  ,r_ref_theta(model_levels)                                      &
                             ! Local dynamic array
  ,r_ref_rho(model_levels)
                             ! Local dynamic array

LOGICAL ::                                                        &
   landmask(1-offx:row_length+offx, 1-offy:rows+offy)             &
  ,L_error                                                        &
  ,L_varres

LOGICAL ::                                                        &
  l_do_boundaries                                                 &
, l_do_halos                                                      &
, l_include_halos

INTEGER ::                                                        &
        ! Mostly loop counters, but j also used for interp points
   i,kk,                                                          &
   j,gi,gj,                                                       &
   j0,j1,k,                                                       &
   level                                                          &
                   ! Used to set up cloud type boundaries
  ,first_row                                                      &
  ,last_row                                                       &
  ,info                                                           &
  , active_levels

INTEGER, PARAMETER :: height_gen_original = 1 ! methods for height
INTEGER, PARAMETER :: height_gen_smooth   = 2 ! generation
INTEGER, PARAMETER :: height_gen_linear   = 3

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER :: i_field

TYPE(swapable_field_pointer_type) :: fields_to_swap(2)

INTEGER :: ozone_levels,bl_levels
REAL    :: z_offset


!----------------------------------------------------------------------

IF (lhook) CALL dr_hook(RoutineName,zhook_in,zhook_handle)

! Allocate workspace arrays for ENDGame model evolution

CALL init_coriolis()
CALL init_ref_pro()
CALL init_depart_pts()
CALL init_helmholtz_const_matrix()
CALL init_gravity()
CALL init_star()
CALL eg_init_helmholtz(row_length,rows,n_rows,model_levels,offx,offy)


! Set flag for orography correction in module solinc_data
L_orog = L_use_orog_corr .OR. L_use_grad_corr
l_skyview = l_orog .AND. l_use_skyview
! Grid type in dump: L_varres = T / F for variable / unif grid
L_varres = .FALSE.
IF (RowDepCStart >= 1) L_varres = .TRUE.
! ----------------------------------------------------------------------
! 0. Set-up vertical co-ordinate arrays
! ----------------------------------------------------------------------

z_offset = planet_radius

DO j=1,rows
  DO i=1,row_length
    orog_halos(i,j)=orog(i,j)
  END DO
END DO

CALL swap_bounds(orog_halos,row_length,rows,1,                    &
                 halo_i,halo_j,fld_type_p, swap_field_is_scalar)

! ----------------------------------------------------------------------
! 0.1 Write out some useful information at start of run
! ----------------------------------------------------------------------

CALL umPrint('',src='setcon_4A')
IF (global_sum_method == global_sum_reprod_orig .OR.              &
    global_sum_method == global_sum_reprod_dd) THEN
  CALL umPrint('**********************  This run uses bit-'//     &
      'reproducible code ********************',src='setcon_4A')
ELSE IF (global_sum_method == global_sum_fast) THEN
  CALL umPrint('***************** This run uses fast code, so'//  &
      ' different processor **************',src='setcon_4A')
  CALL umPrint('*****************     configurations may not '//  &
      'bit-reproduce        **************',src='setcon_4A')
END IF

CALL umPrint('',src='setcona_4A')
CALL umPrint('************************  PE configuration for '//         &
    'this run  ***********************',src='setcon_4A')
IF ( nproc_x > 1 .AND. nproc_y > 1 ) THEN
  WRITE(umMessage,'(A,I5,A,I5,A)') '****  PE configuration: ', nproc_x,  &
      ' processors East-West ', nproc_y, ' processors North-South *'
  CALL umPrint(umMessage,src='setcona_4A')
ELSE IF ( nproc_x > 1 .AND. nproc_y == 1 ) THEN
  WRITE(umMessage,'(A,I5,A,I5,A)') '****  PE configuration: ', nproc_x,  &
      ' processors East-West ', nproc_y, ' processor North-South *'
  CALL umPrint(umMessage,src='setcona_4A')
ELSE IF ( nproc_x == 1 .AND. nproc_y > 1 ) THEN
  WRITE(umMessage,'(A,I5,A,I5,A)') '****  PE configuration:', nproc_x,   &
      ' processor East-West ', nproc_y, ' processors North-South *'
  CALL umPrint(umMessage,src='setcona_4A')
ELSE ! nproc_x = 1 .and. nproc_y =1 ) then
  CALL umPrint('****  PE configuration:  Single processor 1x1 * ',       &
      src='setcona_4A')
END IF ! nproc_x > 1 .and. nproc_y > 1 )

IF ( offx >= halo_i .OR. offy >= halo_j ) THEN
  CALL umPrint('',src='setcona_4A')
  CALL umPrint( '  ****  DANGER  ****   ',src='setcona_4A')
  CALL umPrint( '  *** SMALL halo >= LARGE halo  *****   ',      &
      src='setcona_4A')
  CALL umPrint( '  This could result in overwriting  ',src='setcona_4A')
  icode    = 20
  cmessage ='SETCONA: Large halo size is too small'
  CALL Ereport(RoutineName,icode,Cmessage)
END IF ! offx >= halo_i .or. offy >= halo_j

CALL umPrint('',src='setcona_4A')

IF (model_type == mt_lam) THEN

  L_lbc_old = .TRUE.
  IF ( L_lbc_new ) L_lbc_old = .FALSE.

  CALL umPrint('',src='setcona_4A')
  CALL umPrint(' ***   LAM set-up ***',src='setcona_4A')

  WRITE(umMessage,*)'LBC frame size for solver, ',          &
      ' n_rims_to_do = ', n_rims_to_do
  CALL umPrint(umMessage,src='setcona_4A')
  IF ( L_LBC_balance ) THEN
    WRITE(umMessage,*)'L_LBC_balance = ', L_LBC_balance,    &
        ' Impose vertically balanced Exner pressures',      &
        ' and rho and set w=0 in lbcs '
    CALL umPrint(umMessage,src='setcona_4A')
  ELSE !  L_LBC_balance = .false.
    WRITE(umMessage,*)'L_LBC_balance = ', L_LBC_balance,    &
        ' No balancing of lbcs '
    CALL umPrint(umMessage,src='setcona_4A')
  END IF ! L_LBC_balance
  IF ( L_lbc_new ) THEN
    WRITE(umMessage,*) 'L_lbc_new = ', L_lbc_new            &
,                              'Use new lbc algorithm '
    CALL umPrint(umMessage,src='setcona_4A')
    IF ( L_transparent ) THEN
      WRITE(umMessage,*) 'L_transparent = ', L_transparent  &
,         'Does not apply lbcs to horizontal winds in the solver'
      CALL umPrint(umMessage,src='setcona_4A')
    ELSE !  L_transparent = .false.
      WRITE(umMessage,*) 'L_transparent = ', L_transparent  &
,         'Applies lbcs to horizontal winds in the solver'
      CALL umPrint(umMessage,src='setcona_4A')
    END IF ! L_transparent
  ELSE !  L_lbc_new = .false.
    WRITE(umMessage,*) 'L_lbc_new = ', L_lbc_new            &
,                              'Use old lbc algorithm '
    CALL umPrint(umMessage,src='setcona_4A')
  END IF ! L_lbc_new
  IF ( L_int_uvw_lbc) THEN
    WRITE(umMessage,*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
,          ' Advecting winds interpolated in lateral boundaries '
    CALL umPrint(umMessage,src='setcona_4A')
  ELSE !  L_int_uvw_lbc = .false.
    WRITE(umMessage,*) 'L_int_uvw_lbc =', L_int_uvw_lbc     &
,                ' Extrapolated advecting winds from lbc file '
    CALL umPrint(umMessage,src='setcona_4A')
  END IF ! L_int_uvw_lbc
  CALL umPrint('',src='setcona_4A')

  l_do_boundaries = .FALSE.
  l_do_halos = .TRUE.

  IF ( rimwidtha(rima_type_orog) > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(                                  &
      row_length,rows,halo_i,halo_j,1,fld_type_p,orog_halos,      &
      lenrima(fld_type_p,halo_type_extended,rima_type_orog),      &
      lbc_sizea(1,fld_type_p,halo_type_extended,rima_type_orog),  &
      lbc_starta(1,fld_type_p,halo_type_extended,rima_type_orog), &
      halo_i,halo_j,orog_lbc,rimwidtha(rima_type_orog),           &
      rimwidtha(rima_type_orog),rimweightsa,                      &
      at_extremity,                                               &
      l_do_boundaries,l_do_halos)
  END IF

  IF ( rimwidtha(rima_type_norm) > 0 ) THEN
  ! DEPENDS ON: set_lateral_boundaries
    CALL set_lateral_boundaries(                                  &
      row_length,rows,Offx,Offy,model_levels+1,fld_type_p,        &
      exner_rho_levels,                                           &
      lenrima(fld_type_p,halo_type_extended,rima_type_norm),      &
      lbc_sizea(1,fld_type_p,halo_type_extended,rima_type_norm),  &
      lbc_starta(1,fld_type_p,halo_type_extended,rima_type_norm), &
      halo_i,halo_j,exner_lbc,rimwidtha(rima_type_norm),          &
      rimwidtha(rima_type_norm),rimweightsa,                      &
      at_extremity,                                               &
      l_do_boundaries,l_do_halos)
  END IF

  CALL init_lbc_wghts(row_length,rows,rimwidtha(rima_type_norm),  &
                   halo_i, halo_j, at_extremity, rimweightsa )

END IF ! IF (model_type == mt_lam)

! Set reference height profile

! Allocate arrays for level_heights_mod module
IF (.NOT. ALLOCATED(eta_theta_levels)) THEN
  ALLOCATE (eta_theta_levels(0:model_levels))
END IF
IF (.NOT. ALLOCATED(eta_rho_levels)) THEN
  ALLOCATE (eta_rho_levels(model_levels))
END IF
IF (.NOT. ALLOCATED(r_theta_levels)) THEN
  ALLOCATE (r_theta_levels(1-halo_i:row_length+halo_i,            &
                           1-halo_j:rows+halo_j, 0:model_levels))
END IF
IF (.NOT. ALLOCATED(r_rho_levels)) THEN
  ALLOCATE (r_rho_levels(1-halo_i:row_length+halo_i,              &
                         1-halo_j:rows+halo_j, model_levels))
END IF

eta_theta_levels(0) = eta_theta_in(0)
z_top_theta = z_top_of_model
DO k = 1, model_levels
  eta_theta_levels(k) = eta_theta_in(k)
  eta_rho_levels(k) = eta_rho_in(k)
  r_ref_theta(k) = eta_theta_levels(k) * z_top_of_model
END DO
DO k = 1, model_levels
  r_ref_rho(k) = eta_rho_levels(k) * z_top_of_model
END DO
! set bottom level, ie: orography
DO j = 1-halo_j, rows+halo_j
  DO i= 1-halo_i, row_length+halo_i
    r_theta_levels(i,j,0) = Orog_halos(i,j) + z_offset
  END DO
END DO
! For constant levels set r to be a constant on the level
DO k = first_constant_r_rho_level, model_levels
  DO j = 1-halo_j, rows+halo_j
    DO i= 1-halo_i, row_length+halo_i
      r_theta_levels(i,j,k) = z_offset + r_ref_theta(k)
      r_rho_levels(i,j,k)   = z_offset + r_ref_rho(k)
    END DO
  END DO
END DO

SELECT CASE( height_gen_method )
CASE ( height_gen_original )
  ! The original version of height generation used in the SI dynamics
  !
  ! For boundary layer levels set depth to be constant.
  DO k = 1, boundary_layer_levels
    DO j = 1-halo_j, rows+halo_j
      DO i= 1-halo_i, row_length+halo_i
        r_theta_levels(i,j,k) = r_theta_levels(i,j,0)               &
                                   + r_ref_theta(k)
        r_rho_levels(i,j,k) = r_theta_levels(i,j,0)                 &
                                 + r_ref_rho(k)
      END DO
    END DO
  END DO
  ! For intermediate levels use linear relaxation to constant value.
  ! set orographic heights.
  DO k = boundary_layer_levels+1, first_constant_r_rho_level-1
    DO j = 1-halo_j, rows+halo_j
      DO i= 1-halo_i, row_length+halo_i
        r_rho_levels(i,j,k) =                                       &
          ( r_rho_levels(i,j,first_constant_r_rho_level) -          &
            r_theta_levels(i,j,boundary_layer_levels) ) *           &
          ( eta_rho_levels(k) -                                     &
            eta_theta_levels(boundary_layer_levels) ) /             &
          ( eta_rho_levels(first_constant_r_rho_level) -            &
            eta_theta_levels(boundary_layer_levels) )               &
           +  r_theta_levels(i,j,boundary_layer_levels)
        r_theta_levels(i,j,k) =                                     &
            ( r_rho_levels(i,j,first_constant_r_rho_level) -        &
              r_theta_levels(i,j,boundary_layer_levels) ) *         &
            ( eta_theta_levels(k) -                                 &
              eta_theta_levels(boundary_layer_levels) ) /           &
            ( eta_rho_levels(first_constant_r_rho_level) -          &
              eta_theta_levels(boundary_layer_levels) )             &
             +  r_theta_levels(i,j,boundary_layer_levels)
      END DO
    END DO
  END DO

CASE ( height_gen_smooth )
  ! A smooth quadratic height generation
  DO k = 1, first_constant_r_rho_level-1
    DO j = 1-halo_j, rows+halo_j
      DO i= 1-halo_i, row_length+halo_i

        r_rho_levels(i,j,k) = eta_rho_levels(k) * z_top_of_model +&
         z_offset + Orog_halos(i,j) * (1.0 - eta_rho_levels(k)    &
              /eta_rho_levels(first_constant_r_rho_level))**2

        r_theta_levels(i,j,k) = eta_theta_levels(k) *             &
             z_top_of_model + z_offset + Orog_halos(i,j) *        &
             (1.0 - eta_theta_levels(k) /                         &
              eta_rho_levels(first_constant_r_rho_level))**2

      END DO
    END DO
  END DO

CASE (height_gen_linear)
  ! A linear progression of heights
  DO k = 1, first_constant_r_rho_level-1
    DO j = 1-halo_j, rows+halo_j
      DO i= 1-halo_i, row_length+halo_i
        r_theta_levels(i,j,k) = z_offset + r_ref_theta(k)         &
                     +  Orog_halos(i,j) * (1.0 - eta_theta_levels(k))          
        r_rho_levels(i,j,k)   = z_offset + r_ref_rho(k)           &
                     +  Orog_halos(i,j) * (1.0 - eta_rho_levels(k))
      END DO
    END DO
  END DO

CASE DEFAULT
  icode = 10
  WRITE (Cmessage,*) 'Unrecognised height generation method - ',  &
                     'Dump needs to be reconfigured'
  CALL Ereport( RoutineName, icode, Cmessage )
END SELECT

! ENDGAME-only

CALL swap_bounds(r_theta_levels,                                       &
                 tdims_l%i_len - 2*tdims_l%halo_i,                     &
                 tdims_l%j_len - 2*tdims_l%halo_j,                     &
                 tdims_l%k_len,                                        &
                 tdims_l%halo_i, tdims_l%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)
CALL swap_bounds(r_rho_levels,                                         &
                 pdims_l%i_len - 2*pdims_l%halo_i,                     &
                 pdims_l%j_len - 2*pdims_l%halo_j,                     &
                 pdims_l%k_len,                                        &
                 pdims_l%halo_i, pdims_l%halo_j,                       &
                 fld_type_p, swap_field_is_scalar)

! calculate r_at_u points on rho levels and
! r_at_v points on rho levels.
!

! Allocate r_at_u/v arrays for level_heights_mod module
IF (.NOT. ALLOCATED(r_at_u)) THEN
  ALLOCATE (r_at_u (udims_l%i_start:udims_l%i_end,                &
                    udims_l%j_start:udims_l%j_end,                &
                    udims_l%k_start:udims_l%k_end))
END IF
IF (.NOT. ALLOCATED(r_at_v)) THEN
  ALLOCATE (r_at_v (vdims_l%i_start:vdims_l%i_end,                &
                    vdims_l%j_start:vdims_l%j_end,                &
                    vdims_l%k_start:vdims_l%k_end))
END IF

IF (.NOT. ALLOCATED(r_at_u_w)) THEN
  ALLOCATE (r_at_u_w(udims_l%i_start:udims_l%i_end,               &
                     udims_l%j_start:udims_l%j_end,               &
                           0:model_levels))
END IF

IF (.NOT. ALLOCATED(r_at_v_w)) THEN
  ALLOCATE (r_at_v_w(vdims_l%i_start:vdims_l%i_end,               &
                     vdims_l%j_start:vdims_l%j_end,               &
                            0:model_levels))
END IF

! NOTE: Simple regular grid is assumed
!
DO k=1, model_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      r_at_u(i,j,k) = .5 * (r_rho_levels(i,j,k) +                 &
                            r_rho_levels(i+1,j,k) )
    END DO
  END DO
END DO
DO k = 0, model_levels
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      r_at_u_w(i,j,k) = .5 * (r_theta_levels(i,j,k) +             &
                              r_theta_levels(i+1,j,k) )
    END DO
  END DO
END DO

DO k = 1, model_levels
  DO j = vdims%j_start,vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      r_at_v(i,j,k) = .5 * (r_rho_levels(i,j,k) +                 &
                            r_rho_levels(i,j+1,k) )
    END DO
  END DO
END DO
DO k = 0, model_levels
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      r_at_v_w(i,j,k) = .5 * (r_theta_levels(i,j,k) +             &
                              r_theta_levels(i,j+1,k) )
    END DO
  END DO
END DO

CALL swap_bounds(r_at_u,                                               &
                 udims_l%i_len - 2*udims_l%halo_i,                     &
                 udims_l%j_len - 2*udims_l%halo_j,                     &
                 udims_l%k_len,                                        &
                 udims_l%halo_i, udims_l%halo_j,                       &
                 fld_type_u,swap_field_is_scalar)
CALL swap_bounds(r_at_v,                                               &
                 vdims_l%i_len - 2*vdims_l%halo_i,                     &
                 vdims_l%j_len - 2*vdims_l%halo_j,                     &
                 vdims_l%k_len,                                        &
                 vdims_l%halo_i, vdims_l%halo_j,                       &
                 fld_type_v,swap_field_is_scalar)

CALL swap_bounds(r_at_u_w,row_length,rows,                             &
                 model_levels+1,halo_i,halo_j,                         &
                 fld_type_u,swap_field_is_scalar)
CALL swap_bounds(r_at_v_w,row_length,n_rows,                           &
                 model_levels+1,halo_i,halo_j,                         &
                 fld_type_v,swap_field_is_scalar)
CALL swap_bounds(exner_rho_levels,row_length,rows,                     &
                 model_levels+1,offx,offy,                             &
                 fld_type_p, swap_field_is_scalar)

! End of ENDGAME-only


! ----------------------------------------------------------------------
! Set level heights for physics routines
! ----------------------------------------------------------------------
IF (.NOT. ALLOCATED(r_layer_centres)) THEN
  ALLOCATE ( r_layer_centres(row_length, rows, model_levels+1) )
END IF
IF (.NOT. ALLOCATED(r_layer_boundaries)) THEN
  ALLOCATE ( r_layer_boundaries(row_length, rows, 0:model_levels+1) )
END IF
IF (.NOT. ALLOCATED(d_layer)) THEN
  ALLOCATE ( d_layer(row_length, rows, model_levels+1) )
END IF

! Bottom boundary is the surface
DO j = 1, rows
  DO i = 1, row_length
    r_layer_boundaries(i,j,0) = r_theta_levels(i,j,0)
  END DO
END DO
! The first rho level is not used, otherwise boundaries lie on rho levels
! and centres lie on theta levels
DO k = 1, model_levels - 1
  DO j = 1, rows
    DO i = 1, row_length
      r_layer_centres(i,j,k) = r_theta_levels(i,j,k)
      r_layer_boundaries(i,j,k) = r_rho_levels(i,j,k+1)
      d_layer(i,j,k) = r_layer_boundaries(i,j,k) - r_layer_boundaries(i,j,k-1)
    END DO
  END DO
END DO
! The top boundary is extrapolated using twice the distance to the layer centre
k = model_levels
DO j = 1, rows
  DO i = 1, row_length
    r_layer_centres(i,j,k) = r_theta_levels(i,j,k)
    r_layer_boundaries(i,j,k) = 2.0*r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
    d_layer(i,j,k) = r_layer_boundaries(i,j,k) - r_layer_boundaries(i,j,k-1)
  END DO
END DO
! A layer is specified above the top of the model which may be used by the
! radiation scheme to contain the remainder of the atmosphere. This is taken
! to be the same height as the layer below.
k = model_levels + 1
DO j = 1, rows
  DO i = 1, row_length
    d_layer(i,j,k) = d_layer(i,j,k-1)
    r_layer_centres(i,j,k) = r_layer_centres(i,j,k-1) + d_layer(i,j,k)
    r_layer_boundaries(i,j,k) = r_layer_boundaries(i,j,k-1) + d_layer(i,j,k)
  END DO
END DO
! ----------------------------------------------------------------------


IF ( PrintStatus == PrStatus_Diag .AND. mype == 0 ) THEN
  CALL umPrint('================================================', &
      src='setcona_4A')
  CALL umPrint('level  | r_theta_levls,r_rho_levls,etc @ i=1,j=1', &
      src='setcona_4A')
  CALL umPrint('------------------------------------------------', &
      src='setcona_4A')

  WRITE(umMessage,'(I4,A,E15.5)') 0,'   |', r_theta_levels(1,1,0)-planet_radius
  CALL umPrint(umMessage,src='setcona_4A')

  DO k=1,model_levels
    WRITE(umMessage,'(I4,A,6E15.5)') k,'   |',                         &
                              r_theta_levels(1,1,k)-planet_radius,     &
                              r_rho_levels  (1,1,k)-planet_radius,     &
                              r_at_u        (1,1,k)-planet_radius,     &
                              r_at_v        (1,1,k)-planet_radius,     &
                              r_at_u_w        (1,1,k)-planet_radius,   &
                              r_at_v_w        (1,1,k)-planet_radius
    CALL umPrint(umMessage,src='setcona_4A')
  END DO
END IF
IF ( mype == 0 ) THEN
  CALL umPrint('================================================', &
      src='setcona_4A')
  IF (l_extra_top) THEN
    WRITE(umMessage,'(A,F20.10)') &
      'Top-of-atmosphere: ', r_layer_boundaries(1,1,model_levels+1)
  ELSE
    WRITE(umMessage,'(A,F20.10)') &
      'Top-of-atmosphere: ', r_layer_boundaries(1,1,model_levels)
  END IF
  CALL umPrint(umMessage,src='setcona_4A')
  CALL umPrint('================================================', &
      src='setcona_4A')
END IF

! 0.1 Initialise secondary arrays.
! Exner at p (=rho) levels is obtained from dump.
! calculate p from exner_rho_levels
! [halos required for diagnostic calculations at T+0 in INITDIAG.]

DO k = 1, model_levels+1
  DO j = 1-offy, rows+offy
    DO i = 1-offx, row_length+offx
      p(i,j,k)= (exner_rho_levels(i,j,k) ** recip_kappa)          &
                      * p_zero
    END DO
  END DO
END DO

! calculate exner at theta_levels which is then used to get
! p at theta_levs
l_include_halos = .TRUE.
! DEPENDS ON: calc_exner_at_theta
CALL Calc_Exner_at_theta( r_theta_levels, r_rho_levels,           &
                 exner_rho_levels,                                &
                 row_length, rows, model_levels,                  &
                 offx, offy, halo_i, halo_j,                      &
                 exner_theta_levels,l_include_halos)

! Calculate pressure from Exner at theta levels.

! DEPENDS ON: calc_p_from_exner
CALL Calc_P_from_Exner(                                           &
                      p_theta_levels,                             &
                      row_length, rows,                           &
                      tdims_s%k_len,                              &
                      offx, offy,                                 &
                      exner_theta_levels,l_include_halos)

! calculate p_star using rho (from dump) and p on model levels
! adding lines seen in setcona (ND) so that pstar is initialsed to
! something in case it is requirested at T+0 in a downscaler without IAU.
! DEPENDS ON: calc_p_star
CALL Calc_P_star (r_theta_levels, r_rho_levels, p, rho,           &
                  row_length, rows, model_levels,                 &
                  offx, offy, halo_i, halo_j,                     &
                  p_star)

! ----------------------------------------------------------------------
! Initialise Constants and trigonometric functions.
! ----------------------------------------------------------------------

! Determine coordinate system in use.
! Long-Lat is assumed if not explicitly set as Cartesian.
cartesian_grid = (CoordSystem == fh_CoordSystem_Cartesian)

! This test cannot go in the relevant check_nml routine without RECON ifdefs.
IF (model_type == mt_bi_cyclic_lam .AND. .NOT. cartesian_grid) THEN
  icode = 20
  WRITE(cmessage,'(A,I0)')                                            &
  "Bicyclic LAM model domains are only supported on Cartesian grids." &
   //newline//                                                        &
  "The provided dump's co-ordinate system (FLH(11)) is: ",CoordSystem
  CALL ereport(RoutineName, icode, cmessage)
END IF

! Convert grid-spacing and base values from degrees to radians
IF ( l_regular ) THEN
  base_lambda_wk  = base_lambda_in
  base_phi_wk     = base_phi_in
  delta_lambda_wk = delta_lambda_in
  delta_phi_wk    = delta_phi_in
ELSE
  base_lambda_wk  = lambda_p_in(1)
  base_phi_wk     = phi_p_in(1)
  delta_lambda_wk = (lambda_p_in(global_row_length) -             &
                    lambda_p_in(1))/(global_row_length - 1)
  delta_phi_wk    = (phi_p_in(global_rows) -                      &
                    phi_p_in(1))/(global_rows - 1)
END IF

IF (cartesian_grid) THEN
  ! Physics makes use of angular mesh parameters to compute cell area
  delta_lambda    = delta_lambda_wk * pi_over_180 / planet_radius
  delta_phi       = delta_phi_wk * pi_over_180 / planet_radius
ELSE
  delta_lambda = delta_lambda_wk * pi_over_180
  delta_phi    = delta_phi_wk    * pi_over_180
END IF

base_phi     = base_phi_wk     * pi_over_180
base_lambda  = base_lambda_wk  * pi_over_180
f_plane_rad  = f_plane         * pi_over_180
ff_plane_rad = ff_plane        * pi_over_180
lat_rot_NP   = lat_rot_NP_in   * pi_over_180
long_rot_NP  = long_rot_NP_in  * pi_over_180

! For LAMs, the (xi1,xi2) mesh parameters are read from the dump:
!   - For Long-Lat grids these are in degrees, and the angular measures are
!     converted to radians in set_horiz_grid.
!   - For Cartesian grids these are in metres (and no conversion is needed).
! For global models the mesh parameters are defined in set_horiz_grid.
! Longitudinal mesh parameters for Long-Lat mono-cyclic LAMs are also defined
!   in set_horiz_grid.
base_xi1  = base_lambda_in
base_xi2  = base_phi_in
delta_xi1 = delta_lambda_in
delta_xi2 = delta_phi_in

CALL alloc_grid()

! Allocate arrays for dyn_var_res_mod module
! Allocate full size arrays to avoid out-of-bounds in atm_step
IF (.NOT. ALLOCATED(glambda_p)) THEN
  ALLOCATE (glambda_p ( 1-halo_i : global_row_length+halo_i ))
END IF
IF (.NOT. ALLOCATED(glambda_u)) THEN
  ALLOCATE (glambda_u ( 1-halo_i : global_row_length+halo_i ))
END IF
IF (.NOT. ALLOCATED(gdlambda_p)) THEN
  ALLOCATE (gdlambda_p ( 1-halo_i : global_row_length+halo_i ))
END IF
IF (.NOT. ALLOCATED(gdlambda_u)) THEN
  ALLOCATE (gdlambda_u ( 1-halo_i : global_row_length+halo_i ))
END IF
IF (.NOT. ALLOCATED(grecip_dlamp)) THEN
  ALLOCATE (grecip_dlamp ( 1-halo_i : global_row_length+halo_i))
END IF
IF (.NOT. ALLOCATED(grecip_dlamu)) THEN
  ALLOCATE (grecip_dlamu ( 1-halo_i : global_row_length+halo_i))
END IF

! Also required to avoid passing unallocated arrays to atmos_physics2:
IF (.NOT. ALLOCATED(dphi_p) ) THEN
  ALLOCATE(dphi_p(1,1))
END IF
IF (.NOT. ALLOCATED(wt_phi_p) ) THEN
  ALLOCATE(wt_phi_p(1,1))
END IF
IF (.NOT. ALLOCATED(wt_phi_v) ) THEN
  ALLOCATE(wt_phi_v(1,1))
END IF
IF (.NOT. ALLOCATED(wt_lambda_p) ) THEN
  ALLOCATE(wt_lambda_p(1))
END IF
IF (.NOT. ALLOCATED(wt_lambda_u) ) THEN
  ALLOCATE(wt_lambda_u(1))
END IF

IF ( l_regular ) THEN
  CALL set_horiz_grid( )
ELSE
  CALL set_var_horiz_grid(a_len2_coldepc,a_len2_rowdepc,                  &
                        Lambda_p_in, Lambda_u_in, Phi_p_in, Phi_v_in )
END IF
CALL calc_global_grid_spacing()

!============================================================================= 
!============== INITIALISE ND GRID VARIABLES EXPECTED BY AC SCHEME ===========
!=============================================================================
! AC scheme should be updated so this code is not required in EG case.
! Allocate arrays for dyn_var_res_mod module  
! Allocate full size arrays to avoid out-of-bounds in atm_step   

ALLOCATE (phi_p(1-halo_i:row_length+halo_i,1-halo_j:rows+halo_j)) 
DO i = 1-halo_i, global_row_length+halo_i 
  glambda_p(i) = glob_xi1_p(i) 
END DO 

DO j = 1-halo_j, rows+halo_j 
  DO i = 1-halo_i, row_length+halo_i 
    phi_p(i,j) = xi2_p(j) 
  END DO 
END DO 

!============================================================================= 

! End of ENDGAME-only

!  Set trig fields and Coriolis components

! DEPENDS ON: set_trigs_4A
CALL set_trigs_4A(                                                &
               row_length, rows, n_rows,model_levels,             &
               base_lambda, base_phi,                             &
               delta_lambda, delta_phi,                           &
               base_lambda_wk, base_phi_wk,                       &
               delta_lambda_wk, delta_phi_wk,                     &
               f_plane_rad, ff_plane_rad,                         &
               two_Omega,                                         &
               lat_rot_NP, long_rot_NP,                           &
               lat_rot_NP_in, long_rot_NP_in,                     &
               lat_n, lat_s, long_e, long_w,                      &
               long_e_model, long_w_model, rmdi)

! ----------------------------------------------------------------------
! Set atm_step_consts
! ----------------------------------------------------------------------

CALL atm_step_consts(                                                   &
                     row_length, rows, model_levels,                    &
                     bl_levels, first_constant_r_rho_level )

! ----------------------------------------------------------------------
! 1. Set polar filtering and diffusion
! ----------------------------------------------------------------------

! Allocate arrays for diff_coeff_mod module
IF (.NOT. ALLOCATED(diff_coeff_u)) THEN
  ALLOCATE (diff_coeff_u ( 1-Offx : row_length+Offx,              &
                           1-Offy : rows+Offy ))
END IF
IF (.NOT. ALLOCATED(diff_coeff_v)) THEN
  ALLOCATE (diff_coeff_v ( 1-Offx : row_length+Offx,              &
                           1-Offy : n_rows+Offy ))
END IF
!
!
!  THIS NEEDS OUT_OF_BOUNDS ERROR FIXED FOR ENDGAME
!
!
IF ( L_pofil_new ) THEN
  CALL umPrint(' Newest combi-polar filter used ',src='setcona_4A')
  ! DEPENDS ON: eg_setdiff
  CALL eg_Setdiff(                                                        &
  ! Input size and control variables
              global_rows, row_length, rows, n_rows, model_levels,        &
              at_extremity, datastart, offx, offy, mype, nproc_y,         &
              model_levels_max, max_121_rows, max_sweeps,                 &
  ! other constants
              delta_lambda, delta_phi,                                    &
              polar_cap, scale_ratio,                                     &
              ref_lat_deg, diff_coeff_ref,                                &
              cos_theta_latitude, sin_theta_latitude,                     &
              cos_v_latitude, sin_v_latitude,                             &
  ! Output data
               global_u_filter, global_v_filter,                          &
               u_sweeps, v_sweeps,                                        &
               u_begin, u_end, v_begin, v_end,                            &
               diff_coeff_u, diff_coeff_v, diff_coeff_phi,                &
               diffusion_coefficient_thermo, diffusion_order_thermo,      &
               diffusion_coefficient_wind, diffusion_order_wind,          &
               diff_order_thermo, diff_order_wind,                        &
               diff_timescale_thermo, diff_timescale_wind,                &
               L_sponge, sponge_power, sponge_ew, sponge_ns,              &
               sponge_wts_ew, sponge_wts_ns, L_subfilter_horiz,           &
               L_diffusion, L_cdiffusion, L_filter, L_diff_auto,          &
               L_pfcomb, L_pfexner, L_pftheta, L_pfuv, L_pfw, L_pfincs,   &
               L_diff_thermo, L_diff_wind, L_diff_w, L_diff_incs,         &
               L_diff_wind_ew_only, halo_i, halo_j,xi2_p, xi2_v)
ELSE
  CALL umPrint('WARNING - ',src='setcona_4A')
  CALL umPrint('   only New combi-filter is compatible with ENDGAME',     &
      src='setcona_4A')
END IF

! ----------------------------------------------------------------------
! 1.5 Set filtering control variable
! ----------------------------------------------------------------------
L_diff_active = .FALSE.

IF ( L_polar_filter .OR. L_pfcomb ) THEN
  L_diff_active = .TRUE.
  L_filter = .TRUE.
  IF (model_type /= mt_global) THEN
    WRITE(umMessage,*)' ******   Polar filter NOT ALLOWED in LAMs ****** '
    CALL umPrint(umMessage,src='setcona_4A')
    WRITE(umMessage,'(A)')' ******       RESET NAMELIST       ****** '
    CALL umPrint(umMessage,src='setcona_4A')
    icode    = 1
    cmessage ='SETCONA:Polar filter NOT ALLOWED in LAMs'
    CALL Ereport(RoutineName,icode,Cmessage)
  ELSE IF ( L_pofil_new ) THEN
    WRITE(umMessage,*)' ******  Newest filtering code being used ****** '
    CALL umPrint(umMessage,src='setcona_4A')
  END IF ! model_type  /=  mt_global
ELSE
  WRITE(umMessage,*) ' ******   Polar filter is OFF   ****** '
  CALL umPrint(umMessage,src='setcona_4A')
END IF ! L_polar_filter .or. L_pfcomb
IF ( L_diff_thermo .OR. L_diff_wind .OR. L_diff_w ) THEN
  L_diff_active = .TRUE.
  IF ( L_pofil_new ) THEN
    WRITE(umMessage,*) ' ******  Newest diffusion code being used ****** '
    CALL umPrint(umMessage,src='setcona_4A')
  END IF !  L_pofil_new
  L_filter = .TRUE.
END IF ! L_diff_thermo .or. L_diff_wind .or. L_diff_w
IF ( L_diff_incs .OR. L_pfincs) THEN
  L_filter_incs = .TRUE.
END IF ! L_diff_incs .or. L_pfincs
IF ( L_pfexner  ) THEN
  L_diff_active = .TRUE.
  WRITE(umMessage,*) ' ** Polar filtering of Exner pressure is active ** '
  CALL umPrint(umMessage,src='setcona_4A')
END IF ! L_pfexner
IF ( L_diff_exner  ) THEN
  L_diff_active = .TRUE.
  WRITE(umMessage,*) ' ** Diffusion of Exner pressure is active ** '
  CALL umPrint(umMessage,src='setcona_4A')
END IF ! L_diff_exner

IF ( L_tardiff_q ) THEN
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage,*) '  ***   Targeted diffusion of q is ACTIVE  *** '
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage, 913) w_conv_limit, tardiffq_factor
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage, 914) tardiffq_test
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage, 915) tardiffq_start, tardiffq_end
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage, 916) tar_horizontal
  CALL umPrint(umMessage,src='setcona_4A')
ELSE
  WRITE(umMessage,*) ' '
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage,*) '  ***   Targeted diffusion of q is NOT ACTIVE  *** '
  CALL umPrint(umMessage,src='setcona_4A')
END IF !   L_tardiff_q

! ----------------------------------------------------------------------
! Section 1.6 Set up turb_diff control
! ----------------------------------------------------------------------
IF ( L_subfilter_horiz .OR. L_subfilter_vert ) THEN

  L_diff_active = .TRUE.

  ! DEPENDS ON: init_turb_diff
  CALL init_turb_diff(                                            &
                      model_levels, bl_levels,                    &
                      turb_startlev_horiz, turb_endlev_horiz,     &
                      turb_startlev_vert, turb_endlev_vert,       &
                      global_u_filter, global_v_filter,           &
                      diff_factor, mix_factor,                    &
                      L_subfilter_horiz, L_subfilter_vert)

END IF ! L_subfilter_horiz .OR. L_subfilter_vert

!  ----------------------------------------------------------
!  Print information to tell user what is active
!  ----------------------------------------------------------
!  vertical diffusion coefficients
IF ( L_vdiff_uv ) THEN
  vdiffuv_factor = 0.25 *                                         &
                     (1.0 - EXP( -1.0 / vdiffuv_timescale) )
  WRITE(umMessage,*) 'A factor delta_z**2/timestep is '           &
  ,' allowed for in the application of the vertical diffusion.'
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage,*) 'Vertical diffusion of order 1 with'         &
         , ' vertical diffusion coefficient = ', vdiffuv_factor
  CALL umPrint(umMessage,src='setcona_4A')
  WRITE(umMessage,*) 'Vertical diffusion timescale is '           &
                            , vdiffuv_timescale,' steps '
  CALL umPrint(umMessage,src='setcona_4A')
END IF ! L_vdiff_uv

IF ( top_filt_start > model_levels .OR. top_filt_start==imdi ) THEN
  WRITE(umMessage,*) 'No additional upper-level diffusion'
  CALL umPrint(umMessage,src='setcona_4A')
ELSE  ! top_filt_start <= model_levels
  active_levels = top_filt_end - top_filt_start + 1
  IF ( active_levels > max_updiff_levels) THEN
    top_filt_start  = top_filt_end - max_updiff_levels + 1
    WRITE(umMessage,*) ' Max of uppermost ', max_updiff_levels    &
          , ' levels allowed for upper-level'                     &
          , ' diffusion - start level reset to ', top_filt_start
    CALL umPrint(umMessage,src='setcona_4A')
    active_levels = max_updiff_levels
  END IF ! active_levels > max_updiff_levels
  WRITE(umMessage,*) 'Extra upper-level diffusion applied'
  CALL umPrint(umMessage,src='setcona_4A')
  SCALE = 1.0
  IF ( L_upper_ramp ) THEN
    SCALE = up_diff_scale
    WRITE(umMessage, 911) top_diff, top_filt_end
    CALL umPrint(umMessage,src='setcona_4A')
    WRITE(umMessage, 912) SCALE, top_filt_start
    CALL umPrint(umMessage,src='setcona_4A')
  ELSE  !
    WRITE(umMessage, 910) top_diff, top_filt_start, top_filt_end
    CALL umPrint(umMessage,src='setcona_4A')
  END IF  !  L_upper_ramp
  kk = active_levels
  up_diff(active_levels) = top_diff
  DO k = top_filt_end - 1, top_filt_start, - 1
    kk = kk - 1
    up_diff(kk) = SCALE * up_diff(kk + 1)
  END DO !  k = top_filt_end - 1, top_filt_start, - 1
  DO k = 1, active_levels
    kk = k + top_filt_start - 1
    WRITE(umMessage,*)'Level', kk,' Diffusion factor = ',up_diff(k)
    CALL umPrint(umMessage,src='setcona_4A')
  END DO !  k = 1, active_levels
END IF ! top_filt_start > model_levels

! ----------------------------------------------------------------------
! Section 1.5a . Print FORMATTING
! ----------------------------------------------------------------------

910  FORMAT(' Diffusion coefficient = ',f5.2                           &
           ,' from level ', i3,' to level ', i3)
911  FORMAT(' Diffusion coefficient = ',f5.2,' at level ', i3)
912  FORMAT(' decreasing by ', f6.3, ' at each level down to ', i3)
913  FORMAT('   Vertical velocity test value is  ',f5.2,' m/s and ',   &
                                      'diffusion factor = ',f6.3)
914  FORMAT('   Start testing at level  ',i3)
915  FORMAT('   Apply from level  ',i3,' to level ',i3)
916  FORMAT('   Slope test applied up to level ',i4)

! ----------------------------------------------------------------------
! 2. Semi-Lagrangian scheme options.
! ----------------------------------------------------------------------
! j holds number of points required for interpolation scheme
j = 2
DO i = 1, 3
  IF (high_order_scheme(i)  ==  quinticLagrange) THEN
    j = 3
    IF (halo_j  <   5) THEN
      WRITE(umMessage,*)' error halo_j too small ',halo_j
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
  ELSE
    IF (halo_j  <   4) THEN
      WRITE(umMessage,*)' error halo_j too small ',halo_j
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
  END IF
END DO

! Calculate max cfl possible in LAM
! Y-direction also applies in global model though variable is not
! used in code.

LAM_max_cfl(1) = halo_i - j
LAM_max_cfl(2) = halo_j - j

! ----------------------------------------------------------------------
! Set up data.
! ----------------------------------------------------------------------

! Calculate land_index
land_points = 0
DO j =1, rows
  DO i = 1, row_length
    IF (land_sea_mask(i,j)) THEN
      land_points = land_points + 1
      land_index(land_points) = (j-1)*row_length + i
    END IF
  END DO
END DO

! set-up code for hydrology
soil_points = 0
land_ice_points = 0
l_soil_point(:) = .FALSE.
l_lice_point(:) = .FALSE.
DO i= 1, land_points
  ! Test on soil moisture concentration at saturation
  IF (smvcst(i) >   0.0) THEN       ! Soil points
    soil_points                     = soil_points+1
    soil_index(soil_points)         = i
    l_soil_point(i)                 = .TRUE.
  ELSE IF (smvcst(i)  ==  0.0) THEN   ! Land-ice points
    land_ice_points                 = land_ice_points+1
    land_ice_index(land_ice_points) = i
    l_lice_point(i)                 = .TRUE.
  END IF
END DO

! call swap_bounds to set extra points
i_field = 0
i_field = i_field + 1
fields_to_swap(i_field) % field      => r_at_u(:,:,:)
fields_to_swap(i_field) % field_type = fld_type_u
fields_to_swap(i_field) % levels     = model_levels
fields_to_swap(i_field) % rows       = rows
fields_to_swap(i_field) % vector     = .FALSE.

i_field = i_field + 1
fields_to_swap(i_field) % field      => r_at_v(:,:,:)
fields_to_swap(i_field) % field_type = fld_type_v
fields_to_swap(i_field) % levels     = model_levels
fields_to_swap(i_field) % rows       = n_rows
fields_to_swap(i_field) % vector     = .FALSE.

CALL swap_bounds_mv(fields_to_swap, i_field, row_length,         &
                       halo_i, halo_j)

!  Initialise diagnostic printing items
CALL init_print_diag()

!  Set values for block printing:
!               calls for printing need to be inserted by user
! DEPENDS ON: init_block4_pr
CALL init_block4_pr(                                            &
                    global_row_length, global_rows,             &
                    blockx_in, blocky_in,                       &
                    dom_w_in, dom_e_in, dom_s_in, dom_n_in)

!  ----------------------------------------------------------
!  Set number of levels for non-local PBL scheme
!  ----------------------------------------------------------
IF (i_bl_vn == i_bl_vn_9b .OR. i_bl_vn == i_bl_vn_9c) THEN
  IF (z_nl_bl_levels > 0.0 .AND. z_nl_bl_levels < z_top_of_model) THEN 
  ! nl_bl_levels is the first theta-level at or above z_nl_bl_levels
  ! and no greater than boundary_layer_levels
    k=1
    DO WHILE ( k < boundary_layer_levels .AND.                          &
               r_ref_theta(k) < z_nl_bl_levels )
      k=k+1
    END DO
    nl_bl_levels = k 
    IF (printstatus  >=  prstatus_normal) THEN
      WRITE(umMessage,'(A,ES12.4,A,I4)')                                &
                   ' z_nl_bl_levels=',z_nl_bl_levels,                   &
                   ' which gives nl_bl_levels = ',nl_bl_levels
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
  ELSE
    icode    = 1
    WRITE(cmessage,'(A,ES12.4,A)') 'SETCONA: Parameter z_nl_bl_levels=',&
           z_nl_bl_levels,' is outside allowed range'
    CALL Ereport(RoutineName,icode,cmessage)

  END IF ! test on z_nl_bl_levels
END IF

!  ----------------------------------------------------------
!  Set the values of alpha_cd
!  ----------------------------------------------------------
ALLOCATE(alpha_cd(bl_levels))
alpha_cd(1) = alpha_cd_in(1)
alpha_cd(2:bl_levels) = alpha_cd_in(2)

!  ----------------------------------------------------------
!  Set number of levels for stochastic perturbations
!  ----------------------------------------------------------
IF (i_pert_theta > off) THEN
  IF (zmin_pert_theta >= 0.0 .AND.                                &
      zmin_pert_theta < z_top_of_model) THEN 
    ! minlev_pert_theta is the highest theta-level at or below 
    ! zmin_pert_theta and no greater than model_levels
    k=2
    DO WHILE ( k < model_levels .AND.                             &
               r_ref_theta(k) < zmin_pert_theta )
      k=k+1
    END DO
    minlev_pert_theta = k-1
    IF (printstatus  >=  prstatus_normal) THEN
      WRITE(umMessage,'(A,ES12.4,A,I4)')                          &
          ' zmin_pert_theta=',zmin_pert_theta,                    &
          ' which gives minlev_pert_theta = ',minlev_pert_theta
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
  ELSE
    icode    = 1
    WRITE(cmessage,'(A,ES12.4,A)')                                &
                    'setcona_4A: Parameter zmin_pert_theta=',     &
                    zmin_pert_theta,' is outside allowed range'
    CALL Ereport(RoutineName,icode,cmessage)
  END IF ! test on zmin_pert_theta

  IF (zmax_pert_theta > 0.0 .AND.                                 &
      zmax_pert_theta < z_top_of_model) THEN 
    ! maxlev_pert_theta is the highest theta-level at or below 
    ! zmax_pert_theta and no greater than model_levels
    k=2
    DO WHILE ( k < model_levels .AND.                             &
               r_ref_theta(k) < zmax_pert_theta )
      k=k+1
    END DO
    maxlev_pert_theta = k-1
    IF (printstatus  >=  prstatus_normal) THEN
      WRITE(umMessage,'(A,ES12.4,A,I4)')                          &
             ' zmax_pert_theta=',zmax_pert_theta,                 &
             ' which gives maxlev_pert_theta = ',maxlev_pert_theta
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
  ELSE
    icode    = 1
    WRITE(cmessage,'(A,ES12.4,A)')                                &
                    'setcona_4A: Parameter zmax_pert_theta=',     &
                    zmax_pert_theta,' is outside allowed range'
    CALL Ereport(RoutineName,icode,cmessage)
  END IF ! test on zmax_pert_theta

  IF (maxlev_pert_theta < minlev_pert_theta) THEN 
    ! will give no perturbations so stop
    icode    = 1
    WRITE(cmessage,'(A)')                                        &
        'setcona_4A: maxlev_pert_theta is less than '//          &
        'minlev_pert_theta so no perturbations will be generated'
    CALL Ereport(RoutineName,icode,cmessage)
  END IF ! test on maxlev_pert_theta

END IF ! test on i_pert_theta
!  ----------------------------------------------------------
!  3. Set up cloud type boundaries for low/medium/high cloud.
!  ----------------------------------------------------------
IF (num_cloud_types  >   3) THEN
  icode    = 1
  cmessage = 'SETCONA: Parameter NUM_CLOUD_TYPES exceeds 3'
  WRITE(umMessage,*) 'NUM_CLOUD_TYPES=',num_cloud_types
  CALL umPrint(umMessage,src='setcona_4A')
  CALL Ereport(RoutineName,icode,Cmessage)

END IF

!  Diagnostics of cloud amount for low/medium/high cloud are calculated
!  by finding the max cloud cover over a range of model levels. The
!  ranges are determined by heights [ 1->2:2->3:3->4 = L:M:H ] held in
!  array h_split by comparison with heights of model levels (above
!  surface) in r_ref_theta.

DO kk = 1, num_cloud_types + 1
  level = 1
  !
  !       r_ref_theta turns out to be A_theta(k) - planet_radius at every
  !       model level because z_top_of_model = z_rho_top is chosen here.
  !
  DO WHILE (r_ref_theta(level)  <=  h_split(kk))
    level = level + 1

    IF (level == model_levels) THEN
      EXIT
    END IF
  END DO

  IF (level  >   model_levels) THEN
    icode    = -1
    cmessage ='SETCONA:Error in locating levels for cloud layers'
    CALL Ereport(RoutineName,icode,Cmessage)

  END IF
  cloud_bound(kk) = level

END DO

low_bot_level  = cloud_bound(1)
low_top_level  = cloud_bound(2) - 1
med_bot_level  = cloud_bound(2)
med_top_level  = cloud_bound(3) - 1
high_bot_level = cloud_bound(3)
high_top_level = cloud_bound(4) - 1

IF (low_top_level  >   cloud_levels) THEN
  icode    = 1
  cmessage = 'SETCONA: No of cloud levels less than Top of Low'
  CALL Ereport(RoutineName,icode,Cmessage)

END IF

IF (med_top_level >  cloud_levels) THEN
  icode    = 1
  cmessage = 'SETCONA:  No of cloud levels less than Top of Med'
  CALL Ereport(RoutineName,icode,Cmessage)

END IF

IF (high_top_level  >   cloud_levels) THEN
  icode    = 1
  cmessage = 'SETCONA: No of cloud levels less than Top of High'
  CALL Ereport(RoutineName,icode,Cmessage)

END IF

!  ----------------------------------------------------------
!  4. Set up locally-derived parameters.
!  ----------------------------------------------------------

!     The tropopause diagnosed for radiative purposes
!     divides theta-levels considered to lie in the stratosphere
!     from those considered to lie in the troposphere: the
!     tropopause is therefore taken as a rho-level. This level
!     is constrained to lie between heights of z_min_trop
!     and z_max_trop. The level is used in the setting of the
!     background aerosol climatology and of ozone profiles,
!     subject to the choice of appropriate options; additionally
!     it is used in the calculation of diagnostics defined at
!     the tropopause.
!
!     Start at the second rho-level because the first is omitted
!     from the grid seen in the physics.
min_trop_level=2
DO ; IF ( (r_ref_rho(min_trop_level) >= z_min_trop) .OR.          &
          (min_trop_level == model_levels) ) EXIT
  min_trop_level = min_trop_level+1
END DO
!
max_trop_level=min_trop_level
DO ; IF ( (r_ref_rho(max_trop_level) > z_max_trop) .OR.           &
          (max_trop_level == model_levels) ) EXIT
  max_trop_level = max_trop_level+1
END DO
max_trop_level = max_trop_level-1

!  ----------------------------------------------------------
!  5.0 Set up tracer advection parameters
!  ----------------------------------------------------------

! set up super tracer array size

super_array_size=0

IF (l_CO2_interactive) super_array_size = super_array_size+1
IF (l_Soot) super_array_size = super_array_size+3
IF (l_Biomass) super_array_size = super_array_size+3
IF (l_Sulpc_so2)super_array_size = super_array_size+4
IF (L_sulpc_nh3)super_array_size = super_array_size+1
IF (L_sulpc_dms)super_array_size = super_array_size+1
IF (l_dust) THEN
  IF (l_twobin_dust) THEN
    super_array_size = super_array_size + 2
  ELSE
    super_array_size = super_array_size + 6
  END IF
END IF
IF (L_ocff) super_array_size = super_array_size+3
IF (L_nitrate) super_array_size=super_array_size+2
IF (L_Murk_advect) super_array_size = super_array_size + 1
IF (l_use_cariolle) super_array_size = super_array_size + 1

! ----------------------------------------------------------------------
! Section 5.1   limit for tracer_level required
! ----------------------------------------------------------------------

IF ( tr_levels == model_levels ) THEN

  IF ( tr_vars > 0 ) super_array_size = super_array_size + tr_vars

  IF ( tr_ukca > 0 ) super_array_size = super_array_size + tr_ukca

ELSE   !  tr_levels /= model_levels

  IF (tr_ukca > 0 .OR. tr_vars > 0 ) THEN

    WRITE(umMessage,*) 'Tracer Levels: tr_levels =',tr_levels
    CALL umPrint(umMessage,src='setcona_4A')
    WRITE(umMessage,*) 'Model Levels: model_levels =', model_levels
    CALL umPrint(umMessage,src='setcona_4A')
    icode = 1
    cmessage='Tracer levels must equal model levels, when '  //  &
             'running with Free or UKCA tracers'

    CALL Ereport(RoutineName, icode, cmessage)

  END IF   ! tr_ukca > 0 .or. tr_vars > 0

END IF ! tr_levels == model_levels

IF (mype == 0) THEN
  WRITE(umMessage,*)'tracer super_array_size =  ',super_array_size
  CALL umPrint(umMessage,src='setcona_4A')
END IF

!  ----------------------------------------------------------
!  6. Set up moisture array sizes    always do q,qcl,qcf
!  ----------------------------------------------------------

moisture_array_size = 3
IF (i_cld_vn==i_cld_pc2) moisture_array_size = moisture_array_size + 4
IF (l_mcr_qrain)         moisture_array_size = moisture_array_size + 1
IF (l_mcr_qcf2 )         moisture_array_size = moisture_array_size + 1
IF (l_mcr_qgraup)        moisture_array_size = moisture_array_size + 1
IF (l_conv_prog_group_1) moisture_array_size = moisture_array_size + 1
IF (l_conv_prog_group_2) moisture_array_size = moisture_array_size + 1
IF (l_conv_prog_group_3) moisture_array_size = moisture_array_size + 1
IF (l_conv_prog_precip)  moisture_array_size = moisture_array_size + 1

IF (mype == 0) THEN
  WRITE(umMessage,*)'moisture_array_size =  ',moisture_array_size
  CALL umPrint(umMessage,src='setcona_4A')
END IF

IF (bdy_tke == mymodel3) THEN
  turb_array_size=4
ELSE
  turb_array_size=1
END IF

! ----------------------------------------------------------------
! 7. Calculate the coeffs for the interpolation of the radiation
!    quantities if spatial degradation of the radiation code is on.
!    Calculate the mask which ensures a chequer-board pattern over
!    the whole domain (and not just for a single PE).
! -----------------------------------------------------------------

! Create a land mask with a halo.
IF (L_rad_deg) THEN

  DO j=1,rows
    DO i=1,row_length
      landmask(i,j)=Land_sea_mask(i,j)
    END DO
  END DO

  CALL swap_bounds(landmask,row_length,rows,1,                    &
                   offx,offy,fld_type_p, swap_field_is_scalar)

  first_row=1
  last_row=rows
  IF (model_type == mt_global) THEN
    IF (at_extremity(PNorth)) THEN
      last_row=rows-1
    END IF
    IF (at_extremity(PSouth)) THEN
      first_row=2
    END IF
  END IF

  ! Allocate arrays for rad_mask_trop_mod module
  IF (.NOT. ALLOCATED(es_space_interp)) THEN
    ALLOCATE (es_space_interp(4, row_length, rows))
  END IF
  IF (.NOT. ALLOCATED(rad_mask)) THEN
    ALLOCATE (rad_mask(row_length, rows))
  END IF

  CALL coeffs_degrade(es_space_interp, landmask,                  &
                      first_row, last_row,                        &
                      at_extremity(PNorth), at_extremity(PSouth), &
                      at_extremity(PWest), at_extremity(PEast),   &
                      row_length*rows,row_length, rows,offx, offy)
  CALL rad_degrade_mask(rad_mask, datastart,                      &
                      Ndim_max, first_row, last_row, offx,        &
                      row_length, rows)

ELSE
  IF (.NOT. ALLOCATED(es_space_interp)) THEN
    ALLOCATE (es_space_interp(1,1,1))
  END IF
  IF (.NOT. ALLOCATED(rad_mask)) THEN
    ALLOCATE (rad_mask(1,1))
  END IF
END IF  !  If L_rad_deg loop

!  ----------------------------------------------------------
!  8. Set up dynamical core parameters
!  ----------------------------------------------------------
!     Parameters for dynamical core tests may be set here.

! ----------------------------------------------------------------------
! 9. GCR diagnostics initialisation.
! ----------------------------------------------------------------------
IF (GCR_diagnostics == 3 ) THEN
  L_error = .FALSE.
  DO i = 1, 3
    IF (GCR_its_avg_step(i) == 0 )L_error = .TRUE.
  END DO  !  i = 1, 3
  IF ( L_error ) THEN
    WRITE(umMessage,*) 'WARNING GCR iteration counting at'  &
,         ' timestep 0 or interval of 0 NOT PERMITTED '
    CALL umPrint(umMessage,src='setcona_4A')
    ! Following values will output iteration info at 6, 12 and
    !  1440 timesteps (30 days for a 30 minute timestep)
    GCR_its_avg_step(1) = 6
    GCR_its_avg_step(2) = 12
    GCR_its_avg_step(3) = 1440
    WRITE(umMessage,*)' Iteration count diagnostics reset for timesteps '&
          , GCR_its_avg_step(1), GCR_its_avg_step(2)                     &
                               , GCR_its_avg_step(3)
    CALL umPrint(umMessage,src='setcona_4A')
    WRITE(umMessage,*)  '  and at intervals of ', GCR_its_avg_step(3),   &
            ' timesteps thereafter.'
    CALL umPrint(umMessage,src='setcona_4A')
    WRITE(umMessage,'(A)')   &
            ' Change in namelist if you want different values'
    CALL umPrint(umMessage,src='setcona_4A')
  ELSE ! L_error .false.
    WRITE(umMessage,*)  ' '
    CALL umPrint(umMessage,src='setcona_4A')
    WRITE(umMessage,*)  'Iteration count diagnostics at timesteps '      &
          , GCR_its_avg_step(1), GCR_its_avg_step(2)                     &
          , GCR_its_avg_step(3)
    CALL umPrint(umMessage,src='setcona_4A')
    WRITE(umMessage,*)  ' and at intervals of ', GCR_its_avg_step(3),    &
            ' timesteps thereafter.'
    CALL umPrint(umMessage,src='setcona_4A')
  END IF ! L_error
  GCR_its_switch =  1
  GCR_sum_its = 0
  GCR_min_its = 1000
  GCR_max_its = 0
  GCR_max_time = 0
  GCR_min_time = 0
END IF ! GCR_diagnostics == 3

! ----------------------------------------------------------------------
! 10. Error trapping for advection choices
! ----------------------------------------------------------------------

DO i=1,3
  IF (.NOT. L_mono(i) .AND. .NOT. L_high(i)) THEN
    IF (i == 1) THEN
      WRITE(umMessage,*)                            &
      'WARNING Advection choices incompatible for theta'
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
    IF (i == 2) THEN
      WRITE(umMessage,*)                            &
      'WARNING Advection choices incompatible for moisture'
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
    IF (i == 3) THEN
      WRITE(umMessage,*)                            &
      'WARNING Advection choices incompatible for winds'
      CALL umPrint(umMessage,src='setcona_4A')
    END IF
    WRITE(umMessage,*)                              &
      'Both L_high and L_mono switches set to false '
    CALL umPrint(umMessage,src='setcona_4A')

  END IF
END DO

! ----------------------------------------------------------------------
! 11. Orographic correction parameters for the radiation code.
! ----------------------------------------------------------------------

IF (l_orog) THEN

  ! Calculate (rotated) grid coordinates with a halo out to the
  ! horizon limit, to be used in calculation of gridpoint separation.
  IF (l_regular) THEN
    DO j = 1-horiz_limit, rows+horiz_limit
      gj = datastart(2) + j - 1
      phi_horiz(j) = base_phi+(gj-1)*delta_phi
    END DO
    dphi_horiz=delta_phi
    dlambda_horiz=delta_lambda
  ELSE
    DO i = 1-horiz_limit, row_length+horiz_limit
      gi = MAX(MIN(datastart(1)+i-1,global_row_length-1),1)
      dlambda_horiz(i) = (lambda_p_in(gi+1)-lambda_p_in(gi))      &
        * pi_over_180
    END DO
    DO j = 1-horiz_limit, rows+horiz_limit
      gj = MAX(MIN(datastart(2)+j-1,global_rows-1),1)
      dphi_horiz(j) = (phi_p_in(gj+1)-phi_p_in(gj))*pi_over_180
    END DO
    DO j = 1-horiz_limit, rows+horiz_limit
      gj = datastart(2) + j - 1
      IF (gj >= 1 .AND. gj <= global_rows) THEN
        phi_horiz(j) = phi_p_in(gj)*pi_over_180
      ELSE IF (gj < 1) THEN
        phi_horiz(j) = phi_p_in(1)*pi_over_180 +                  &
          (gj-1)*dphi_horiz(j)
      ELSE
        phi_horiz(j) = phi_p_in(global_rows)*pi_over_180 +        &
          (gj-global_rows)*dphi_horiz(j)
      END IF
    END DO
  END IF


  ! Calculate local bearing of the rotated pole.
  IF (model_type == mt_global) THEN
    bear_rot_NP = 0.0
  ELSE
    CALL pole_bearing(row_length, rows,                           &
      lat_rot_NP, long_rot_NP, bear_rot_NP)
  END IF


  IF (l_skyview) THEN
    ! Set orographic heights for calculation of horizon angles
    ! (using unsmoothed height data if available).
    orog_horiz=0.0
    IF (l_orog_unfilt) THEN
      ! DEPENDS ON: all_gather_field
      CALL all_gather_field(orog_unfilt, orog_global,             &
        row_length, rows, global_row_length, global_rows,         &
        fld_type_p, halo_type_no_halo, gc_all_proc_group,         &
        icode, cmessage)
    ELSE
      CALL all_gather_field(orog, orog_global,                    &
        row_length, rows, global_row_length, global_rows,         &
        fld_type_p, halo_type_no_halo, gc_all_proc_group,         &
        icode, cmessage)
    END IF
    DO j = 1-horiz_limit, rows+horiz_limit
      gj = datastart(2) + j - 1
      IF (gj >= 1 .AND. gj <= global_rows) THEN
        DO i = 1-horiz_limit, row_length+horiz_limit
          gi = datastart(1) + i - 1
          IF (gi >= 1 .AND. gi <= global_row_length) THEN
            orog_horiz(i,j) = orog_global(gi,gj)
          ELSE IF (model_type == mt_global) THEN
            IF (gi < 1) THEN
              orog_horiz(i,j)=orog_global(gi+global_row_length,gj)
            ELSE
              orog_horiz(i,j)=orog_global(gi-global_row_length,gj)
            END IF
          END IF
        END DO
      END IF
    END DO
  END IF


  ! Calculate slope aspect and angle.
  IF (l_use_grad_corr) THEN
    ! Use ancillary slopes calculated from the raw orography.
    CALL aspang_ancil(row_length, rows, land_points,              &
      land_sea_mask, grad_x, grad_y, bear_rot_NP)
  ELSE IF (l_skyview) THEN
    ! Use the same orographic heights as used for the
    ! calculation of the horizon angles.
    CALL aspang(row_length, rows, bear_rot_NP,                    &
      phi_horiz, dphi_horiz, dlambda_horiz,                       &
      orog_horiz(0:row_length+1,0:rows+1))
  ELSE
    ! Use the model orography (generally smoothed).
    CALL aspang(row_length, rows, bear_rot_NP,                    &
      phi_horiz, dphi_horiz, dlambda_horiz,                       &
      orog_halos(0:row_length+1,0:rows+1))
  END IF

  IF (model_type == mt_global) THEN
    IF (at_extremity(PNorth)) slope_angle(:,rows) = 0.0
    IF (at_extremity(PSouth)) slope_angle(:,1) = 0.0
  END IF


  IF (l_skyview) THEN
    ! Calculate horizon angles and sky-view correction factor.
    CALL skyview(row_length, rows, bear_rot_NP,                   &
      phi_horiz, dphi_horiz, dlambda_horiz, orog_horiz)
  ELSE
    ALLOCATE(horiz_ang(1,1,1,1))
    ALLOCATE(horiz_aspect(1,1,1))
  END IF

ELSE

  ALLOCATE(slope_aspect(1,1))
  ALLOCATE(slope_angle(1,1))
  ALLOCATE(horiz_ang(1,1,1,1))
  ALLOCATE(horiz_aspect(1,1,1))

END IF


! Additional idealised grid variables
IF ( l_idealised_data ) THEN

  IF (.NOT. ALLOCATED(z_ref_theta)) THEN
    ALLOCATE (z_ref_theta(0:model_levels))
  END IF
  IF (.NOT. ALLOCATED(z_ref_rho)) THEN
    ALLOCATE (z_ref_rho(model_levels))
  END IF

  z_ref_theta(:) = z_theta_in(:)
  z_ref_rho(:) = z_rho_in(:)

END IF

! Determine if this is a zero-orography run for idealised models:
IF (problem_number /= standard) THEN
  max_orog = MAXVAL(orog)
  CALL gcg_rmax(1, gc_all_proc_group, icode, max_orog)
  IF (max_orog < EPSILON(1.0)) THEN
    zero_orography = .TRUE.
  END IF
END IF

CALL eg_init_lam_domain_kind(halo_i, halo_j, row_length, rows,       &
                        datastart, global_rows, global_row_length)

! Initialise timestep and timestep parameters
!
WRITE(umMessage,*)  ' '
CALL umPrint(umMessage,src='setcona_4A')

CALL atm_step_time_init

! Set constants for sisl

CALL eg_sisl_setcon(                                                 &
       row_length, rows, n_rows, model_levels, halo_i, halo_j,       &
       offx, offy,  datastart,                                       &
       l_shallow, l_const_grav,                                      &
       z_top_of_model, g_theta,  g_rho,                              &
       pole_consts)

IF (lhook) CALL dr_hook(RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Setcona_4A
