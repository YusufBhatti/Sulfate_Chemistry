! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Perform idealised initialisations

MODULE rcf_ideal_initialisation_mod

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_INITIALISATION_MOD'

CONTAINS

! Subroutine rcf_ideal_initialisation
!
! Description:
!   Main control routine for initialising idealised dumps. From here
!   the appropriate initialisation routines are chosen, depending
!   on the details of the idealised problem selected.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_initialisation ( fields_out, field_count_out, hdr_out)

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Locate_Mod, ONLY: &
    rcf_locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    rcf_alloc_field,           &
    rcf_dealloc_field

USE Rcf_Write_Field_Mod, ONLY: &
    rcf_write_field

USE Rcf_Generate_Heights_Mod, ONLY: &
    rcf_generate_heights

USE rcf_interp_weights_mod, ONLY: &
    intw_rho2w,                   &
    intw_w2rho

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE rcf_ideal_baroclinic_mod, ONLY: &
    rcf_ideal_baroclinic

USE rcf_ideal_baroclinic_channel_mod, ONLY: &
    rcf_ideal_baroclinic_channel

USE rcf_ideal_rotating_solid_body_mod, ONLY: &
    rcf_ideal_rotating_solid_body

USE rcf_ideal_deep_baroclinic_mod, ONLY: &
    rcf_ideal_deep_baroclinic

USE rcf_ideal_generic_temperature_profile_mod, ONLY: &
    rcf_ideal_generic_temperature_profile

USE rcf_ideal_bubble_mod, ONLY: &
    rcf_ideal_bubble

USE rcf_ideal_isothermal_profile_mod, ONLY: &
    rcf_ideal_isothermal_profile

USE rcf_ideal_random_perturbation_mod, ONLY: &
    rcf_ideal_random_perturbation

USE rcf_ideal_thermodynamic_profile_mod, ONLY: &
    rcf_ideal_thermodynamic_profile

USE rcf_ideal_hydrostatic_balance_mod, ONLY: &
    rcf_ideal_hydrostatic_balance

USE rcf_ideal_vapour_profile_mod, ONLY: &
    rcf_ideal_vapour_profile

USE rcf_ideal_vapour_calc_mod, ONLY: &
    rcf_ideal_vapour_calc

USE rcf_recompute_wet_rho_mod, ONLY: &
    rcf_recompute_wet_rho

USE Rcf_theta_t_Convs_mod, ONLY: &
    Rcf_Conv_t_theta

USE decomp_params, ONLY: &
    decomp_rcf_output

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,     &
    stashcode_orog,         &
    stashcode_u,            &
    stashcode_v,            &
    stashcode_w,            &
    stashcode_theta,        &
    stashcode_thetavd,      &
    stashcode_dry_rho,      &
    stashcode_exner,        &
    stashcode_exner_surf,   &
    stashcode_mv

USE cppxref_mod, ONLY: &
    ppx_atm_cuall,     &
    ppx_atm_tall

USE stparam_mod, ONLY:      &
    st_levels_model_rho,    &
    st_levels_model_theta

USE nlstcall_mod, ONLY: &
    LTimer

USE rcf_nlist_recon_idealised_mod, ONLY: &
    base_xi1,                            &
    base_xi2,                            &
    delta_xi1,                           &
    delta_xi2,                           &
    idl_bubble_option,                   &
    initial_profile,                     &
    l_perturb_q,                         &
    l_perturb_t,                         &
    mv_init_data,                        &
    mv_init_field_type,                  &
    mv_init_height,                      &
    num_mv_init_heights,                 &
    num_theta_init_heights,              &
    num_uv_init_heights,                 &
    perturb_magnitude_q,                 &
    perturb_magnitude_t,                 &
    pressure_balance,                    &
    surface_type,                        &
    theta_init_data,                     &
    theta_init_field_type,               &
    theta_init_height,                   &
    tprofile_number,                     &
    uv_init_height,                      &
    u_init_data,                         &
    v_init_data

USE idealise_run_mod, ONLY: &
    l_const_grav,           &
    l_shallow,              &
    p_surface,              &
    t_surface

USE prof_interp_mod, ONLY: &
    extrap_constant,       &
    interp_linear,         &
    prof_interp

USE profiles_mod, ONLY:     &
    mv_init,                &
    num_data_max,           &
    theta_init,             &
    u_init,                 &
    v_init

USE copy_profile_mod, ONLY: &
    copy_profile

USE model_domain_mod, ONLY: &
    l_regular,              &
    model_type,             &
    mt_global,              &
    l_cartesian,            &
    mt_cyclic_lam

USE dynamics_testing_mod, ONLY: &
    l_dry,                      &
    problem_number

USE problem_mod, ONLY: &
    dynamical_core,    &
    idealised_planet,  &
    idealised_problem

USE planet_constants_mod, ONLY: &
    g,                          &
    kappa,                      &
    planet_radius,              &
    p_zero,                     &
    r,                          &
    recip_epsilon

USE rcf_ideal_tprofile_mod, ONLY: &
    tp_isothermal

USE rcf_ideal_pprofile_constants_mod, ONLY: &
    moist_hydro,                            &
    no_balance

USE rcf_ideal_surface_mod, ONLY: &
    surface_zero

USE rcf_ideal_initial_profiles_mod, ONLY: &
    no_preset,                            &
    baro_inst,                            &
    rot_solid_body,                       &
    deep_baro_inst

USE missing_data_mod, ONLY: &
    imdi,                   &
    rmdi

USE ereport_mod, ONLY: &
    ereport

USE errormessagelength_mod, ONLY: &
    errormessagelength

USE umPrintMgr, ONLY: &
    umPrint,          &
    umMessage,        &
    PrintStatus,      &
    PrStatus_Diag

USE um_parcore, ONLY: &
    mype

USE lam_config_inputs_mod, ONLY: &
    frstlona,    frstlata,       &
    delta_lon,   delta_lat

USE Rcf_Setup_RealC_Mod, ONLY: &
    Rcf_Setup_RealC

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER,            INTENT(IN)       :: field_count_out  ! no. of output fields
TYPE( field_type ), POINTER          :: fields_out( : )  ! ptr to output fields
TYPE( um_header_type ), INTENT(IN)   :: hdr_out          ! output header

! Local variables
INTEGER                            :: i
INTEGER                            :: k
INTEGER                            :: k_top
INTEGER                            :: icode
INTEGER                            :: ErrorStatus
INTEGER                            :: pos_orog
INTEGER                            :: pos_u
INTEGER                            :: pos_v
INTEGER                            :: pos_theta
INTEGER                            :: pos_thetavd
INTEGER                            :: pos_dry_rho
INTEGER                            :: pos_exner
INTEGER                            :: pos_exner_surf
INTEGER                            :: pos_mv
REAL, ALLOCATABLE                  :: xi3_at_rho(:,:)
REAL, ALLOCATABLE                  :: xi3_at_theta(:,:)
REAL, ALLOCATABLE                  :: xi3_at_u(:,:)
REAL, ALLOCATABLE                  :: xi3_at_v(:,:)
REAL, ALLOCATABLE                  :: g_theta(:,:)
REAL, ALLOCATABLE                  :: z_theta_1d(:)! Height of theta_levels(1,:)
                                                   ! above surface
REAL, ALLOCATABLE                  :: z_rho_1d(:)  ! Height of rho_levels(1,:)
                                                   ! above surface
REAL, ALLOCATABLE                  :: g_theta_1d(:)  ! g_theta(1,:)
REAL, ALLOCATABLE                  :: theta_1d(:)    ! theta(1,:)
REAL, ALLOCATABLE                  :: thetavd_1d(:)  ! thetavd(1,:)
REAL, ALLOCATABLE                  :: u_1d(:)    ! Sample 1d profile for u
REAL, ALLOCATABLE                  :: v_1d(:)    ! Sample 1d profile for v
REAL, ALLOCATABLE                  :: exner_1d(:)! Sample 1d profile for exner
REAL, ALLOCATABLE                  :: mv_1d(:)   ! Sample 1d profile for mv
REAL, ALLOCATABLE                  :: dry_rho_1d(:) ! Sample 1d profile 
                                                    ! for dry rho
REAL                               :: p1
REAL                               :: p2
REAL                               :: p_on_orog  ! surface pressure isothermally
                                                 ! adjusted for orography
LOGICAL                            :: initialised = .FALSE.

! Initial profiles specified at a single time.
INTEGER, PARAMETER                 :: num_init_times          = 1
REAL, PARAMETER                    :: init_time(num_data_max) = 0.0

CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_IDEAL_INITIALISATION'
CHARACTER (LEN=errormessagelength) :: Cmessage

TYPE( field_type ), POINTER        :: orog_out   ! pointer to output orography

INTEGER(KIND=jpim), PARAMETER      :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER      :: zhook_out = 1
REAL(KIND=jprb)                    :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (LTimer) CALL Timer( RoutineName, 3)

WRITE(umMessage, '(A)') 'Initialising ideal fields'
CALL umPrint(umMessage, src=RoutineName)

! --------------------------------------------
! Find/calculate all the fields we will need.
! --------------------------------------------

CALL rcf_locate (stashcode_prog_sec, stashcode_orog,          &
                 fields_out, field_count_out, pos_orog, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_u,             &
                 fields_out, field_count_out, pos_u, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_v,             &
                 fields_out, field_count_out, pos_v, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_theta,         &
                 fields_out, field_count_out, pos_theta, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_thetavd,       &
                 fields_out, field_count_out, pos_thetavd, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_dry_rho,       &
                 fields_out, field_count_out, pos_dry_rho, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_exner,         &
                 fields_out, field_count_out, pos_exner, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_exner_surf,    &
                 fields_out, field_count_out, pos_exner_surf, .TRUE.)
CALL rcf_locate (stashcode_prog_sec, stashcode_mv,            &
                 fields_out, field_count_out, pos_mv, .TRUE.)

! Exner, thetavd, dry_rho and u, v winds are needed in all cases.
IF (pos_exner == 0) THEN
  ErrorStatus = 10
  cmessage = 'Could not find exner in output dump'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF
IF (pos_theta == 0) THEN
  ErrorStatus = 15
  cmessage = 'Could not find theta in output dump'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF
IF (pos_thetavd == 0) THEN
  ErrorStatus = 20
  cmessage = 'Could not find thetavd in output dump'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF
IF (pos_dry_rho == 0) THEN
  ErrorStatus = 30
  cmessage = 'Could not find dry rho in output dump'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF
IF (pos_u == 0) THEN
  ErrorStatus = 40
  cmessage = 'Could not find u wind in output dump'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF
IF (pos_v == 0) THEN
  ErrorStatus = 50
  cmessage = 'Could not find v wind in output dump'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF
IF (pos_mv == 0) THEN
  ErrorStatus = 60
  cmessage = 'Could not find mv in output dump'
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

! Rho and theta heights are also likely to be needed each time.
! Orog must be in the output dump.
CALL rcf_locate( stashcode_prog_sec, stashcode_orog,           &
               fields_out, field_count_out, pos_orog)
orog_out => fields_out(pos_orog)
CALL rcf_alloc_field( orog_out )
CALL rcf_read_field( orog_out, hdr_out, decomp_rcf_output )

ALLOCATE (xi3_at_rho  (output_grid % loc_p_row_length *        &
                       output_grid % loc_p_rows,               &
                       0:output_grid % model_levels+1))
ALLOCATE (xi3_at_theta(output_grid % loc_p_row_length *        &
                       output_grid % loc_p_rows,               &
                       0:output_grid % model_levels+1))

CALL rcf_generate_heights(output_grid, orog_out, ppx_atm_tall,   &
                          st_levels_model_theta, xi3_at_theta,   &
                          fields_out(pos_thetavd) % level_size)
CALL rcf_generate_heights(output_grid, orog_out, ppx_atm_tall,   &
                          st_levels_model_rho, xi3_at_rho,       &
                          fields_out(pos_dry_rho) % level_size)


! Calculate gravity on theta levels
ALLOCATE (g_theta(fields_out(pos_thetavd) % level_size,  &
                0:fields_out(pos_thetavd) % levels) )

IF (l_cartesian .OR. l_shallow .OR. l_const_grav) THEN

  DO k = 0, fields_out(pos_thetavd) % levels
    DO i = 1, fields_out(pos_thetavd) % level_size
      g_theta(i,k) = g
    END DO
  END DO

ELSE

  DO k = 0, fields_out(pos_thetavd) % levels
    DO i = 1, fields_out(pos_thetavd) % level_size
      g_theta(i,k) = g * (planet_radius/xi3_at_theta(i,k))**2
    END DO
  END DO

END IF 

! Set up interpolation weights.
! Allocate space:
IF ( .NOT. ALLOCATED( intw_w2rho) ) THEN
  ALLOCATE( intw_w2rho(output_grid % model_levels,2))
END IF
IF ( .NOT. ALLOCATED( intw_rho2w) ) THEN
  ALLOCATE( intw_rho2w(output_grid % model_levels,2))
END IF

! Define weights:
DO k = 1, output_grid % model_levels
  intw_w2rho(k,1) = ( output_grid % eta_rho_levels(k)          &
                     - output_grid % eta_theta_levels(k-1) ) / &
                    ( output_grid % eta_theta_levels(k)        &
                     - output_grid % eta_theta_levels(k-1) )
  intw_w2rho(k,2) = 1.0 - intw_w2rho(k,1)
END DO

DO k = 1, output_grid % model_levels-1
  intw_rho2w(k,1) = ( output_grid % eta_theta_levels(k)        &
                     - output_grid % eta_rho_levels(k) ) /     &
                    ( output_grid % eta_rho_levels(k+1)        &
                     - output_grid % eta_rho_levels(k) )
  intw_rho2w(k,2) = 1.0 - intw_rho2w(k,1)
END DO
intw_rho2w(output_grid % model_levels,1) = 0.5
intw_rho2w(output_grid % model_levels,2) = 0.5

! ---------------------------------------
! Begin processing the idealised options
! ---------------------------------------

SELECT CASE (problem_number)

CASE (dynamical_core)

  ! More height calculations. u is needed each time, v in most.
  ALLOCATE (xi3_at_u (output_grid % loc_u_row_length *    &
                      output_grid % loc_u_rows,           &
                      0:output_grid % model_levels+1))
  ALLOCATE (xi3_at_v (output_grid % loc_v_row_length *    &
                      output_grid % loc_v_rows,           &
                      0:output_grid % model_levels+1))

  CALL rcf_generate_heights(output_grid, orog_out, ppx_atm_cuall, &
                            st_levels_model_rho, xi3_at_u,        &
                            fields_out(pos_u) % level_size )
  CALL rcf_generate_heights(output_grid, orog_out, ppx_atm_cuall, &
                            st_levels_model_rho, xi3_at_v ,       &
                            fields_out(pos_v) % level_size )

  CALL rcf_alloc_field(fields_out(pos_exner))

  ! Mark the top exner level as unset
  k_top = fields_out(pos_exner) % levels
  DO i = 1, fields_out(pos_exner) % level_size
    fields_out(pos_exner) % data(i,k_top) = rmdi
  END DO

  SELECT CASE(initial_profile)

  CASE(baro_inst)

    SELECT CASE (model_type)

    CASE (mt_global)
      CALL rcf_ideal_baroclinic(hdr_out, xi3_at_theta, xi3_at_rho, xi3_at_u, &
                                g_theta, fields_out, field_count_out)
      initialised = .TRUE.

    CASE (mt_cyclic_lam)
      CALL rcf_ideal_baroclinic_channel(hdr_out, xi3_at_theta, xi3_at_rho,   &
                                xi3_at_u, g_theta, fields_out, field_count_out)
      initialised = .TRUE.
    CASE DEFAULT
      icode = 80
      cmessage = 'Baroclinic profile initialisation (initial_profile = 1) '  &
              // 'is only supported by global and E-W cyclic LAM models.'
      CALL ereport(RoutineName, icode, cmessage)
    END SELECT

  CASE (rot_solid_body)

    CALL rcf_ideal_rotating_solid_body(hdr_out, xi3_at_theta, xi3_at_rho,    &
                          xi3_at_u, xi3_at_v, fields_out, field_count_out)
    initialised = .TRUE.

  CASE (deep_baro_inst)

    CALL rcf_ideal_deep_baroclinic(hdr_out, xi3_at_theta, xi3_at_rho,        &
                      xi3_at_u, xi3_at_v, g_theta, fields_out, field_count_out)
    initialised = .TRUE.

  CASE (no_preset)

    IF (tprofile_number /= imdi) THEN

      CALL rcf_ideal_generic_temperature_profile(hdr_out, xi3_at_theta,      &
                        xi3_at_rho, g_theta, fields_out, field_count_out)
      initialised = .TRUE.

    END IF

  END SELECT

  IF (idl_bubble_option(1) > 0) THEN
    CALL rcf_ideal_bubble(hdr_out, xi3_at_theta, fields_out, field_count_out)
  END IF

  ! exner data array may have been written and deallocated, in which case we
  ! need to reallocate space and read in the current contents.
  ! (If it hasn't been deallocated then we should keep the current values and
  ! not read the contents of the provisional output dump.)
  IF (.NOT. ALLOCATED(fields_out(pos_exner) % data)) THEN
    CALL rcf_alloc_field(fields_out(pos_exner))
    CALL rcf_read_field(fields_out(pos_exner), hdr_out, decomp_rcf_output)
  END IF

  ! Populate the top exner level if it's still undefined.
  ! N.B. pos_exner defined over 1:model_levels+1.
  !      xi3_at_theta and xi3_at_rho defined over 0:model_levels+1.
  IF (fields_out(pos_exner) % data(1,k_top) == rmdi) THEN
    k = output_grid % model_levels
    DO i = 1, output_grid % loc_p_field
      p1 = (xi3_at_theta(i,k) - xi3_at_rho(i,k-1))                           &
         / (xi3_at_rho(i,k)   - xi3_at_rho(i,k-1))
      p2 = (xi3_at_rho(i,k) - xi3_at_theta(i,k))                             &
         / (xi3_at_rho(i,k) - xi3_at_rho(i,k-1))
      fields_out(pos_exner) % data(i,k_top) =                                &
                 2.0 * fields_out(pos_exner) % data(i,k)**p1                 &
                     * fields_out(pos_exner) % data(i,k-1)**p2               &
                     - fields_out(pos_exner) % data(i,k)
    END DO
    CALL rcf_write_field (fields_out(pos_exner), hdr_out, decomp_rcf_output)
  END IF
  CALL rcf_dealloc_field (fields_out(pos_exner))

  DEALLOCATE(xi3_at_v)
  DEALLOCATE(xi3_at_u)

  IF (.NOT. initialised) THEN
    ErrorStatus = 100
    CALL ereport ('rcf_ideal_initialisation', ErrorStatus,                   &
    'Dynamical core-type idealised initialisation requested '                &
    // 'but no valid setup procedure specified.')
  END IF


CASE (idealised_problem)

  ! For non-zero orography, profiles will have to be set by looping
  ! over the relevant horizontal grids for scalars, u and v.
  IF (surface_type == surface_zero) THEN

    CALL rcf_alloc_field( fields_out(pos_thetavd) )
    CALL rcf_read_field( fields_out(pos_thetavd), hdr_out, decomp_rcf_output )

    ALLOCATE(z_theta_1d(0:output_grid % model_levels))
    ALLOCATE(g_theta_1d(0:output_grid % model_levels))
    ALLOCATE(theta_1d(0:output_grid % model_levels))
    ALLOCATE(thetavd_1d(0:output_grid % model_levels))
    ALLOCATE(z_rho_1d(1:output_grid % model_levels))

! Note potential mismatches, e.g. g_theta(:,0:fields_out(pos_thetavd) % levels)
    DO k = 0, output_grid % model_levels
      z_theta_1d(k) = xi3_at_theta(1,k) - planet_radius
      g_theta_1d(1) = g_theta(1,k)
    END DO
    DO k = 1, output_grid % model_levels
      z_rho_1d(k) = xi3_at_rho(1,k) - planet_radius
    END DO

! Adjust pressure over orography using isothermal assumption
! Note: if support for non-zero orography is added this will need to be
!       evaluated at each horizontal grid-point
    CALL rcf_alloc_field( fields_out(pos_exner_surf) )
    p_on_orog = p_surface * EXP(-z_theta_1d(0) * g / (r * t_surface) )
    DO i = 1, output_grid % loc_p_field
      fields_out(pos_exner_surf) % data(i,1) = (p_on_orog / p_zero)**kappa
    END DO

! Set up the required profiles. The equivalent process in the atmosphere model
! requires more profiles for use during the model run. Here, only these are
! needed:
    ! Needed by rcf_thermodynamic_profile
    CALL copy_profile(num_init_times, num_theta_init_heights,                &
                      init_time, theta_init_height, theta_init_data,         &
                      theta_init, field_type = theta_init_field_type)

    CALL copy_profile(num_init_times, num_uv_init_heights,                   &
                      init_time, uv_init_height, u_init_data, u_init)

    CALL copy_profile(num_init_times, num_uv_init_heights,                   &
                      init_time, uv_init_height, v_init_data, v_init)


    ALLOCATE(u_1d(1:output_grid % model_levels))
    ALLOCATE(v_1d(1:output_grid % model_levels))

    CALL prof_interp(interp_linear, extrap_constant,                         &
                     u_init % n_heights, u_init % height, u_init % vprof,    &
                     output_grid % model_levels, z_rho_1d, u_1d)

    CALL prof_interp(interp_linear, extrap_constant,                         &
                     v_init % n_heights, v_init % height, v_init % vprof,    &
                     output_grid % model_levels, z_rho_1d, v_1d)


    ALLOCATE(exner_1d(0:output_grid % model_levels+1))

    CALL rcf_ideal_thermodynamic_profile(output_grid % model_levels,         &
                                         p_on_orog, z_theta_1d, z_rho_1d,    &
                                         g_theta_1d, theta_1d, exner_1d)

! Hydrostatic balance assuming dry atmosphere
    ALLOCATE(mv_1d(0:output_grid % model_levels))
    ALLOCATE(dry_rho_1d(1:output_grid % model_levels))
    DO k = 0, output_grid % model_levels
      mv_1d(k) = 0.0
    END DO
    IF (pressure_balance == no_balance) THEN
      ! thetavd_1d(:) = theta_1d(:)
      ! It's not clear whether this option is still used, but a stub for
      ! it is being retained until there's certainty it can be removed.
      errorstatus = 200
      cmessage = "The mystery is of wet or dry. And where does the " //     &
                 "solution lie? A suitable pressure_balance is required."
      CALL ereport(RoutineName, errorstatus, cmessage)
    ELSE
      CALL rcf_ideal_hydrostatic_balance(output_grid % model_levels,        &
                                         p_surface, z_theta_1d, z_rho_1d,   &
                                         g_theta, theta_1d, mv_1d,          &
                                         thetavd_1d, exner_1d, dry_rho_1d)
    END IF

! Set moisture profile, if required and reset hydrostatic balance.
! NB. Need to look at combining temperature and moisture initialisation.
!     If specifiying relative humidity, need to get temperature, pressure
!     and moisture simultaneously.

    IF (.NOT. l_dry) THEN
      CALL copy_profile(num_init_times, num_mv_init_heights,                 &
                        init_time, mv_init_height, mv_init_data,             &
                        mv_init, field_type = mv_init_field_type)

      CALL rcf_ideal_vapour_profile(output_grid % model_levels, z_theta_1d,   &
                                    theta_1d, exner_1d, mv_1d)

      IF (pressure_balance == moist_hydro) THEN
        CALL rcf_ideal_hydrostatic_balance(output_grid % model_levels,        &
                                           p_surface, z_theta_1d, z_rho_1d,   &
                                           g_theta, theta_1d, mv_1d,          &
                                           thetavd_1d, exner_1d, dry_rho_1d)
      END IF
    END IF



    CALL rcf_alloc_field( fields_out(pos_exner) )
    CALL rcf_alloc_field( fields_out(pos_u) )
    CALL rcf_alloc_field( fields_out(pos_v) )
    CALL rcf_alloc_field( fields_out(pos_dry_rho) )
    CALL rcf_alloc_field( fields_out(pos_mv) )

    ! Transfer 1d values back to primary data arrays:
    DO k = 1, output_grid % model_levels + 1
      DO i = 1, output_grid % loc_p_field
        fields_out(pos_exner) % data(i,k) = exner_1d(k)
        ! 1d arrays over 0:model_levels, data array over 1:model_levels+1
        fields_out(pos_thetavd) % data(i,k) = thetavd_1d(k-1)
        fields_out(pos_mv) % data(i,k) = mv_1d(k-1)
      END DO
    END DO
  
    DO k = 1, output_grid % model_levels
  
      DO i = 1, fields_out(pos_u) % level_size
        fields_out(pos_u) % data(i,k) = u_1d(k)
      END DO
  
      DO i = 1, fields_out(pos_v) % level_size
        fields_out(pos_v) % data(i,k) = v_1d(k)
      END DO
  
      DO i = 1, fields_out(pos_dry_rho) % level_size
        fields_out(pos_dry_rho) % data(i,k) = dry_rho_1d(k)
      END DO
  
    END DO

   
    ! Print the current setup:
    IF (mype == 0 .AND. PrintStatus >= PrStatus_Diag) THEN
      WRITE(umMessage, '(A)') 'Idealised problem 1d profile output:'
      CALL umPrint(umMessage, src=RoutineName)
 
      WRITE(umMessage, '(A,A)')                                               &
      'k    z_rho(k)      exner(1,k)       dry rho(1,k)     U(1,k)',          &
      '          V(1,k)'
      CALL umPrint(umMessage, src=RoutineName)
      DO k = 1, output_grid % model_levels
        WRITE(umMessage, '(I3,F9.1,4(1X,E16.8))') k, z_rho_1d(k), &
        fields_out(pos_exner) % data(1,k), fields_out(pos_dry_rho) % data(1,k),&
        fields_out(pos_u) % data(1,k), fields_out(pos_v) % data(1,k)
        CALL umPrint(umMessage, src=RoutineName)
      END DO
  
      WRITE(umMessage, '(A,A)') 'k  z_theta(k)     theta_1d(k)   ',        &
                                ' thetavd(1,k)      mv(1,k) '
      CALL umPrint(umMessage, src=RoutineName)
      DO k = 1, output_grid % model_levels
        WRITE(umMessage, '(I3,F10.2,3(1X,E16.8))') k, z_theta_1d(k),         &
          theta_1d(k),fields_out(pos_thetavd) % data(1,k),                   &
          fields_out(pos_mv) % data(1,k)
        CALL umPrint(umMessage, src=RoutineName)
      END DO
  
    END IF

    DEALLOCATE(dry_rho_1d)
    DEALLOCATE(mv_1d)
    DEALLOCATE(exner_1d)
    DEALLOCATE(v_1d)
    DEALLOCATE(u_1d)
    DEALLOCATE(z_rho_1d)
    DEALLOCATE(thetavd_1d)
    DEALLOCATE(theta_1d)
    DEALLOCATE(g_theta_1d)
    DEALLOCATE(z_theta_1d)

    ! Write data back to output dump:
    CALL rcf_write_field  (fields_out(pos_thetavd), hdr_out, decomp_rcf_output)
    CALL rcf_dealloc_field( fields_out(pos_thetavd) )
    CALL rcf_write_field  (fields_out(pos_exner),   hdr_out, decomp_rcf_output)
    CALL rcf_dealloc_field( fields_out(pos_exner) )
    CALL rcf_write_field  (fields_out(pos_u),       hdr_out, decomp_rcf_output)
    CALL rcf_dealloc_field( fields_out(pos_u) )
    CALL rcf_write_field  (fields_out(pos_v),       hdr_out, decomp_rcf_output)
    CALL rcf_dealloc_field( fields_out(pos_v) )
    CALL rcf_write_field  (fields_out(pos_dry_rho), hdr_out, decomp_rcf_output)
    CALL rcf_dealloc_field( fields_out(pos_dry_rho) )
    CALL rcf_write_field  (fields_out(pos_mv),      hdr_out, decomp_rcf_output)
    CALL rcf_dealloc_field( fields_out(pos_mv) )
    CALL rcf_write_field  (fields_out(pos_exner_surf),hdr_out,decomp_rcf_output)
    CALL rcf_dealloc_field( fields_out(pos_exner_surf) )

  ELSE   ! Non-zero orography

    errorstatus = 300
    cmessage = 'Idealised problem initialisation only supports zero orography.'
    CALL ereport(RoutineName, errorstatus, cmessage)
  END IF  ! surface type


  ! Apply random perturbations:
  IF (l_perturb_t) THEN
    CALL rcf_ideal_random_perturbation(hdr_out, fields_out, field_count_out, &
                                stashcode_thetavd, perturb_magnitude_t, .TRUE.)
  END IF
  IF (l_perturb_q) THEN
    CALL rcf_ideal_random_perturbation(hdr_out, fields_out, field_count_out, &
                                  stashcode_mv, perturb_magnitude_q, .TRUE.)
  END IF

  ! Apply bubble perturbations:
  IF (ANY(idl_bubble_option > 0)) THEN
    CALL rcf_ideal_bubble(hdr_out, xi3_at_theta, fields_out, field_count_out)
  END IF


CASE (idealised_planet)

  CALL rcf_alloc_field(fields_out(pos_exner))

  ! Mark the top exner level as unset
  k_top = fields_out(pos_exner) % levels
  DO i = 1, fields_out(pos_exner) % level_size
    fields_out(pos_exner) % data(i,k_top) = rmdi
  END DO

  IF (tprofile_number /= imdi) THEN

    CALL rcf_ideal_isothermal_profile(hdr_out, xi3_at_theta, xi3_at_rho,     &
                                      g_theta, fields_out, field_count_out,  &
                                      .TRUE.)
    initialised = .TRUE.

  ELSE

    ErrorStatus = 300
    WRITE(cmessage, '(A,I0)') 'Unrecognised value of tprofile_number: ',    &
                              tprofile_number
    CALL ereport('rcf_ideal_initialisation', ErrorStatus, cmessage)

  END IF 

  ! exner data array may have been written and deallocated, in which case we
  ! need to reallocate space and read in the current contents.
  ! (If it hasn't been deallocated then we should keep the current values and
  ! not read the contents of the provisional output dump.)
  IF (.NOT. ALLOCATED(fields_out(pos_exner) % data)) THEN
    CALL rcf_alloc_field(fields_out(pos_exner))
    CALL rcf_read_field(fields_out(pos_exner), hdr_out, decomp_rcf_output)
  END IF

  ! Populate the top exner level if it's still undefined.
  ! N.B. pos_exner defined over 1:model_levels+1.
  !      xi3_at_theta and xi3_at_rho defined over 0:model_levels+1.
  IF (fields_out(pos_exner) % data(1,k_top) == rmdi) THEN
    k = output_grid % model_levels
    DO i = 1, output_grid % loc_p_field
      p1 = (xi3_at_theta(i,k) - xi3_at_rho(i,k-1))                           &
         / (xi3_at_rho(i,k)   - xi3_at_rho(i,k-1))
      p2 = (xi3_at_rho(i,k) - xi3_at_theta(i,k))                             &
         / (xi3_at_rho(i,k) - xi3_at_rho(i,k-1))
      fields_out(pos_exner) % data(i,k_top) =                                &
                 2.0 * fields_out(pos_exner) % data(i,k)**p1                 &
                     * fields_out(pos_exner) % data(i,k-1)**p2               &
                     - fields_out(pos_exner) % data(i,k)
    END DO
    CALL rcf_write_field (fields_out(pos_exner), hdr_out, decomp_rcf_output)
  END IF

  ! Calculate vapour mixing ratio
  CALL rcf_ideal_vapour_calc( hdr_out, xi3_at_theta, fields_out,             &
                              field_count_out )

  CALL rcf_dealloc_field (fields_out(pos_exner))

  IF (.NOT. initialised) THEN
    ErrorStatus = 310
    CALL ereport ('rcf_ideal_initialisation', ErrorStatus,                   &
    'Idealised planet-type idealised initialisation requested '              &
    // 'but no valid setup procedure specified.')
  END IF

END SELECT

DEALLOCATE(intw_rho2w)
DEALLOCATE(intw_w2rho)

! Recalculate (wet) theta and rho, which still potentially contain garbage
! from the initial dump:
CALL rcf_alloc_field( fields_out(pos_theta) )
CALL rcf_alloc_field( fields_out(pos_mv) )
CALL rcf_read_field( fields_out(pos_mv),      hdr_out, decomp_rcf_output )
CALL rcf_alloc_field( fields_out(pos_thetavd) )
CALL rcf_read_field( fields_out(pos_thetavd), hdr_out, decomp_rcf_output )

DO k = 1, output_grid % model_levels + 1
  DO i = 1, output_grid % loc_p_field
    fields_out(pos_theta) % data(i,k) = fields_out(pos_thetavd) % data(i,k)   &
             / (1.0 + fields_out(pos_mv) % data(i,k) * recip_epsilon)
  END DO
END DO

CALL rcf_write_field  (fields_out(pos_theta), hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field( fields_out(pos_theta) )
CALL rcf_dealloc_field( fields_out(pos_thetavd) )
CALL rcf_dealloc_field( fields_out(pos_mv) )

CALL rcf_recompute_wet_rho( fields_out, field_count_out, hdr_out)

DEALLOCATE(xi3_at_theta)
DEALLOCATE(xi3_at_rho)

IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE  rcf_ideal_initialisation
END MODULE rcf_ideal_initialisation_mod
