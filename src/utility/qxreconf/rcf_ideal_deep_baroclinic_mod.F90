! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Idealised deep baroclinic setup

MODULE rcf_ideal_deep_baroclinic_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_DEEP_BAROCLINIC_MOD'

CONTAINS

! Subroutine rcf_ideal_deep_baroclinic
!
! Description:
!   Sets wind field, pressure, potential temperature and density
!   for baroclinic wave solution of deep atmosphere equations.
!   Isothermal case only. Version based on Staniforth deep tests.
!
! Method:
!   Simply set fields using analytic formulae. May need to balance
!   pressure field with respect to model discretisation?
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_deep_baroclinic ( hdr_out, xi3_at_theta, xi3_at_rho, &
                     xi3_at_u, xi3_at_v,  g_theta, fields_out, field_count_out)

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type
                                 
USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_Locate_Mod, ONLY: &
    rcf_locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    rcf_alloc_field,           &
    rcf_dealloc_field

USE Rcf_Write_Field_Mod, ONLY: &
    rcf_write_field

USE rcf_calc_coords_mod, ONLY: &
    rcf_calc_coords

USE rcf_scatter_field_mod, ONLY: &
    rcf_scatter_field_real

USE rcf_ideal_1d_profiles_mod, ONLY: &
    rcf_ideal_1d_profiles

USE rcf_interp_weights_mod, ONLY: &
    intw_w2rho

USE decomp_params, ONLY: &
    decomp_rcf_output

USE um_parcore, ONLY: &
    mype,             &
    nproc

USE um_parvars, ONLY: &
    gc_all_proc_group

USE rcf_nlist_recon_idealised_mod, ONLY: &
    b_const,                &
    grid_np_lat,            &
    grid_np_lon,            &
    k_const,                &
    l_baro_perturbed,       &
    l_rotate_grid,          &
    t0_e,                   &
    t0_p

USE idealise_run_mod, ONLY: &
    l_shallow

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,     &
    stashcode_exner,        &
    stashcode_exner_surf,   &
    stashcode_dry_rho,      &
    stashcode_thetavd,      &
    stashcode_u,            &
    stashcode_v

USE rcf_ideal_deep_baroclinic_constants_mod, ONLY: &
    b,                          &
    c,                          &
    eg_lapse,                   &
    h,                          &
    t0,                         &
    taper_top,                  &
    tau1,                       &
    tau2,                       &
    tau1_int,                   &
    tau2_int

USE rcf_ideal_deep_baro_pert_mod, ONLY: &
    deep_baro_pert

USE planet_constants_mod, ONLY: &
    g,                          &
    r,                          &
    omega,                      &
    two_omega,                  &
    kappa,                      &
    planet_radius,              &
    p_zero

USE conversions_mod, ONLY: &
    pi,                    &
    pi_over_180

USE rcf_ideal_tprofile_mod, ONLY: &
    tp_BV_isoth

USE model_domain_mod, ONLY: &
    model_type,             &
    mt_global

USE lam_config_inputs_mod, ONLY: &
    delta_lon

USE rcf_eg_poles_mod, ONLY: &
    rcf_eg_poles

USE umPrintMgr, ONLY:       &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE(um_header_type), INTENT(IN) :: hdr_out      ! Output dump header
REAL, INTENT(IN)          :: xi3_at_theta(output_grid % loc_p_field,         &
                                        0:output_grid % model_levels+1)
                                          ! heights at theta levels
REAL, INTENT(IN)          :: xi3_at_rho(output_grid % loc_p_field,           &
                                      0:output_grid % model_levels+1)
                                          ! heights at rho levels
REAL, INTENT(IN)          :: xi3_at_u(output_grid % loc_u_field,             &
                                    0:output_grid % model_levels+1)
                                          ! heights at u levels
REAL, INTENT(IN)          :: xi3_at_v(output_grid % loc_v_field,             &
                                    0:output_grid % model_levels+1)
                                          ! heights at v levels
REAL, INTENT(IN)          :: g_theta(:,:)        ! gravity at theta_levels
TYPE(field_type), POINTER :: fields_out(:)       ! output fields
INTEGER, INTENT(IN)       :: field_count_out     ! no. of output fields

! Local variables
INTEGER :: i
INTEGER :: j
INTEGER :: k

REAL    :: t_guess   ! First guess at temp. sol'n.
REAL    :: p_guess   ! First guess at pressure sol'n.
REAL    :: xi1_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi2_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi1_u_global(output_grid % glob_u_row_length, &
                        output_grid % glob_u_rows)
REAL    :: xi2_v_global(output_grid % glob_v_row_length, &
                        output_grid % glob_v_rows)
REAL    :: xi1_p(output_grid % loc_p_row_length,         &
                 output_grid % loc_p_rows)
REAL    :: xi2_p(output_grid % loc_p_row_length,         &
                 output_grid % loc_p_rows)
REAL    :: xi1_u(output_grid % loc_u_row_length,         &
                 output_grid % loc_u_rows)
REAL    :: xi2_v(output_grid % loc_v_row_length,         &
                 output_grid % loc_v_rows)
REAL    :: xi1_p_1d(output_grid % loc_p_field)
REAL    :: xi2_p_1d(output_grid % loc_p_field)
REAL    :: xi1_u_1d(output_grid % loc_u_field)
REAL    :: xi2_v_1d(output_grid % loc_v_field)
REAL    :: u_out_2d(output_grid % loc_u_row_length, output_grid % loc_u_rows)
REAL    :: v_out_2d(output_grid % loc_v_row_length, output_grid % loc_v_rows)
REAL    :: exner_ext(0:output_grid % model_levels+1) ! Extended 1D exner profile

! Indices for locating fields
INTEGER  :: pos_thetavd
INTEGER  :: pos_dry_rho
INTEGER  :: pos_exner
INTEGER  :: pos_exner_surf
INTEGER  :: pos_u
INTEGER  :: pos_v
INTEGER  :: pe

INTEGER  :: i_err

! Pointers to output fields:
TYPE( field_type ), POINTER        :: exner_out
TYPE( field_type ), POINTER        :: thetavd_out
TYPE( field_type ), POINTER        :: dry_rho_out
TYPE( field_type ), POINTER        :: exner_surf_out
TYPE( field_type ), POINTER        :: u_out
TYPE( field_type ), POINTER        :: v_out

! Variables for the rotated grid:
REAL :: lon_p
REAL :: lat_p
REAL :: sin_lat_p
REAL :: cos_lat_p
REAL :: xi1_rot(output_grid % loc_p_field)
REAL :: xi2_rot(output_grid % loc_p_field)
REAL :: xx
REAL :: yy
REAL :: sn_lonmlonp_u(output_grid % loc_u_field)
REAL :: cs_lonmlonp_u(output_grid % loc_u_field)
REAL :: sn_rot_u(output_grid % loc_u_field)
REAL :: cs_rot_u(output_grid % loc_u_field)
REAL :: sn_lonmlonp_v(output_grid % loc_v_field)
REAL :: cs_lonmlonp_v(output_grid % loc_v_field)
REAL :: sn_rot_v(output_grid % loc_v_field)
REAL :: cs_rot_v(output_grid % loc_v_field)

REAL :: r_height
REAL :: taper
REAL :: u_at_v
REAL :: v_at_u
REAL :: tau1_u
REAL :: tau2_u
REAL :: tau2_int_u
REAL :: tau2_term
REAL :: t_temp

REAL :: divergence
REAL :: div_max
REAL :: div_min
REAL :: vorticity
REAL :: vort_max
REAL :: vort_min

INTEGER, PARAMETER :: pert_type = 3  ! Fixed choice of perturbation option

! Values required for rcf_ideal_1d_profs
REAL,    PARAMETER :: dtheta_dz1(3) = (/0.0, 0.0, 0.0/)
REAL,    PARAMETER :: height_dz1(2) = (/0.0, 0.0/)
LOGICAL, PARAMETER :: l_constant_dz = .TRUE.

CHARACTER (LEN=*),  PARAMETER :: routinename='RCF_IDEAL_DEEP_BAROCLINIC'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Extract some fields
CALL rcf_locate( stashcode_prog_sec, stashcode_thetavd,         &
               fields_out, field_count_out, pos_thetavd)
thetavd_out => fields_out(pos_thetavd)
CALL rcf_alloc_field( thetavd_out )

CALL rcf_locate( stashcode_prog_sec, stashcode_dry_rho,         &
               fields_out, field_count_out, pos_dry_rho)
dry_rho_out => fields_out(pos_dry_rho)
CALL rcf_alloc_field( dry_rho_out )

CALL rcf_locate( stashcode_prog_sec, stashcode_exner,           &
               fields_out, field_count_out, pos_exner)
exner_out => fields_out(pos_exner)
CALL rcf_alloc_field( exner_out )

CALL rcf_locate( stashcode_prog_sec, stashcode_exner_surf,      &
               fields_out, field_count_out, pos_exner_surf)
exner_surf_out => fields_out(pos_exner_surf)
CALL rcf_alloc_field( exner_surf_out )

CALL rcf_locate( stashcode_prog_sec, stashcode_u,               &
               fields_out, field_count_out, pos_u)
u_out => fields_out(pos_u)
CALL rcf_alloc_field( u_out )

CALL rcf_locate( stashcode_prog_sec, stashcode_v,               &
               fields_out, field_count_out, pos_v)
v_out => fields_out(pos_v)
CALL rcf_alloc_field( v_out )

! Calculate global grid info and scatter to generate local data
CALL rcf_calc_coords(hdr_out, output_grid, xi1_p_global, xi2_p_global, &
                     xi1_u_global, xi2_v_global)

pe = 0

CALL rcf_scatter_field_real(xi1_p, xi1_p_global,                 &
                            output_grid % loc_p_row_length,      &
                            output_grid % loc_p_rows,            &
                            output_grid % glob_p_row_length,     &
                            output_grid % glob_p_rows, pe,       &
                            gc_all_proc_group )

CALL rcf_scatter_field_real(xi2_p, xi2_p_global,                 &
                            output_grid % loc_p_row_length,      &
                            output_grid % loc_p_rows,            &
                            output_grid % glob_p_row_length,     &
                            output_grid % glob_p_rows, pe,       &
                            gc_all_proc_group )

CALL rcf_scatter_field_real(xi1_u, xi1_u_global,                 &
                            output_grid % loc_u_row_length,      &
                            output_grid % loc_u_rows,            &
                            output_grid % glob_u_row_length,     &
                            output_grid % glob_u_rows, pe,       &
                            gc_all_proc_group )

CALL rcf_scatter_field_real(xi2_v, xi2_v_global,                 &
                            output_grid % loc_v_row_length,      &
                            output_grid % loc_v_rows,            &
                            output_grid % glob_v_row_length,     &
                            output_grid % glob_v_rows, pe,       &
                            gc_all_proc_group )


! Idealised code expects 1D arrays.
xi1_p_1d = RESHAPE(xi1_p, (/output_grid % loc_p_field/))
xi2_p_1d = RESHAPE(xi2_p, (/output_grid % loc_p_field/))
xi1_u_1d = RESHAPE(xi1_u, (/output_grid % loc_u_field/))
xi2_v_1d = RESHAPE(xi2_v, (/output_grid % loc_v_field/))


! Begin deep baroclinic test case
!
! Set up some constants:
t0 = 0.5 * (t0_e + t0_p)
b = (t0 - t0_p) / (t0 * t0_p)
c = (k_const + 2.0) / 2.0 * (t0_e - t0_p) / (t0_e * t0_p)
h = r * t0/g

IF (l_rotate_grid) THEN
  lon_p = grid_np_lon * pi_over_180
  lat_p = grid_np_lat * pi_over_180
  sin_lat_p = SIN(lat_p)
  cos_lat_p = COS(lat_p)
END IF

IF (.NOT. ALLOCATED(tau1)) THEN
  ALLOCATE (tau1(0:output_grid % model_levels))
END IF
IF (.NOT. ALLOCATED(tau2)) THEN
  ALLOCATE (tau2(0:output_grid % model_levels))
END IF
IF (.NOT. ALLOCATED(tau1_int)) THEN
  ALLOCATE (tau1_int(0:output_grid % model_levels))
END IF
IF (.NOT. ALLOCATED(tau2_int)) THEN
  ALLOCATE (tau2_int(0:output_grid % model_levels))
END IF

! Put each column into hydrostatic balance with temperature of
! analytic steady state.
t_guess = t0
p_guess = 1.0e5

DO i = 1, output_grid % loc_p_field
  IF (l_rotate_grid) THEN
    xi2_rot(i) = ASIN(SIN(xi2_p_1d(i)) * sin_lat_p                           &
                    + COS(xi2_p_1d(i)) * COS(xi1_p_1d(i)) * cos_lat_p)
  ELSE
    xi2_rot(i) = xi2_p_1d(i)
  END IF

  ! Set up tau fields
  DO k = 0, output_grid % model_levels
    IF (k == 0) THEN
      r_height = planet_radius
    ELSE
      r_height = xi3_at_theta(i,k)
    END IF

    CALL calculate_tau_terms(r_height, tau2(k), tau1(k), tau2_int(k),        &
                             tau1_int(k))
  END DO

  CALL rcf_ideal_1d_profiles( exner_ext, dry_rho_out % data(i,:),            &
                              thetavd_out % data(i,:), tp_BV_isoth,          &
                              t_guess, p_guess, xi2_rot(i),                  &
                              xi3_at_theta(i,:), xi3_at_rho(i,:),            &
                              output_grid % z_top_of_model, planet_radius,   &
                              dtheta_dz1, height_dz1, l_constant_dz,         &
                              g_theta(i,:), output_grid % model_levels)

!  Copy relevant part of extended exner array back to the actual field
  DO k = 1, output_grid % model_levels+1
    exner_out % data(i, k) = exner_ext(k)
  END DO
  exner_surf_out % data(i, 1) = exner_ext(0)

END DO

! Balanced pressure field
IF (l_baro_perturbed) THEN

! Note repeated potential array bounds violations, as everything is calculated
! over loc_p_field. The same logic was contained in the original atmosphere
! code this was derived from.
! Set up xi1 metrics. xi2_rot was already calculated earlier.
  IF (l_rotate_grid) THEN
    DO i = 1, output_grid % loc_p_field
      xx = -sin_lat_p * COS(xi2_p_1d(i)) * SIN(xi1_p_1d(i))
      yy = SIN(xi2_rot(i)) * cos_lat_p - COS(xi2_p_1d(i)) * COS(xi1_p_1d(i))
      IF (ABS(yy) < EPSILON(1.0)) THEN
        IF (ABS(xx) < EPSILON(1.0)) THEN
          xi1_rot(i) = lon_p + pi
        ELSE
          xi1_rot(i) = lon_p + 2.0*pi
        END IF
      ELSE
        xi1_rot(i) = lon_p + ATAN2(xx,yy) + pi
      END IF
    END DO
  ELSE
    ! Note recon u grid is offset by 1 relative to its position in the model,
    ! and the point on the final boundary is not stored.
    DO i = 1, output_grid % loc_p_field - 1
      xi1_rot(i) = xi1_u_1d(i+1)
    END DO
    i = output_grid % loc_p_field
    IF (model_type == mt_global) THEN
      xi1_rot(i) = xi1_u_1d(1)
    ELSE
      xi1_rot(i) = xi1_u_1d(i-1) + delta_lon
    END IF
  END IF

  DO k = 1, output_grid % model_levels
    DO i = 1, output_grid % loc_p_field

      CALL calculate_tau_terms(xi3_at_rho(i,k), tau2_u, tau1_u)

      IF (l_shallow) THEN
        r_height = planet_radius
      ELSE
        r_height = xi3_at_rho(i,k)
      END IF

      ! Unperturbed wind field
      tau2_term = (r_height / planet_radius * COS(xi2_rot(i)))**k_const      &
                - (k_const / (k_const + 2.0)) *                              &
                   (r_height / planet_radius * COS(xi2_rot(i)))**(k_const + 2.0)

      t_temp = (planet_radius / r_height)**2 / (tau1_u - tau2_u * tau2_term)

      IF (xi3_at_rho(i,k) - planet_radius <= taper_top) THEN
        taper = 1.0                                                          &
              - 3.0 * ((xi3_at_rho(i,k) - planet_radius) / taper_top)**2     &
              + 2.0 * ((xi3_at_rho(i,k) - planet_radius) / taper_top)**3

        exner_out % data(i,k) = exner_out % data(i,k) * EXP(-taper           &
            * deep_baro_pert(xi1_rot(i), xi2_rot(i), pert_type, 3)           &
            * two_omega * SIN(xi2_rot(i)) / (r * t_temp) )

        ! thetavd_out declared over 1:model_levels+1, not 0:model_levels, so +1:
        dry_rho_out % data(i,k) = p_zero                                     &
               / ( r * (intw_w2rho(k,1) * thetavd_out % data(i,k+1)          &
                    + intw_w2rho(k,2) * thetavd_out % data(i,k)) )           &
               * exner_out % data(i,k)**((1.0 - kappa) / kappa)
      END IF
    END DO
  END DO
END IF ! IF (l_baro_perturbed)

! Wind components.
! More potential bounds violations follow (looping over loc_u_field).
! Set metrics:
IF (l_rotate_grid) THEN
  DO i = 1, output_grid % loc_u_field
    xi2_rot(i) = ASIN(SIN(xi2_p_1d(i)) * sin_lat_p                           &
                    + COS(xi2_p_1d(i)) * COS(xi1_u_1d(i)) * cos_lat_p)
    xx = -sin_lat_p * COS(xi2_p_1d(i)) * SIN(xi1_u_1d(i))
    yy = SIN(xi2_rot(i)) * cos_lat_p - COS(xi2_p_1d(i)) * COS(xi1_u_1d(i))
    IF (ABS(yy) < EPSILON(1.0)) THEN
      IF (ABS(xx) < EPSILON(1.0)) THEN
        xi1_rot(i) = lon_p + pi
      ELSE
        xi1_rot(i) = lon_p + 2.0*pi
      END IF
    ELSE
      xi1_rot(i) = lon_p + ATAN2(xx,yy) + pi
    END IF
  END DO
ELSE
  DO i = 1, output_grid % loc_u_field
    xi1_rot(i) = xi1_u_1d(i)
  END DO
END IF

DO k = 1, output_grid % model_levels
  DO i = 1, output_grid % loc_u_field

    CALL calculate_tau_terms(xi3_at_u(i,k), tau2_u, tau1_u, tau2_int_u)

    IF (l_shallow) THEN
      r_height = planet_radius
    ELSE
      r_height = xi3_at_u(i,k)
    END IF

    ! Unperturbed wind field
    tau2_term = (r_height / planet_radius * COS(xi2_rot(i)))**k_const        &
              - (k_const / (k_const + 2.0))                                  &
                * (r_height / planet_radius * COS(xi2_rot(i)))**(k_const + 2.0)

    t_temp = (planet_radius / r_height)**2 / (tau1_u - tau2_u * tau2_term)

    ! Wind proxy:
    u_out % data(i,k) = (g / planet_radius) * k_const * tau2_int_u           &
           * ( (r_height / planet_radius * COS(xi2_rot(i)))**(k_const - 1)   &
             - (r_height / planet_radius * COS(xi2_rot(i)))**(k_const + 1) ) &
           * t_temp

    u_out % data(i,k) = - omega * r_height * COS(xi2_rot(i))                 &
                        + SQRT( (omega * r_height * COS(xi2_rot(i)))**2      &
                              + r_height * COS(xi2_rot(i)) * u_out % data(i,k))

    IF (l_baro_perturbed) THEN
      ! Add perturbation
      IF (xi3_at_u(i,k) - planet_radius <= taper_top) THEN
        taper = 1.0                                                          &
              - 3.0 * ((xi3_at_u(i,k) - planet_radius) / taper_top)**2       &
              + 2.0 * ((xi3_at_u(i,k) - planet_radius) / taper_top)**3

        u_out % data(i,k) = u_out % data(i,k) + taper                        &
                * deep_baro_pert(xi1_rot(i), xi2_rot(i), pert_type, 1)

      END IF
    END IF

  END DO
END DO

DO k = 1, output_grid % model_levels
  DO i = 1, output_grid % loc_v_field
    v_out % data(i,k) = 0.0
  END DO
END DO

IF (l_baro_perturbed .AND. pert_type /= 1) THEN

! Yet more potential bounds violations follow (looping over loc_v_field).
! Set metrics:
  IF (l_rotate_grid) THEN
    DO i = 1, output_grid % loc_v_field
      xi2_rot(i) = ASIN(SIN(xi2_v_1d(i)) * sin_lat_p                         &
                      + COS(xi2_v_1d(i)) * COS(xi1_p_1d(i)) * cos_lat_p)
      xx = -sin_lat_p * COS(xi2_v_1d(i)) * SIN(xi1_p_1d(i))
      yy = SIN(xi2_rot(i)) * cos_lat_p - COS(xi2_v_1d(i)) * COS(xi1_p_1d(i))
      IF (ABS(yy) < EPSILON(1.0)) THEN
        IF (ABS(xx) < EPSILON(1.0)) THEN
          xi1_rot(i) = lon_p + pi
        ELSE
          xi1_rot(i) = lon_p + 2.0*pi
        END IF
      ELSE
        xi1_rot(i) = lon_p + ATAN2(xx,yy) + pi
      END IF
    END DO
  ELSE
    DO i = 1, output_grid % loc_v_field
      xi1_rot(i) = xi1_p_1d(i)
    END DO
  END IF

  DO k = 1, output_grid % model_levels
    DO i = 1, output_grid % loc_v_field
      IF (xi3_at_v(i,k) - planet_radius <= taper_top) THEN
        taper = 1.0                                                          &
              - 3.0 * ((xi3_at_v(i,k) - planet_radius) / taper_top)**2       &
              + 2.0 * ((xi3_at_v(i,k) - planet_radius) / taper_top)**3

        v_out % data(i,k) =                                                  &
             taper * deep_baro_pert(xi1_rot(i), xi2_rot(i), pert_type, 2)
      END IF
    END DO
  END DO

END IF


! Rotate wind vectors
IF (l_rotate_grid) THEN

  ! Metrics for u winds:
  DO i = 1, output_grid % loc_u_field
    xi2_rot(i) = ASIN(SIN(xi2_p_1d(i)) * sin_lat_p                           &
                    + COS(xi2_p_1d(i)) * COS(xi1_u_1d(i)) * cos_lat_p)

    sn_lonmlonp_u(i) = -COS(xi2_p_1d(i)) * SIN(xi1_u_1d(i)) / COS(xi2_rot(i))
    cs_lonmlonp_u(i) = (SIN(xi2_rot(i)) * cos_lat_p                          &
                     - COS(xi2_p_1d(i)) * COS(xi1_u_1d(i)))                  &
                     / (COS(xi2_rot(i)) * sin_lat_p)

    cs_rot_u(i) = - COS(xi1_u_1d(i)) * cs_lonmlonp_u(i)                      &
                  - SIN(xi1_u_1d(i)) * sn_lonmlonp_u(i) * sin_lat_p
    sn_rot_u(i) = COS(xi1_u_1d(i)) * SIN(xi2_rot(i)) * sn_lonmlonp_u(i)      &
                - SIN(xi1_u_1d(i)) * (COS(xi2_rot(i)) * cos_lat_p            &
                             + SIN(xi2_rot(i)) * sin_lat_p * cs_lonmlonp_u(i))
  END DO

  ! Add perturbation to u
  IF (l_baro_perturbed) THEN
    DO k = 1, output_grid % model_levels
      DO i = 1, output_grid % loc_u_field
        v_at_u = 0.0
        IF (xi3_at_u(i,k) - planet_radius <= taper_top) THEN
          xx = -sin_lat_p * COS(xi2_p_1d(i)) * SIN(xi1_u_1d(i))
          yy = SIN(xi2_rot(i)) * cos_lat_p - COS(xi2_p_1d(i)) * COS(xi1_u_1d(i))
          IF (ABS(yy) < EPSILON(1.0)) THEN
            IF (ABS(xx) < EPSILON(1.0)) THEN
              xi1_rot(i) = lon_p + pi
            ELSE
              xi1_rot(i) = lon_p + 2.0*pi
            END IF
          ELSE
            xi1_rot(i) = lon_p + ATAN2(xx,yy) + pi
          END IF

          taper = 1.0                                                        &
                - 3.0 * ((xi3_at_u(i,k) - planet_radius) / taper_top)**2     &
                + 2.0 * ((xi3_at_u(i,k) - planet_radius) / taper_top)**3
          v_at_u = taper                                                     &
                   * deep_baro_pert(xi1_rot(i), xi2_rot(i), pert_type, 2)
        END IF

        u_out % data(i,k) = u_out % data(i,k) * cs_rot_u(i)                  &
                          + v_at_u * sn_rot_u(i)
      END DO
    END DO
  ELSE
    DO k = 1, output_grid % model_levels
      DO i = 1, output_grid % loc_u_field
        u_out % data(i,k) = u_out % data(i,k) * cs_rot_u(i)
      END DO
    END DO
  END IF

  ! Metrics for v winds:
  DO i = 1, output_grid % loc_v_field
    xi2_rot(i) = ASIN(SIN(xi2_v_1d(i)) * sin_lat_p                           &
                    + COS(xi2_v_1d(i)) * COS(xi1_p_1d(i)) * cos_lat_p)
    sn_lonmlonp_v(i) = -COS(xi2_v_1d(i)) * SIN(xi1_p_1d(i)) / COS(xi2_rot(i))
    cs_lonmlonp_v(i) = (SIN(xi2_rot(i)) * cos_lat_p                          &
                     - COS(xi2_v_1d(i)) * COS(xi1_p_1d(i)))                  &
                     / (COS(xi2_rot(i)) * sin_lat_p)

    cs_rot_v(i) = - COS(xi1_p_1d(i)) * cs_lonmlonp_v(i)                      &
                  - SIN(xi1_p_1d(i)) * sn_lonmlonp_v(i) * sin_lat_p
    sn_rot_v(i) = COS(xi1_p_1d(i)) * SIN(xi2_rot(i)) * sn_lonmlonp_v(i)      &
                - SIN(xi1_p_1d(i)) * (COS(xi2_rot(i)) * cos_lat_p            &
                             + SIN(xi2_rot(i)) * sin_lat_p * cs_lonmlonp_v(i))
  END DO

  ! Add perturbation to v
  IF (l_baro_perturbed) THEN
    DO i = 1, output_grid % loc_v_field
      xx = -sin_lat_p * COS(xi2_v_1d(i)) * SIN(xi1_p_1d(i))
      yy = SIN(xi2_rot(i)) * cos_lat_p - COS(xi2_v_1d(i)) * COS(xi1_p_1d(i))
      IF (ABS(yy) < EPSILON(1.0)) THEN
        IF (ABS(xx) < EPSILON(1.0)) THEN
          xi1_rot(i) = lon_p + pi
        ELSE
          xi1_rot(i) = lon_p + 2.0*pi
        END IF
      ELSE
        xi1_rot(i) = lon_p + ATAN2(xx,yy) + pi
      END IF
    END DO
  END IF

  DO k = 1, output_grid % model_levels
    DO i = 1, output_grid % loc_v_field
      CALL calculate_tau_terms(xi3_at_v(i,k), tau2_u, tau1_u, tau2_int_u)

      IF (l_shallow) THEN
        r_height = planet_radius
      ELSE
        r_height = xi3_at_v(i,k)
      END IF

      ! Unperturbed wind field:
      tau2_term = (r_height / planet_radius * COS(xi2_rot(i)))**k_const      &
                - (k_const / (k_const + 2.0)) * (r_height                    &
                        / planet_radius * COS(xi2_rot(i)))**(k_const + 2.0)

      t_temp = (planet_radius / r_height)**2 / (tau1_u - tau2_u * tau2_term)

      ! Wind proxy:
      u_at_v = (g / planet_radius) * k_const * tau2_int_u                    &
              * ((r_height / planet_radius * COS(xi2_rot(i)))**(k_const-1)   &
               - (r_height / planet_radius * COS(xi2_rot(i)))**(k_const+1))  &
              * t_temp
      u_at_v = - omega * r_height * COS(xi2_rot(i))                          &
               + SQRT( (omega * r_height * COS(xi2_rot(i)))**2               &
                      + r_height * COS(xi2_rot(i)) * u_at_v)

      ! Add perturbation
      IF (l_baro_perturbed) THEN
        IF (xi3_at_v(i,k) - planet_radius <= taper_top) THEN
          taper = 1.0                                                        &
                - 3.0 * ((xi3_at_v(i,k) - planet_radius) / taper_top)**2     &
                + 2.0 * ((xi3_at_v(i,k) - planet_radius) / taper_top)**3
          u_at_v = u_at_v + taper                                            &
                    * deep_baro_pert(xi1_rot(i), xi2_rot(i), pert_type, 1)
        END IF
      END IF

      v_out % data(i,k) = - u_at_v * sn_rot_v(i)                             &
                          + v_out % data(i,k) * cs_rot_v(i)
    END DO
  END DO


  ! Set pole winds correctly
  IF (model_type == mt_global) THEN
    CALL rcf_eg_poles(fields_out, field_count_out, hdr_out)
  END IF

END IF  ! IF (l_rotate_grid)

! Print max/min initial divergence /vorticity on level 1, if required.
! Note the reconfiguration does not store the winds on the north or east
! boundaries, hence the final row and column are missing from this calculation.
IF (PrintStatus >= PrStatus_Normal) THEN
  k = 1
  div_max  = -1.0
  div_min  =  1.0
  vort_max = -1.0
  vort_min =  1.0

  ! These calculations are much easier to process in two dimensions:
  u_out_2d = RESHAPE(u_out % data(:,1), (/output_grid % loc_u_row_length,    &
                                          output_grid % loc_u_rows/))
  v_out_2d = RESHAPE(v_out % data(:,1), (/output_grid % loc_v_row_length,    &
                                          output_grid % loc_v_rows/))

  ! u and v grids start at index 1 in the reconfiguration, not 0.
  DO j = 1, output_grid % loc_u_rows - 1
    DO i = 1, output_grid % loc_u_row_length - 1
      divergence = (u_out_2d(i+1,j) - u_out_2d(i,j))                         &
                      / (xi1_u(i+1,j) - xi1_u(i,j))                          &
                 + (COS(xi2_v(i,j+1)) * v_out_2d(i,j+1)                      &
                      - COS(xi2_v(i,j)) * v_out_2d(i,j))                     &
                   / (xi2_v(i,j+1) - xi2_v(i,j))                             &
                   * 1.0 / (planet_radius * COS(xi2_p(i,j)))
      div_max = MAX(divergence, div_max)
      div_min = MIN(divergence, div_min)
    END DO
  END DO

  DO j = 1, output_grid % loc_v_rows - 1
    DO i = 1, output_grid % loc_v_row_length - 1
      vorticity = (v_out_2d(i+1,j+1) - v_out_2d(i,j+1)) /                    &
                   (xi1_p(i,j+1) - xi1_p(i,j))                               &
                + (COS(xi2_p(i,j+1)) * u_out_2d(i+1,j+1)                     &
                     - COS(xi2_p(i,j)) * u_out_2d(i+1,j))                    &
                  / (xi2_p(i,j+1) - xi2_p(i,j))
      vort_max = MAX(vorticity, vort_max)
      vort_min = MIN(vorticity, vort_min)
    END DO
  END DO

  CALL gc_rmax(1,nproc,i_err,div_max)
  CALL gc_rmin(1,nproc,i_err,div_min)
  CALL gc_rmax(1,nproc,i_err,vort_max)
  CALL gc_rmin(1,nproc,i_err,vort_min)

  IF (mype == 0) THEN
    WRITE(umMessage,'(A,2E16.8)') 'Max/Min initial divergence: ',            &
                                  div_max, div_min
    CALL umPrint(umMessage, src='rcf_ideal_deep_baroclinic')
    WRITE(umMessage,'(A,2E16.8)') 'Max/Min initial vorticity:  ',            &
                                  vort_max, vort_min
    CALL umPrint(umMessage, src='rcf_ideal_deep_baroclinic')
  END IF
END IF


! Write out and deallocate the modified fields
CALL rcf_write_field   (fields_out(pos_thetavd),    hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_thetavd))
CALL rcf_write_field   (fields_out(pos_dry_rho),    hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_dry_rho))
CALL rcf_write_field   (fields_out(pos_exner),      hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_exner))
CALL rcf_write_field   (fields_out(pos_exner_surf), hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_exner_surf))
CALL rcf_write_field   (fields_out(pos_u),          hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_u))
CALL rcf_write_field   (fields_out(pos_v),          hdr_out, decomp_rcf_output)
CALL rcf_dealloc_field (fields_out(pos_v))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN

CONTAINS

SUBROUTINE calculate_tau_terms(r_height, tau2, tau1, tau2_int, tau1_int)

IMPLICIT NONE

REAL, INTENT(IN)            :: r_height
REAL, INTENT(OUT)           :: tau2
REAL, INTENT(OUT)           :: tau1
REAL, INTENT(OUT), OPTIONAL :: tau2_int
REAL, INTENT(OUT), OPTIONAL :: tau1_int
REAL                        :: temp

CHARACTER (LEN=*),  PARAMETER :: routinename='CALCULATE_TAU_TERMS'
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

temp = ((r_height - planet_radius) / (b_const * h))**2

tau2 = c * (1.0 - 2.0*temp) * EXP(-temp)
tau1 = 1.0 / t0 * EXP(eg_lapse / t0 * (r_height - planet_radius))       &
          + b * (1.0 - 2.0*temp) * EXP(-temp)

IF (PRESENT(tau2_int)) THEN
  tau2_int = c * (r_height - planet_radius) * EXP(-temp)
END IF

IF (PRESENT(tau1_int)) THEN
  tau1_int = 1.0 / eg_lapse                                             &
            * (EXP(eg_lapse / t0 * (r_height - planet_radius)) - 1.0)   &
            + b * (r_height - planet_radius) * EXP(-temp)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE calculate_tau_terms

END SUBROUTINE rcf_ideal_deep_baroclinic
END MODULE rcf_ideal_deep_baroclinic_mod
