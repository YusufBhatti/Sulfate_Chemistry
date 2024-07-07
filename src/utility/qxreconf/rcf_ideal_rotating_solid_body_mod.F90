! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Idealised rotating solid body setup

MODULE rcf_ideal_rotating_solid_body_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_ROTATING_SOLID_BODY_MOD'

CONTAINS

! Subroutine rcf_ideal_rotating_solid_body
!
! Description:
!   Sets wind field, pressure, potential temperature and density
!   for exact solid body rotation solution of deep atmosphere equations.
!   Isothermal case only. Version based on Williamson shallow water tests.
!
! Method:
!   Simply set fields using analytic formulae. May need to balance
!   pressure field with respect to model discretisation?

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_rotating_solid_body ( hdr_out, xi3_at_theta, xi3_at_rho, &
                             xi3_at_u, xi3_at_v, fields_out, field_count_out )

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

USE um_parvars, ONLY: &
    gc_all_proc_group

USE model_domain_mod, ONLY: &
    l_cartesian

USE rcf_nlist_recon_idealised_mod, ONLY: &
    l_rotate_grid,                       &
    grid_np_lon,                         &
    grid_np_lat,                         &
    aa_jet_u0,                           &
    aa_jet_a,                            &
    aa_jet_m,                            &
    aa_jet_n

USE idealise_run_mod, ONLY: &
    l_const_grav,           &
    l_shallow,              &
    t_surface,              &
    p_surface

USE rcf_ideal_tprofile_mod, ONLY: &
    tp_isothermal

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,     &
    stashcode_exner,        &
    stashcode_exner_surf,   &
    stashcode_dry_rho,      &
    stashcode_thetavd,      &
    stashcode_u,            &
    stashcode_v

USE planet_constants_mod, ONLY: &
    g,                          &
    kappa,                      &
    p_zero,                     &
    planet_radius,              &
    r,                          &
    two_omega

USE conversions_mod, ONLY: &
    pi_over_180

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
TYPE(field_type), POINTER :: fields_out(:)       ! output fields
INTEGER, INTENT(IN)       :: field_count_out     ! no. of output fields

! Local variables
INTEGER :: i
INTEGER :: k

REAL    :: xi1_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi2_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi1_u_global(output_grid % glob_u_row_length, &
                        output_grid % glob_u_rows)
REAL    :: xi1_p(output_grid % loc_p_row_length,         &
                 output_grid % loc_p_rows)
REAL    :: xi2_p(output_grid % loc_p_row_length,         &
                 output_grid % loc_p_rows)
REAL    :: xi1_u(output_grid % loc_u_row_length,         &
                 output_grid % loc_u_rows)
REAL    :: xi1_p_1d(output_grid % loc_p_field)
REAL    :: xi2_p_1d(output_grid % loc_p_field)
REAL    :: xi1_u_1d(output_grid % loc_u_field)

REAL    :: exner_ext(0:output_grid % model_levels+1) ! Extended 1D exner profile
REAL    :: dummy_dry_rho(output_grid % model_levels)   ! Dummy routines for 
REAL    :: dummy_thetavd(0:output_grid % model_levels) ! rcf_ideal_1d_profiles

! Variables for the rotated grid:
REAL :: lon_p
REAL :: lat_p
REAL :: sin_lon_p
REAL :: cos_lon_p
REAL :: sin_lat_p
REAL :: cos_lat_p
REAL :: sin_lon
REAL :: cos_lon
REAL :: sin_lat
REAL :: cos_lat
REAL :: s       ! = for deep or shallow atmosphere cos_lat derivation

REAL :: grav_1d(0:output_grid % model_levels)  ! Gravity profile
REAL :: t_surf
REAL :: p_surf
REAL :: u0_const
REAL :: f_sb
REAL :: dz

! Indices for locating fields
INTEGER  :: pos_dry_rho
INTEGER  :: pos_exner
INTEGER  :: pos_exner_surf
INTEGER  :: pos_thetavd
INTEGER  :: pos_u
INTEGER  :: pos_v

INTEGER  :: pe

! Trigonometric terms:
REAL     :: csxi1_p(output_grid % loc_p_field)
REAL     :: snxi1_p(output_grid % loc_p_field)
REAL     :: csxi2_p(output_grid % loc_p_field)
REAL     :: snxi2_p(output_grid % loc_p_field)
REAL     :: csxi1_u(output_grid % loc_u_field)
REAL     :: snxi1_u(output_grid % loc_u_field)


! Pointers to output fields:
TYPE( field_type ), POINTER        :: exner_out
TYPE( field_type ), POINTER        :: thetavd_out
TYPE( field_type ), POINTER        :: dry_rho_out
TYPE( field_type ), POINTER        :: exner_surf_out
TYPE( field_type ), POINTER        :: u_out
TYPE( field_type ), POINTER        :: v_out

! Values required for rcf_ideal_1d_profs
REAL, PARAMETER :: depth_factor = 1.0
REAL,    PARAMETER :: dtheta_dz1(3) = (/0.0, 0.0, 0.0/)
REAL,    PARAMETER :: height_dz1(2) = (/0.0, 0.0/)
LOGICAL, PARAMETER :: l_constant_dz = .TRUE.

CHARACTER (LEN=*),  PARAMETER :: RoutineName='RCF_IDEAL_ROTATING_SOLID_BODY'
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
                     xi1_u_global)

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


! Idealised code expects 1D arrays.
xi1_p_1d = RESHAPE(xi1_p, (/output_grid % loc_p_field/))
xi2_p_1d = RESHAPE(xi2_p, (/output_grid % loc_p_field/))
xi1_u_1d = RESHAPE(xi1_u, (/output_grid % loc_u_field/))

! Terms adapted from set_metric_terms_4A:
IF (l_cartesian) THEN

  DO i = 1, output_grid % loc_p_field
    csxi1_p(i) = 1.0/planet_radius
    snxi1_p(i) = 0.0
    csxi2_p(i) = 1.0/planet_radius
    snxi2_p(i) = 0.0
  END DO
  DO i = 1, output_grid % loc_u_field
    csxi1_u(i) = 1.0/planet_radius
    snxi1_u(i) = 0.0
  END DO

ELSE

  DO i = 1, output_grid % loc_p_field
    csxi1_p(i) = COS(xi1_p_1d(i))
    snxi1_p(i) = SIN(xi1_p_1d(i))
    csxi2_p(i) = COS(xi2_p_1d(i))
    snxi2_p(i) = SIN(xi2_p_1d(i))
  END DO
  DO i = 1, output_grid % loc_u_field
    csxi1_u(i) = COS(xi1_u_1d(i))
    snxi1_u(i) = SIN(xi1_u_1d(i))
  END DO

END IF


! Begin rotating solid body test case
!
! Set constant or variable gravity
IF (l_const_grav) THEN
  DO k = 0, output_grid % model_levels
    grav_1d(k) = g
  END DO
ELSE
  DO k = 0, output_grid % model_levels
    grav_1d(k) = g*( planet_radius/xi3_at_theta(1,k) )**2
  END DO
END IF

! Make copies of initial temp and pressure variables that can be modified.
t_surf = t_surface
p_surf = p_surface

CALL rcf_ideal_1d_profiles( exner_ext, dummy_dry_rho, dummy_thetavd,         &
                            tp_isothermal, t_surf, p_surf, 1.0,              &
                            xi3_at_theta(1,:), xi3_at_rho(1,:),              &
                            depth_factor, planet_radius,                     &
                            dtheta_dz1, height_dz1, l_constant_dz,           &
                            grav_1d, output_grid % model_levels)


u0_const = aa_jet_u0 * ( aa_jet_u0 + two_omega*planet_radius)/(r*t_surface)

IF (l_rotate_grid) THEN

  lon_p = grid_np_lon * pi_over_180
  lat_p = grid_np_lat * pi_over_180

  sin_lon_p = SIN(lon_p)
  cos_lon_p = COS(lon_p)
  sin_lat_p = SIN(lat_p)
  cos_lat_p = COS(lat_p)

! N.B. thetavd_out is allocated over 1:model_levels+1, not 0:model_levels,
!      so all second dimension indices for thetavd appear out by +1.
  DO i = 1, output_grid % loc_p_field
    sin_lon = snxi1_p(i) * cos_lon_p - csxi1_p(i) * sin_lon_p
    cos_lon = csxi1_p(i) * cos_lon_p + snxi1_p(i) * sin_lon_p
    sin_lat = snxi2_p(i) * sin_lat_p - cos_lon * csxi2_p(i) * cos_lat_p
    cos_lat = SQRT(1.0 - sin_lat * sin_lat)

    ! Thetavd at surface
    f_sb = 0.5 * u0_const * cos_lat**2
    ! Modify p_surface to take account of orography
    exner_surf_out % data(i,1) = p_surface * EXP(f_sb) *                     &
           EXP(-(xi3_at_theta(i,0) - planet_radius) * g/(r*t_surface) )
    thetavd_out % data(i,1) = t_surface /                                    &
                            (exner_surf_out % data(i,1) / p_zero)**kappa

    ! Modify hydrostatic pressure, thetavd and dry rho
    DO k = 1, output_grid % model_levels
      IF (l_shallow) THEN
        s = cos_lat
      ELSE
        s = xi3_at_rho(i,k) * cos_lat / planet_radius
      END IF
      f_sb = 0.5 * u0_const * s**2
      dz = xi3_at_rho(i,k) - planet_radius
      exner_out % data(i,k) = exner_ext(k) * EXP(kappa * f_sb)
      thetavd_out % data(i,k+1) = 2.0 * t_surface / exner_out % data(i,k)    &
                                 - thetavd_out % data(i,k)
      dry_rho_out % data(i,k) = p_zero                                       &
           * exner_out % data(i,k)**((1.0-kappa)/kappa)                      &
           / (r * (intw_w2rho(k,2) * thetavd_out % data(i,k)                 &
                 + intw_w2rho(k,1) * thetavd_out % data(i,k+1) ) )
    END DO
  END DO

  DO i = 1, output_grid % loc_u_field
    sin_lon = snxi1_u(i) * cos_lon_p - csxi1_u(i) * sin_lon_p
    cos_lon = csxi1_u(i) * cos_lon_p + snxi1_u(i) * sin_lon_p

    DO k = 1, output_grid % model_levels
      u_out % data(i,k) = aa_jet_u0 * xi3_at_u(i,k)                          &
           * (csxi2_p(i) * sin_lat_p + cos_lon * snxi2_p(i) * cos_lat_p)     &
           / planet_radius
    END DO

  END DO

  DO i = 1, output_grid % loc_v_field
    sin_lon = snxi1_p(i) * cos_lon_p - csxi1_p(i) * sin_lon_p
    cos_lon = csxi1_p(i) * cos_lon_p + snxi1_p(i) * sin_lon_p

    DO k = 1, output_grid % model_levels
      v_out % data(i,k) = aa_jet_u0 * xi3_at_v(i,k) * sin_lon * cos_lat_p    &
                          / planet_radius
    END DO

  END DO

ELSE  ! (l_rotate_grid = .FALSE.)

  ! Modify hydrostatic pressure, thetavd and dry rho
  DO i = 1, output_grid % loc_p_field
    f_sb = (aa_jet_u0 * csxi2_p(i)**aa_jet_m)**2  *                            &
           (1.0 / (2 * aa_jet_m) - 2.0 * (aa_jet_a * csxi2_p(i))**aa_jet_n     &
             / (2 * aa_jet_m + aa_jet_n) + 0.5 *                               &
             (aa_jet_a * csxi2_p(i))**(2 * aa_jet_n) / (aa_jet_m + aa_jet_n) ) &
         + two_omega * planet_radius * aa_jet_u0 * csxi2_p(i)**(aa_jet_m + 1)  &
          * (1.0 / (aa_jet_m + 1) - (aa_jet_a * csxi2_p(i))**aa_jet_n          &
          / (aa_jet_m + aa_jet_n + 1) )
    f_sb = f_sb / (r * t_surface)

    ! Modify p_surface to take account of orography
    exner_surf_out % data(i,1) = p_surface * EXP(f_sb) *                     &
                   EXP(-(xi3_at_theta(i,0) - planet_radius) * g/(r*t_surface))
    thetavd_out % data(i,1) = t_surface /                                    &
                            (exner_surf_out % data(i,1) / p_zero)**kappa
  END DO

  DO k = 1, output_grid % model_levels
    DO i = 1, output_grid % loc_p_field
      IF (l_shallow) THEN
        s = csxi2_p(i)
      ELSE
        s = xi3_at_rho(i,k) * csxi2_p(i) / planet_radius
      END IF
      f_sb = (aa_jet_u0 * s**aa_jet_m)**2 *                                  &
             (1.0 / (2 * aa_jet_m) - 2*(aa_jet_a * s)**aa_jet_n              &
               / (2 * aa_jet_m + aa_jet_n) + 0.5 *                           &
                 (aa_jet_a * s)**(2 * aa_jet_n) / (aa_jet_m + aa_jet_n) )    &
           + two_omega * planet_radius * aa_jet_u0 * s**(aa_jet_m + 1)       &
             * (1.0 / (aa_jet_m + 1) - (aa_jet_a * s)**aa_jet_n              &
               / (aa_jet_m + aa_jet_n + 1) )
      f_sb = f_sb / (r * t_surface)
      dz = xi3_at_rho(i,k) - planet_radius

      exner_out % data(i,k) = exner_ext(k) * EXP(kappa * f_sb)
      thetavd_out % data(i,k+1) = 2.0 * t_surface / exner_out % data(i,k)    &
                              - thetavd_out % data(i,k)
      dry_rho_out % data(i,k) = p_zero                                       &
             * exner_out % data(i,k)**((1.0 - kappa) / kappa)                &
             / (r * 0.5 * (thetavd_out % data(i,k) + thetavd_out % data(i,k+1)))
    END DO
  END DO

  ! Set winds.
  ! This is potentially unsafe as csxi2_p is defined on loc_p_field, whereas
  ! the usage is over loc_u_field. The same flaw was contained in the
  ! original atmosphere code this was derived from.
  IF (l_shallow) THEN
    DO k = 1, output_grid % model_levels
      DO i = 1, output_grid % loc_u_field
        u_out % data(i,k) = aa_jet_u0 * csxi2_p(i)**aa_jet_m                 &
                          * (1.0 - (aa_jet_a * csxi2_p(i))**aa_jet_n)
      END DO
    END DO
  ELSE
    DO k = 1, output_grid % model_levels
      DO i = 1, output_grid % loc_u_field
        ! As above but replacing csxi2_p with:
        s = xi3_at_u(i,k) * csxi2_p(i) / planet_radius
        u_out % data(i,k) = aa_jet_u0 * s**aa_jet_m                          &
                          * (1.0 - (aa_jet_a * s)**aa_jet_n)
      END DO
    END DO
  END IF

  DO k = 1, output_grid % model_levels
    DO i = 1, output_grid % loc_v_field
      v_out % data(i,k) = 0.0
    END DO
  END DO


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
END SUBROUTINE rcf_ideal_rotating_solid_body
END MODULE rcf_ideal_rotating_solid_body_mod
