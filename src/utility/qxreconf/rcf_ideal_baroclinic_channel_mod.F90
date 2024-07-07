! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Idealised baroclinic channel setup

MODULE rcf_ideal_baroclinic_channel_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_BAROCLINIC_CHANNEL_MOD'

CONTAINS

! Subroutine rcf_ideal_baroclinic_channel
!
! Description:
!   Sets wind field, pressure, potential temperature and density
!   for a baroclinic wave test in a Cartesian channel.
!   Isothermal case only. Version based on Ullrich & Jablonowski test.
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

SUBROUTINE rcf_ideal_baroclinic_channel ( hdr_out, xi3_at_theta, xi3_at_rho, &
                                xi3_at_u, g_theta, fields_out, field_count_out)

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

USE decomp_params, ONLY: &
    decomp_rcf_output

USE um_parvars, ONLY: &
    gc_all_proc_group

USE rcf_nlist_recon_idealised_mod, ONLY: &
    l_baro_perturbed

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,     &
    stashcode_exner,        &
    stashcode_exner_surf,   &
    stashcode_dry_rho,      &
    stashcode_thetavd,      &
    stashcode_u

USE rcf_ideal_baroclinic_constants_mod, ONLY: &
    baro_phi0,                                &
    beta0,                                    &
    f0,                                       &
    l_channel,                                &
    Lx,                                       &
    Ly,                                       &
    Lp,                                       &
    t0

USE rcf_ideal_baro_eta_conv_mod, ONLY: &
    baro_eta_conv

USE rcf_ideal_baro_u_p_mod, ONLY: &
    baro_u_p

USE lam_config_inputs_mod, ONLY: &
    delta_lat,                   &
    delta_lon

USE planet_constants_mod, ONLY: &
    kappa,                      &
    omega,                      &
    planet_radius

USE conversions_mod, ONLY: &
    pi_over_180

USE rcf_ideal_tprofile_mod, ONLY: &
    tp_BruntV

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY: &
    umPrint,          &
    umMessage,        &
    newline,          &
    PrintStatus,      &
    PrStatus_Diag

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
REAL, INTENT(IN)          :: g_theta(:,:)        ! gravity at theta_levels
TYPE(field_type), POINTER :: fields_out(:)       ! output fields
INTEGER, INTENT(IN)       :: field_count_out     ! no. of output fields

! Local variables
INTEGER :: i
INTEGER :: k

REAL    :: x0
REAL    :: y0 
REAL    :: t_guess   ! First guess at temp. sol'n.
REAL    :: p_guess   ! First guess at pressure sol'n.
REAL    :: eta_surface
REAL    :: eta_p_u(output_grid % model_levels)
REAL    :: xi2_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi1_u_global(output_grid % glob_u_row_length, &
                        output_grid % glob_u_rows)
REAL    :: xi2_v_global(output_grid % glob_v_row_length, &
                        output_grid % glob_v_rows)
REAL    :: xi2_p(output_grid % loc_p_row_length,         &
                 output_grid % loc_p_rows)
REAL    :: xi1_u(output_grid % loc_u_row_length,         &
                 output_grid % loc_u_rows)
REAL    :: xi2_p_1d(output_grid % loc_p_field)
REAL    :: xi1_u_1d(output_grid % loc_u_field)
REAL    :: exner_ext(0:output_grid % model_levels+1) ! Extended 1D exner profile

CHARACTER (LEN=*), PARAMETER :: RoutineName = 'RCF_IDEAL_BAROCLINIC_CHANNEL'

! Indices for locating fields
INTEGER                            :: pos_thetavd
INTEGER                            :: pos_dry_rho
INTEGER                            :: pos_exner
INTEGER                            :: pos_exner_surf
INTEGER                            :: pos_u
INTEGER                            :: pe

! Pointers to output fields:
TYPE( field_type ), POINTER        :: exner_out
TYPE( field_type ), POINTER        :: thetavd_out
TYPE( field_type ), POINTER        :: dry_rho_out
TYPE( field_type ), POINTER        :: exner_surf_out
TYPE( field_type ), POINTER        :: u_out

! Values required for rcf_ideal_1d_profs
REAL,    PARAMETER :: dtheta_dz1(3) = (/0.0, 0.0, 0.0/)
REAL,    PARAMETER :: height_dz1(2) = (/0.0, 0.0/)
LOGICAL, PARAMETER :: l_constant_dz = .TRUE.

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

! Calculate global grid info and scatter to generate local data
CALL rcf_calc_coords(hdr_out, output_grid, xi2_p_global=xi2_p_global, &
                     xi1_u_global=xi1_u_global, xi2_v_global=xi2_v_global)

pe = 0

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

! No need to process xi2_v_global here.

! Idealised code expects 1D arrays.
xi2_p_1d = RESHAPE(xi2_p, (/output_grid % loc_p_field/))
xi1_u_1d = RESHAPE(xi1_u, (/output_grid % loc_u_field/))

! Begin baroclinic channel test case
!
! Set up constants held in the baroclinic module.
! delta_lat and delta_lon needed because e.g. xi1_u_global(glob_u_row_length, 1)
! here is equivalent to xi1_u_global(glob_u_row_length-1, 1) in the model.
l_channel = .TRUE.
Lx = xi1_u_global(output_grid % glob_u_row_length, 1) + delta_lon    &
     - xi1_u_global(1, 1)
Ly = xi2_v_global(1, output_grid % glob_v_rows) + delta_lat - xi2_v_global(1, 1)
Lp = 600.0 * 1000.0
f0 = 2.0 * omega * SIN(baro_phi0 * pi_over_180)
beta0 = 0.0   ! beta plane not implemented in ENDGame, otherwise would be:
              ! 2.0 * omega/planet_radius * COS(baro_phi0 * pi_over_180)
x0 = 2000.0 * 1000.0
y0 = 2500.0 * 1000.0

IF (PrintStatus >= PrStatus_Diag .AND. mype == 0) THEN
  CALL umPrint('Baroclinic Channel wave test parameters:', src=RoutineName)
  WRITE(umMessage,'(7(A,E12.5,A))')  &
    'Lx    = ', Lx,        newline,  &
    'Ly    = ', Ly,        newline,  &
    'phi0  = ', baro_phi0, newline,  &
    'f0    = ', f0,        newline,  &
    'beta0 = ', beta0,     newline,  &
    'x0    = ', x0,        newline,  &
    'y0    = ', y0
  CALL umPrint(umMessage, src=RoutineName)
END IF

! Initial guesses for temp and pressure
t_guess = t0
p_guess = 1.0e5

DO i = 1, thetavd_out % level_size
  CALL rcf_ideal_1d_profiles( exner_ext, dry_rho_out % data(i,:),            &
                              thetavd_out % data(i,:), tp_BruntV,            &
                              t_guess, p_guess, xi2_p_1d(i),                 &
                              xi3_at_theta(i,:), xi3_at_rho(i,:),            &
                              output_grid % z_top_of_model, planet_radius,   &
                              dtheta_dz1, height_dz1, l_constant_dz,         &
                              g_theta(i,:), output_grid % model_levels)

!  Copy relevant part of extended exner array back to the actual field
  DO k = 1, output_grid % model_levels+1
    exner_out % data(i, k) = exner_ext(k)
  END DO

! eta_surface should be 1, but calculate it for safety:
  eta_surface = baro_eta_conv(xi2_p_1d(i), xi3_at_theta(i,0) - planet_radius)
  exner_surf_out % data(i,1) = eta_surface**kappa
END DO

! Calculate winds:
DO k = 1, output_grid % model_levels
  DO i = 1, u_out % level_size
    eta_p_u(k) = baro_eta_conv(xi2_p_1d(i), xi3_at_u(i,k) - planet_radius)
    u_out % data(i,k) = baro_u_p(xi2_p_1d(i), eta_p_u(k))
! Add perturbation
   IF (l_baro_perturbed) THEN
    u_out % data(i,k) = u_out % data(i,k) + &
          EXP( -(  (xi1_u_1d(i) - x0)**2 + (xi2_p_1d(i) - y0)**2 ) / Lp**2 )
   END IF
  END DO
END DO

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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_baroclinic_channel
END MODULE rcf_ideal_baroclinic_channel_mod
