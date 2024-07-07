! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Idealised isothermal setup

MODULE rcf_ideal_isothermal_profile_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_IDEAL_ISOTHERMAL_PROFILE_MOD'

CONTAINS

! Subroutine rcf_ideal_isothermal_profile
!
! Description:
!   A simple wrapper routine to set up a series of 1D profiles, assuming
!   isothermal pressure.
!   Winds are not initialised and are assumed to have been set elsewhere.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_ideal_isothermal_profile ( hdr_out, xi3_at_theta, xi3_at_rho, &
                                        g_theta, fields_out, field_count_out,&
                                        print_first)

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

USE um_stashcode_mod, ONLY: &
    stashcode_prog_sec,     &
    stashcode_exner,        &
    stashcode_exner_surf,   &
    stashcode_dry_rho,      &
    stashcode_thetavd

USE rcf_nlist_recon_idealised_mod, ONLY: &
    dtheta_dz1,                          &
    height_dz1,                          &
    l_constant_dz,                       &
    tprofile_number

USE idealise_run_mod, ONLY: &
    p_surface,              &
    t_surface

USE planet_constants_mod, ONLY: &
    g,                          &
    kappa,                      &
    p_zero,                     &
    planet_radius,              &
    r

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
REAL, INTENT(IN)          :: g_theta(:,:)        ! gravity at theta_levels
TYPE(field_type), POINTER :: fields_out(:)       ! output fields
INTEGER, INTENT(IN)       :: field_count_out     ! no. of output fields
LOGICAL, INTENT(IN), OPTIONAL :: print_first     ! If TRUE, print first profile

! Local variables
INTEGER :: i
INTEGER :: k

REAL    :: p_on_orog
REAL    :: xi2_p_global(output_grid % glob_p_row_length, &
                        output_grid % glob_p_rows)
REAL    :: xi2_p(output_grid % loc_p_row_length,         &
                 output_grid % loc_p_rows)
REAL    :: xi2_p_1d(output_grid % loc_p_field)
REAL    :: exner_ext(0:output_grid % model_levels+1) ! Extended 1D exner profile

LOGICAL :: print_profile   ! Flag for passing to rcf_ideal_1d_profiles.
                           ! Set from print_first.

! Indices for locating fields
INTEGER :: pos_thetavd
INTEGER :: pos_dry_rho
INTEGER :: pos_exner
INTEGER :: pos_exner_surf

INTEGER :: pe

! Pointers to output fields:
TYPE( field_type ), POINTER        :: exner_out
TYPE( field_type ), POINTER        :: exner_surf_out
TYPE( field_type ), POINTER        :: thetavd_out
TYPE( field_type ), POINTER        :: dry_rho_out

CHARACTER (LEN=*),  PARAMETER :: RoutineName='RCF_IDEAL_ISOTHERMAL_PROFILE'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Process optional argument
IF (PRESENT(print_first)) THEN
  print_profile = print_first
ELSE
  print_profile = .FALSE.
END IF

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

! Calculate global grid info and scatter to generate local data
CALL rcf_calc_coords(hdr_out, output_grid, xi2_p_global=xi2_p_global)

pe = 0

CALL rcf_scatter_field_real(xi2_p, xi2_p_global,                 &
                            output_grid % loc_p_row_length,      &
                            output_grid % loc_p_rows,            &
                            output_grid % glob_p_row_length,     &
                            output_grid % glob_p_rows, pe,       &
                            gc_all_proc_group )


! Idealised code expects 1D arrays.
xi2_p_1d = RESHAPE(xi2_p, (/output_grid % loc_p_field/))


! Set up each 1D profile:
DO i = 1, thetavd_out % level_size
  ! Adjust pressure over orography using isothermal assumption:
  p_on_orog = p_surface * EXP(-(xi3_at_theta(i,0) - planet_radius) * g       &
                               / (r * t_surface))

  exner_surf_out % data(i,1) = (p_on_orog / p_zero)**kappa

  ! Find thermodynamics profiles:
  CALL rcf_ideal_1d_profiles( exner_ext, dry_rho_out % data(i,:),            &
                              thetavd_out % data(i,:), tprofile_number,      &
                              t_surface, p_on_orog, xi2_p_1d(i),             &
                              xi3_at_theta(i,:), xi3_at_rho(i,:),            &
                              1.0, planet_radius,                            &
                              dtheta_dz1, height_dz1, l_constant_dz,         &
                              g_theta(i,:), output_grid % model_levels,      &
                              print_profile)

  !  Copy relevant part of extended exner array back to the actual field.
  DO k = 1, output_grid % model_levels+1
    exner_out % data(i, k) = exner_ext(k)
  END DO

  ! Only print the first profile
  print_profile = .FALSE.
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

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE rcf_ideal_isothermal_profile
END MODULE rcf_ideal_isothermal_profile_mod
