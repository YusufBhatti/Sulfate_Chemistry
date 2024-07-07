! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Read the recon_idealised namelist

MODULE rcf_readnl_recon_idealised_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_READNL_RECON_IDEALISED_MOD'

CONTAINS

! Subroutine rcf_readnl_recon_idealised
!
! Description:
!   Read, print and check the recon_idealised namelist.
!   Update various other settings as appropriate.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.

SUBROUTINE rcf_readnl_recon_idealised (unit_no)

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_generate_heights_mod, ONLY: &
    height_gen_linear,              &
    height_gen_smooth

USE umPrintMgr,   ONLY: &
    PrintStatus,        &
    PrStatus_Oper,      &
    umMessage,          &
    umPrint

USE UM_ParCore,   ONLY: &
    mype

USE rcf_nlist_recon_idealised_mod, ONLY: &
    read_nml_recon_idealised,            &
    print_nlist_recon_idealised,         &
    check_nlist_recon_idealised,         &
    grid_flat,                           &
    grid_number,                         &
    height_domain,                       &
    big_layers,                          &
    transit_layers,                      &
    base_xi1,                            &
    base_xi2,                            &
    delta_xi1,                           &
    delta_xi2,                           &
    l_init_idealised,                    &
    first_constant_r_rho_level_new

USE rcf_ideal_vgrid_mod, ONLY:  &
    vert_idl_um_grid,           &
    vert_quadratic_uvtheta,     &
    vert_stretch_plus_regular

USE rcf_ideal_gflat_mod, ONLY: &
    gflat_linear,              &
    gflat_quadratic

USE rcf_ideal_set_eta_levels_mod, ONLY: &
    rcf_ideal_set_eta_levels

USE dynamics_testing_mod, ONLY: &
    l_idealised_data,           &
    problem_number

USE problem_mod, ONLY: &
    idealised_planet

USE lam_config_inputs_mod, ONLY: &
    delta_lat,                   &
    delta_lon,                   &
    frstlata,                    &
    frstlona,                    &
    target_delta_lat,            &
    target_delta_lon,            &
    target_frstlata,             &
    target_frstlona

USE model_domain_mod, ONLY: &
    l_regular,              &
    l_cartesian,            &
    model_type,             &
    mt_global,              &
    mt_cyclic_lam

USE ereport_mod,  ONLY: &
    ereport

USE yomhook,   ONLY: &
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: unit_no   ! Unit no. for file containing recon_idealised

! Local variables
INTEGER :: i
INTEGER :: fcrrl ! Original value of first_constant_r_rho_level
INTEGER :: icode ! Error reporting
LOGICAL :: update_grid   ! Reset grid size and position
CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_READNL_RECON_IDEALISED'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL read_nml_recon_idealised(unit_no)

IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
  CALL print_nlist_recon_idealised()
END IF

CALL check_nlist_recon_idealised()


! Option 11 vert_idl_um_grid makes the idealised vertical level setting 
! behave like all other UM configurations i.e. not using the idealised
! namelist values of grid_flat, height_domain and first_constant_r_rho_level_new
! All other options may overide the standard behaviour in some way.

IF ( grid_number == vert_idl_um_grid) THEN

  ! Just set height generation method to value used standard UM models
  Output_Grid % height_gen_method = height_gen_smooth

ELSE

  ! Reset some values from the recon_idealised namelist
  ! --------------------------------------------
  ! first_constant_r_rho_level:
  fcrrl = Output_Grid % first_constant_r_rho_level

  IF ( grid_number  ==  vert_stretch_plus_regular .AND. big_layers  >   0) THEN
    Output_Grid % first_constant_r_rho_level =                           &
      MIN(Output_Grid % first_constant_r_rho_level,                      &
          Output_Grid % model_levels - big_layers - transit_layers)
  ELSE
    IF (first_constant_r_rho_level_new > 0) THEN
      Output_Grid % first_constant_r_rho_level = MIN(  &
      Output_Grid % first_constant_r_rho_level, first_constant_r_rho_level_new)
    END IF
  END IF

  IF (Output_Grid % first_constant_r_rho_level /= fcrrl) THEN
    IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
      WRITE(umMessage,'(2(A,I0))')                                           &
                                'Resetting first_constant_r_rho_level from ',&
                      fcrrl, ' to ', Output_Grid % first_constant_r_rho_level
      CALL umPrint(umMessage, src=RoutineName)
    END IF
  END IF

  ! Overwrite output dump values of eta_theta_levels & eta_rho_levels:
  CALL rcf_ideal_set_eta_levels (Output_Grid % model_levels,             &
                                 Output_Grid % eta_theta_levels,         &
                                 Output_Grid % eta_rho_levels)

  ! Reset z_top_of_model
  Output_Grid % z_top_of_model = height_domain
  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE(umMessage,'((A,E12.5))') 'z_top_of_model reset to height_domain = ', &
                                   height_domain
    CALL umPrint(umMessage, src=RoutineName)
  END IF

  ! Set the height generation method.
  ! The idealised quadratic structure is equivalent to 
  ! the reconfiguration's "smooth" height generation.

  IF (grid_number <= vert_quadratic_uvtheta) THEN
    SELECT CASE(grid_flat)

    CASE(gflat_linear)
      Output_Grid % height_gen_method = height_gen_linear

    CASE(gflat_quadratic)
      Output_Grid % height_gen_method = height_gen_smooth

    CASE DEFAULT
      icode = 100
      WRITE(umMessage,'(A,I0)') 'Unsupported value of grid_flat: ',grid_flat
      CALL ereport(RoutineName,icode,umMessage)

    END SELECT

  ELSE
    icode = 200
    WRITE(umMessage,'(A,I0)') 'Unsupported value of grid_number: ',grid_number
    CALL ereport(RoutineName,icode,umMessage)
  END IF
END IF

! Under certain circumstances the basic grid structure/position is taken from
! the recon_idealised namelist. 
update_grid = .FALSE.
IF (l_idealised_data)  THEN

  IF (l_regular) THEN

    IF ( (model_type /= mt_global .AND. model_type /= mt_cyclic_lam)     &
       .OR. l_cartesian) THEN

      update_grid = .TRUE.

      ! Save the original target LAM co-ordinates for possible use with
      ! interpolation
      IF (l_cartesian) THEN
        target_frstlona  = frstlona
        target_delta_lon = delta_lon
        target_frstlata  = frstlata
        target_delta_lat = delta_lat
      END IF

      ! Set the new LAM co-ordinates from the recon_idealised namelist
      frstlona  = base_xi1
      delta_lon = delta_xi1
      IF (model_type /= mt_global) THEN
        frstlata  = base_xi2
        delta_lat = delta_xi2
      END IF

    END IF ! l_cartesian or not (global or cyclic)

  END IF ! l_regular

ELSE    ! Not l_idealised_data but starting from a Cartesian grid
        ! and wanting to go to a new Cartesian grid

  IF (l_cartesian) THEN
    ! Set the new LAM co-ordinates from the recon_idealised namelist
    frstlona  = base_xi1
    delta_lon = delta_xi1
    frstlata  = base_xi2
    delta_lat = delta_xi2
  END IF
  
END IF ! idealised data/planet

IF (update_grid) THEN
  IF (PrintStatus >= PrStatus_Oper .AND. mype == 0) THEN
    WRITE(umMessage,'(A)') &
        'Resetting grid size and position from recon_idealised namelist:'
    CALL umPrint(umMessage, src=RoutineName)
    WRITE(umMessage, '(2(A,E12.5))') &
      ' first longitude: ',frstlona, ', first latitude: ',frstlata
    CALL umPrint(umMessage, src=RoutineName)
    WRITE(umMessage, '(2(A,E12.5))') &
      ' column spacing: ',delta_lon, ', row spacing: ',delta_lat
    CALL umPrint(umMessage, src=RoutineName)
  END IF

  ! Idealised model only supports ENDGame grids.
  DO i = 1, output_grid % glob_p_row_length
    output_grid % lambda_p(i) = frstlona + (i-0.5)*delta_lon
  END DO
  DO i = 1, output_grid % glob_u_row_length
    output_grid % lambda_u(i) = frstlona + (i-1)*delta_lon
  END DO
  DO i = 1, output_grid % glob_p_rows
    output_grid % phi_p(i)    = frstlata + (i-0.5)*delta_lat
  END DO
  DO i = 1, output_grid % glob_v_rows
    output_grid % phi_v(i)    = frstlata + (i-1)*delta_lat
  END DO
END IF

! If start dump contents are to be overwritten with an idealised setup,
! avoid rotating the output grid.
IF (l_init_idealised) THEN
  output_grid % rotated = .FALSE.
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE  rcf_readnl_recon_idealised
END MODULE rcf_readnl_recon_idealised_mod
