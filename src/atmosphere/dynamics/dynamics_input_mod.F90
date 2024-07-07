! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module containing input switches/settings as used by the
!   semi_lagrangian code and the check_run_dyn routine for logic
!   checking the selected settings.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics

! Method:
!   Switches are initialised to false and read in from the
!   namelists. The module may then be used directly where the switches
!   are needed within the dynamics code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE dynamics_input_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE gcr_input_mod, ONLY:                                         &
    L_GCR_cycle_opt,GCR_zero_init_guess,GCR_use_residual_tol,    &
    GCR_adi_add_full_soln,GCR_use_tol_abs,L_gcr_fast_x,          &
    GCR_max_iterations,GCR_diagnostics,GCR_precon_option,        &
    GCR_n_ADI_pseudo_timesteps,GCR_Restart_value,GCR_tol_res,    &
    GCR_tol,GCR_tol_res2,GCR_tol_abs,GCR_tol_abs2,GCR_tol2,      &
    GCR_ADI_pseudo_timestep,G_term_tol,GCR_its_avg_step

USE lbc_input_mod, ONLY:                                         &
    L_LBC_balance,L_lbc_new,L_transparent,L_lbc_old

USE eg_alpha_ramp_mod, ONLY:                                     &
     alpha_relax_type

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! ------------------------------------------------
! Parameters not read in by the run_dyn namelist.
! ------------------------------------------------


! T: reset polar values to mean value every polar_reset_timesteps
LOGICAL, PARAMETER :: L_polar_reset   = .FALSE.
LOGICAL, PARAMETER :: L_interp_depart = .FALSE.
                                      ! interpolate to find u,v departure pts
LOGICAL, PARAMETER :: L_qwaterload    = .TRUE.
                                      ! true if using adding waterloading terms
LOGICAL, PARAMETER :: L_fint_theta    = .FALSE. ! true:  fully-interpolating
                                      ! semi-lagrangian theta advection will be
                                      ! used false: standard non-interpolating
                                      ! in the vertical

! interval for resetting polar to mean
INTEGER, PARAMETER :: polar_reset_timesteps = 1

REAL, PARAMETER :: tol_sc_fact  = 1.0      ! linear solver scaling factor
REAL, PARAMETER :: T_surf_ref   = 290.0    ! SI Reference surface temperature
REAL, PARAMETER :: p_surf_ref   = 1.0e5    ! SI Reference surface pressure
REAL, PARAMETER :: extrp_weight = 1.5      ! time interpolation weight for first
                                           ! SL iteration.
REAL, PARAMETER :: eccentricity = 0.0      ! Eccentricity of planet, considered
                                           ! as an oblate spheroid

INTEGER :: NumCycles    = 2    !  Parameters for dynamics-physics cycling
                               ! only valid value for ENDGame.

! ===============================================
REAL :: alpha_1    = rmdi                                 ! (To be retired)
REAL :: alpha_2    = rmdi                                 ! (To be retired)
REAL :: alpha_3    = rmdi                                 ! (To be retired)
REAL :: alpha_4    = rmdi                                 ! (To be retired)
REAL :: alpha_1_2  = rmdi                                 ! (To be retired)
REAL :: alpha_2_2  = rmdi                                 ! (To be retired)
REAL :: alpha_3_2  = rmdi                                 ! (To be retired)
REAL :: alpha_4_2  = rmdi                                 ! (To be retired)
LOGICAL :: L_new_tdisc    = .FALSE.                      ! (To be retired)
LOGICAL :: L_thmono_fixed = .FALSE.                      ! (To be retired)
INTEGER :: i_ND_solver_vn = imdi                          ! (To be retired)
! ===============================================


! -------------------------------
! Items read in by run_dyn namelist.
! -------------------------------

LOGICAL :: L_check_moist_inc = .FALSE.  ! T: Check various physics schemes
                                        ! conserve moisture and print global
                                        ! errors. (Only use on short runs.)

! ENDGAME only namelist items
LOGICAL :: l_simple_coriolis    =.FALSE. ! option to use a simplified
                                         ! representation of the Coriolis terms
LOGICAL :: L_eg_dry_static_adj  =.FALSE. ! Enforces static stability on
                                          ! input theta field

LOGICAL :: L_inc_solver         =.FALSE. ! Flag for incremental solver
LOGICAL :: L_conserv_smooth_lap =.FALSE. ! to use the option of redistributing
                                          ! mass based on only the Laplacian in
                                          ! the mass conservation routine

LOGICAL :: L_accel_convergence  =.FALSE. ! accelerated convergence switches
LOGICAL :: L_init_Fnm1          =.FALSE. ! Fnm1 initialisation
LOGICAL :: l_fast_vert_adjust   =.FALSE. ! Switch to use simplified
                                         ! adjust_vert_bound code
LOGICAL :: l_sl_bc_correction   =.FALSE. ! Flag for using corrected treatment
                                         ! of upper and lower boundaries in
                                         ! semi-Lagrangian advection scheme.
! Logicals and options for ZLF scheme
! also ensure there is error checking if needed on these variables
! in chk_run_dyn 
!
! The following logicals/integer parameters are for controlling the ZLF
! scheme to restore mass conservation at the advection stage.
! 
! L_conservation_(*)_zlf=T/F = use/not ZLF scheme for (*) variables
!
! zlf_conservation_(*)_option number specifies the mass conservation
! scheme to be used with the ZLF scheme
! zlf_conservation_(*)_option = 0 (ZLF consevation off )
! zlf_conservation_(*)_option = 1 (OCF scheme )
! zlf_conservation_(*)_option = 2 (ADAS scheme)
! .
INTEGER, PARAMETER :: apply_ocf  = 1
INTEGER, PARAMETER :: apply_adas = 2
INTEGER :: zlf_conservation_moist_option = imdi
INTEGER :: zlf_conservation_tracers_option = imdi
INTEGER :: zlf_conservation_theta_option = imdi

REAL    :: zlf_maximum_height            = rmdi

INTEGER :: eg_vert_damp_profile = imdi    ! Vertical damping flag
REAL :: eta_s                   = rmdi    ! Height (in eta) above
                                          ! which to apply damping

REAL :: Ih                      = rmdi    ! Hydrostatic switch
                                          !(1=nonhydrostatic; 0=hydrostatic)

LOGICAL :: L_RK_dps             =.TRUE.  ! Runge-Kutta departure point scheme
LOGICAL :: L_eliminate_rho      =.TRUE.  ! Diagnostic rho
INTEGER :: InnIts               = 2       ! Number of inner iterations
REAL :: eg_vert_damp_coeff      = rmdi
REAL :: damp_height             = rmdi

LOGICAL :: l_viscosity          = .FALSE. ! Include explicit Viscosity
REAL :: horiz_viscosity         = rmdi    ! Horizontal Viscosity Coeff
REAL :: vert_viscosity          = rmdi    ! Vertical Viscosity Coeff

! Conservation schemes for total mass of dry air
INTEGER :: conserve_dry_mass       = imdi
INTEGER, PARAMETER ::                                                          &
         not_conserved           = 0 ! No conservation of dry air
INTEGER, PARAMETER ::                                                          &
         constant_factor         = 1 ! Adjust dryrho by a constant factor
INTEGER, PARAMETER ::                                                          &
         linear_factor_IE        = 2 ! Adjust dryrho and theta_vd by a linear
                                     ! function of height, preserving total
                                     ! gravitational potential energy and
                                     ! total internal energy (approximately).
INTEGER, PARAMETER ::                                                          &
         linear_factor           = 3 ! Apply linear function to dryrho only.

! ----------------------
! namelist run_dyn
! ----------------------

NAMELIST/RUN_Dyn/                                                &
 L_check_moist_inc,                                              &
 GCR_max_iterations,GCR_precon_option,                           &
 GCR_tol,                                                        &
 conserve_dry_mass, l_sl_bc_correction,                          &
 alpha_relax_type,                                               &
 eg_vert_damp_profile, eta_s, damp_height, eg_vert_damp_coeff,   &
 l_viscosity, horiz_viscosity, vert_viscosity,                   &
 zlf_conservation_moist_option,                                  &
 zlf_conservation_tracers_option,                                &
 zlf_conservation_theta_option,                                  &
 zlf_maximum_height

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DYNAMICS_INPUT_MOD'

CONTAINS

SUBROUTINE check_run_dyn()

! Description:
!   Subroutine to apply logic checks and set variables based on the
!   options selected in the run_dyn namelist.
USE ereport_mod, ONLY: ereport

USE precon_constants_mod, ONLY:                                            &
    vert_plus_xyz_ADI_precon,xyz_ADI_precon,                               &
    vert_plus_xz_ADI_precon,xz_ADI_precon,no_precon

USE bl_option_mod,     ONLY: l_quick_ap2, i_bl_vn, i_bl_vn_1a
USE umPrintMgr, ONLY: umMessage, umPrint

IMPLICIT NONE

REAL(KIND=jprb)               :: zhook_handle

INTEGER                       :: ErrorStatus   ! used for ereport
CHARACTER (LEN=errormessagelength)     :: cmessage      ! used for ereport
CHARACTER (LEN=*), PARAMETER  :: RoutineName = 'CHECK_RUN_DYN'
INTEGER  :: zlf_opt1, zlf_opt2, zlf_opt3

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ---------------------------------------------
! Dynamics logic control - New Dynamics and ENDGAME
! ---------------------------------------------

! Skip all calculations in atmos_physics2 which are identical between
! the first and subsequent iterations
! Not compatable with 1a BL scheme as advection of TKE will differ
! between each outer loop
IF (i_bl_vn /= i_bl_vn_1a) l_quick_ap2 = .TRUE.

IF (GCR_use_residual_Tol) THEN
  GCR_use_tol_abs=.FALSE.
  GCR_tol_res = GCR_tol
  GCR_tol_res2 = GCR_tol2
ELSE
  GCR_use_tol_abs=.TRUE.
  GCR_tol_abs = GCR_tol
  GCR_tol_abs2 = GCR_tol2
END IF


IF (GCR_precon_option /= no_precon   .AND.                  &
    GCR_precon_option /= vert_plus_xyz_ADI_precon   .AND.   &
    GCR_precon_option /= vert_plus_xz_ADI_precon) THEN
  cmessage =                                                 &
  'An invalid preconditioner option for ENDGame has been requested' 
  ErrorStatus = 25
  CALL ereport(RoutineName, ErrorStatus, cmessage)
END IF

zlf_opt1= zlf_conservation_moist_option*(2-zlf_conservation_moist_option)
zlf_opt2= zlf_conservation_tracers_option*(2-zlf_conservation_tracers_option)
zlf_opt3= zlf_conservation_theta_option*(2-zlf_conservation_theta_option)
IF ( zlf_opt1 > 0 ) THEN
 WRITE(umMessage,'(A)') &
 ' Invalid zlf_conservation_moist_option - value reset to 0'
 CALL umPrint(umMessage,src='CHECK_RUN_DYN')
 zlf_conservation_moist_option = 0
END IF
IF ( zlf_opt2 > 0 ) THEN
 WRITE(umMessage,'(A)')  &
 'Invalid zlf_conservation_tracers_option - value reset to 0'
 CALL umPrint(umMessage,src='CHECK_RUN_DYN')
 zlf_conservation_tracers_option = 0
END IF
IF ( zlf_opt3 > 0 ) THEN
 WRITE(umMessage,'(A)') &
 'Invalid zlf_conservation_theta_option - value reset to 0'
 CALL umPrint(umMessage,src='CHECK_RUN_DYN')
 zlf_conservation_theta_option = 0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE check_run_dyn

SUBROUTINE print_nlist_run_dyn()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_DYN'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_dyn', &
    src='dynamics_input_mod')

WRITE(lineBuffer,*) ' L_check_moist_inc = ',L_check_moist_inc
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,*) ' GCR_max_iterations = ',GCR_max_iterations
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,*) ' GCR_precon_option = ',GCR_precon_option
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,*) ' GCR_tol = ',GCR_tol
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,'(A,I0)') ' conserve_dry_mass = ',conserve_dry_mass
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,*) ' l_sl_bc_correction = ',l_sl_bc_correction
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,*) ' alpha_relax_type = ',alpha_relax_type
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,*) ' eg_vert_damp_profile = ',eg_vert_damp_profile
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,'(A,ES17.10)') ' eg_vert_damp_coeff = ',eg_vert_damp_coeff
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,"(A,F16.3)") ' damp_height = ',damp_height
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,'(A,L1)') ' l_viscosity = ',l_viscosity
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,"(A,F16.3)") ' horiz_viscosity = ',horiz_viscosity
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,"(A,F16.3)") ' vert_viscosity = ',vert_viscosity
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,'(A,I0)') 'zlf_conservation_moist_option = ', &
                            zlf_conservation_moist_option
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,'(A,I0)') 'zlf_conservation_tracers_option = ', &
                            zlf_conservation_tracers_option
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,'(A,I0)') 'zlf_conservation_theta_option = ', &
                            zlf_conservation_theta_option
CALL umPrint(lineBuffer,src='dynamics_input_mod')
WRITE(lineBuffer,'(A,E16.3)') 'zlf_maximum_height = ',zlf_maximum_height
CALL umPrint(lineBuffer,src='dynamics_input_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='dynamics_input_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_dyn

#if !defined(LFRIC)
SUBROUTINE read_nml_run_dyn(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) :: unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_DYN'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int  = 8
INTEGER, PARAMETER :: n_real = 7
INTEGER, PARAMETER :: n_log  = 3

TYPE my_namelist
  SEQUENCE
  INTEGER :: GCR_max_iterations
  INTEGER :: GCR_precon_option
  INTEGER :: alpha_relax_type
  INTEGER :: eg_vert_damp_profile
  INTEGER :: zlf_conservation_moist_option
  INTEGER :: zlf_conservation_tracers_option
  INTEGER :: zlf_conservation_theta_option
  INTEGER :: conserve_dry_mass
  REAL :: eg_vert_damp_coeff
  REAL :: GCR_tol
  REAL :: horiz_viscosity
  REAL :: eta_s
  REAL :: damp_height
  REAL :: vert_viscosity
  REAL :: zlf_maximum_height
  LOGICAL :: L_check_moist_inc
  LOGICAL :: l_sl_bc_correction
  LOGICAL :: l_viscosity
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Dyn, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Dyn", iomessage)

  my_nml % GCR_max_iterations = GCR_max_iterations
  my_nml % GCR_precon_option  = GCR_precon_option
  my_nml % alpha_relax_type   = alpha_relax_type
  my_nml % eg_vert_damp_profile = eg_vert_damp_profile
  my_nml % zlf_conservation_moist_option = zlf_conservation_moist_option
  my_nml % zlf_conservation_tracers_option = zlf_conservation_tracers_option
  my_nml % zlf_conservation_theta_option = zlf_conservation_theta_option
  my_nml % conserve_dry_mass  = conserve_dry_mass
  ! end of integers
  my_nml % eg_vert_damp_coeff = eg_vert_damp_coeff
  my_nml % GCR_tol   = GCR_tol
  my_nml % horiz_viscosity = horiz_viscosity
  my_nml % eta_s        = eta_s
  my_nml % damp_height  = damp_height
  my_nml % vert_viscosity = vert_viscosity
  my_nml % zlf_maximum_height = zlf_maximum_height
  ! end of reals
  my_nml % L_check_moist_inc     = L_check_moist_inc
  my_nml % l_sl_bc_correction    = l_sl_bc_correction
  my_nml % l_viscosity    = l_viscosity

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  GCR_max_iterations   = my_nml % GCR_max_iterations
  GCR_precon_option    = my_nml % GCR_precon_option
  alpha_relax_type     = my_nml % alpha_relax_type
  eg_vert_damp_profile = my_nml % eg_vert_damp_profile
  zlf_conservation_moist_option = my_nml % zlf_conservation_moist_option
  zlf_conservation_tracers_option = my_nml % zlf_conservation_tracers_option
  zlf_conservation_theta_option = my_nml % zlf_conservation_theta_option
  conserve_dry_mass    = my_nml % conserve_dry_mass
  ! end of integers
  eg_vert_damp_coeff = my_nml % eg_vert_damp_coeff
  GCR_tol   = my_nml % GCR_tol
  horiz_viscosity = my_nml % horiz_viscosity
  eta_s        = my_nml % eta_s
  damp_height  = my_nml % damp_height
  vert_viscosity = my_nml % vert_viscosity
  zlf_maximum_height = my_nml % zlf_maximum_height
  ! end of reals
  L_check_moist_inc     = my_nml % L_check_moist_inc
  l_sl_bc_correction    = my_nml % l_sl_bc_correction
  l_viscosity    = my_nml % l_viscosity  

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_dyn
#endif

END MODULE  dynamics_input_mod
