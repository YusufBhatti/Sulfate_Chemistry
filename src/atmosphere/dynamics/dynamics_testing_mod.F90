! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for dynamics.

! Description:
!   Module containing input switches/settings and the check_run_dyntest
!   routine for logic checking the selected settings. Used by the
!   dynamics to enable testing and some idealised running.
!   Each switch would usually not be used by the everyday user
!   and is defaulted to the common setting.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics
!
! Method:
!   Switches are initialised to false and read in from the
!   namelists. The module may then be used directly where the switches
!   are needed within the dynamics code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE dynamics_testing_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE gcr_input_mod, ONLY:                                                     &
    GCR_Diagnostics, GCR_its_avg_step

USE var_end_mod, ONLY:                                                       &
    lambda_p_end,phi_p_end,lambda_u_end,phi_v_end,                           &
    dlambda_p_end,dphi_p_end,dlambda_u_end,dphi_v_end

USE var_input_mod, ONLY:                                                     &
    lam_var, phi_var,var_ratio,lam_ratio,phi_ratio,                          &
    phi_frac, lam_frac

USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! ===============================================
LOGICAL :: L_adjust_wet   = .FALSE.           ! (To be retired)
LOGICAL :: L_free_slip    = .FALSE.           ! (To be retired)
REAL    :: uv_limit       = rmdi               ! (To be retired)  
INTEGER :: trap_option    = imdi               ! (To be retired)  
LOGICAL :: L_trap_theta   = .FALSE.           ! (To be retired)  
LOGICAL :: L_trap_w       = .FALSE.           ! (To be retired) 
REAL    :: Cw_max         = rmdi               ! (To be retired) 
INTEGER :: Cw_test_lev    = imdi               ! (To be retired) 
REAL    :: max_thinc      = rmdi               ! (To be retired) 
! ===============================================

! the following general users need not alter.
! power users may like to alter these via the namelist run_dyntest

LOGICAL :: L_Physics           = .FALSE.  ! T: physics to be included
LOGICAL :: L_dynamics_only     = .FALSE.  ! T: run without physics
LOGICAL :: L_Run_With_Physics2 = .FALSE.  ! T: physics2 to be included
LOGICAL :: L_exclude_Physics2  = .FALSE.  ! T: physics2 to be excluded
LOGICAL :: L_perturb_IC_theta  = .FALSE.  ! T: perturb theta on ts1
LOGICAL :: L_Backwards         = .FALSE.  ! F Integrate backwards without
                                        ! physics
LOGICAL :: L_dry               = .FALSE.  ! T run with no moisture

LOGICAL :: L_idealised_data    = .FALSE.  ! T run idealised problem

LOGICAL :: L_hydrostatic_EG    = .FALSE.  ! choose to run hydrostatic in EG.

LOGICAL :: L_prognostic_level0 = .FALSE.  ! F: Use constant extrapolation
                                          !    in vertical for thetavd and
                                          !    moisture on level 0.

! trapping
LOGICAL :: L_trap_uv    = .FALSE.  ! .true. trap excessive u,v 

! Type of problem to be solved
INTEGER :: problem_number = imdi
! Seed for random numbers used to perturb initial condition
! cf. L_perturb_IC_theta
INTEGER :: IntRand_seed = imdi

NAMELIST/RUN_Dyntest/                                                        &
    L_dynamics_only, L_Backwards, L_hydrostatic_EG,                          &
    L_exclude_Physics2, L_perturb_IC_theta,                                  &
    L_dry,                                                                   &
    L_trap_uv, L_idealised_data,                                             &
    GCR_Diagnostics, GCR_its_avg_step,                                       &
    lambda_p_end,  phi_p_end, lambda_u_end,  phi_v_end,                      &
    dlambda_p_end,  dphi_p_end, dlambda_u_end,  dphi_v_end,                  &
    lam_var, phi_var,                                                        &
    var_ratio, lam_ratio, phi_ratio, lam_frac, phi_frac,                     &
    problem_number, IntRand_seed

! DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DYNAMICS_TESTING_MOD'

CONTAINS

SUBROUTINE check_run_dyntest()

! Description:
!   Subroutine to apply logic checks and set variables based on the
!   options selected in the run_dyntest namelist.

USE dynamics_input_mod,     ONLY: Ih
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE
CHARACTER (LEN=errormessagelength)            :: cmessage
INTEGER                                       :: icode

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_RUN_DYNTEST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ---------------------------------------------
! Dynamics logic control - New Dynamics and ENDGAME
! ---------------------------------------------
IF (l_dynamics_only) THEN
  l_physics           = .FALSE.
ELSE
  l_physics           = .TRUE.
END IF

IF (l_exclude_physics2) THEN
  l_run_with_physics2 = .FALSE.
ELSE
  l_run_with_physics2 = .TRUE.
END IF

IF (L_Backwards) L_Physics = .FALSE.

! ---------------------------------------------
! Dynamics logic control - ENDGAME Only
! ---------------------------------------------

IF (l_hydrostatic_EG) THEN
  Ih =0.0
ELSE
  Ih =1.0
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_dyntest

SUBROUTINE print_nlist_run_dyntest()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_DYNTEST'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_dyntest',                           &
    src='dynamics_testing_mod')

WRITE(lineBuffer,'(''L_dynamics_only = '',L1)') L_dynamics_only
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''L_Backwards = '',L1)') L_Backwards
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''L_hydrostatic_EG = '',L1)') L_hydrostatic_EG
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''L_exclude_Physics2 = '',L1)') L_exclude_Physics2
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''L_perturb_IC_theta = '',L1)') L_perturb_IC_theta
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''L_dry = '',L1)') L_dry
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''L_trap_uv = '',L1)') L_trap_uv
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''L_idealised_data = '',L1)') L_idealised_data
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''GCR_Diagnostics = '',I0)') GCR_Diagnostics
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''GCR_its_avg_step ='',3(1X,I0))') GCR_its_avg_step
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''lambda_p_end = '',E15.7)') lambda_p_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''phi_p_end = '',E15.7)') phi_p_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''lambda_u_end = '',E15.7)') lambda_u_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''phi_v_end = '',E15.7)') phi_v_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''dlambda_p_end = '',E15.7)') dlambda_p_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''dphi_p_end = '',E15.7)') dphi_p_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''dlambda_u_end = '',E15.7)') dlambda_u_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''dphi_v_end = '',E15.7)') dphi_v_end
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''lam_var = '',I0)') lam_var
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''phi_var = '',I0)') phi_var
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''var_ratio = '',E15.7)') var_ratio
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''lam_ratio = '',E15.7)') lam_ratio
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''phi_ratio = '',E15.7)') phi_ratio
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''lam_frac = '',E15.7)') lam_frac
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''phi_frac = '',E15.7)') phi_frac
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''problem_number = '',I0)') problem_number
CALL umPrint(lineBuffer,src='dynamics_testing_mod')
WRITE(lineBuffer,'(''IntRand_seed = '',I0)') IntRand_seed
CALL umPrint(lineBuffer,src='dynamics_input_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -',                    &
    src='dynamics_testing_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_dyntest

#if !defined(LFRIC)
SUBROUTINE read_nml_run_dyntest(unit_in)

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

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_DYNTEST'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 5 + 3
INTEGER, PARAMETER :: n_real = 13
INTEGER, PARAMETER :: n_log = 8

TYPE my_namelist
  SEQUENCE
  INTEGER :: GCR_Diagnostics
  INTEGER :: GCR_its_avg_step(3)
  INTEGER :: lam_var
  INTEGER :: phi_var
  INTEGER :: problem_number
  INTEGER :: IntRand_seed
  REAL :: lambda_p_end
  REAL :: phi_p_end
  REAL :: lambda_u_end
  REAL :: phi_v_end
  REAL :: dlambda_p_end
  REAL :: dphi_p_end
  REAL :: dlambda_u_end
  REAL :: dphi_v_end
  REAL :: var_ratio
  REAL :: lam_ratio
  REAL :: phi_ratio
  REAL :: lam_frac
  REAL :: phi_frac
  LOGICAL :: L_dynamics_only
  LOGICAL :: L_Backwards
  LOGICAL :: L_hydrostatic_EG
  LOGICAL :: L_exclude_Physics2
  LOGICAL :: L_perturb_IC_theta
  LOGICAL :: L_dry
  LOGICAL :: L_trap_uv
  LOGICAL :: L_idealised_data
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=RUN_Dyntest, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_Dyntest", iomessage)

  my_nml % GCR_Diagnostics  = GCR_Diagnostics
  my_nml % GCR_its_avg_step = GCR_its_avg_step
  my_nml % lam_var          = lam_var
  my_nml % phi_var          = phi_var
  my_nml % problem_number   = problem_number
  my_nml % IntRand_seed     = IntRand_seed
  ! end of integers
  my_nml % lambda_p_end  = lambda_p_end
  my_nml % phi_p_end     = phi_p_end
  my_nml % lambda_u_end  = lambda_u_end
  my_nml % phi_v_end     = phi_v_end
  my_nml % dlambda_p_end = dlambda_p_end
  my_nml % dphi_p_end    = dphi_p_end
  my_nml % dlambda_u_end = dlambda_u_end
  my_nml % dphi_v_end    = dphi_v_end
  my_nml % var_ratio     = var_ratio
  my_nml % lam_ratio     = lam_ratio
  my_nml % phi_ratio     = phi_ratio
  my_nml % lam_frac      = lam_frac
  my_nml % phi_frac      = phi_frac
  ! end of reals
  my_nml % L_dynamics_only    = L_dynamics_only
  my_nml % L_Backwards        = L_Backwards
  my_nml % L_hydrostatic_EG   = L_hydrostatic_EG
  my_nml % L_exclude_Physics2 = L_exclude_Physics2
  my_nml % L_perturb_IC_theta = L_perturb_IC_theta
  my_nml % L_dry              = L_dry
  my_nml % L_trap_uv          = L_trap_uv
  my_nml % L_idealised_data   = L_idealised_data

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  GCR_Diagnostics  = my_nml % GCR_Diagnostics
  GCR_its_avg_step = my_nml % GCR_its_avg_step
  lam_var          = my_nml % lam_var
  phi_var          = my_nml % phi_var
  problem_number   = my_nml % problem_number
  IntRand_seed     = my_nml % IntRand_seed
  ! end of integers
  lambda_p_end  = my_nml % lambda_p_end
  phi_p_end     = my_nml % phi_p_end
  lambda_u_end  = my_nml % lambda_u_end
  phi_v_end     = my_nml % phi_v_end
  dlambda_p_end = my_nml % dlambda_p_end
  dphi_p_end    = my_nml % dphi_p_end
  dlambda_u_end = my_nml % dlambda_u_end
  dphi_v_end    = my_nml % dphi_v_end
  var_ratio     = my_nml % var_ratio
  lam_ratio     = my_nml % lam_ratio
  phi_ratio     = my_nml % phi_ratio
  lam_frac      = my_nml % lam_frac
  phi_frac      = my_nml % phi_frac
  ! end of reals
  L_dynamics_only    = my_nml % L_dynamics_only
  L_Backwards        = my_nml % L_Backwards
  L_hydrostatic_EG   = my_nml % L_hydrostatic_EG
  L_exclude_Physics2 = my_nml % L_exclude_Physics2
  L_perturb_IC_theta = my_nml % L_perturb_IC_theta
  L_dry              = my_nml % L_dry
  L_trap_uv          = my_nml % L_trap_uv
  L_idealised_data   = my_nml % L_idealised_data

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_dyntest
#endif

END MODULE dynamics_testing_mod
