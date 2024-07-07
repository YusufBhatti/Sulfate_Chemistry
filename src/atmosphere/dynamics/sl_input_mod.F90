! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Input control for SL dynamics.

! Description:
!   Module containing input switches/settings as used by the
!   semi_lagrangian code and the check_run_sl routine
!   for logic checking the selected settings.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics

! Method:
!   Switches are initialised to false and read in from the
!   namelist. The module may then be used directly where the switches
!   are needed within the dynamics semi_lagrangian code.

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3

MODULE sl_input_mod

USE missing_data_mod, ONLY: rmdi, imdi
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE var_look_mod, ONLY:                    &
    max_look, halo_lam ,halo_phi,          &
    recip_dlam,recip_dphi,                 &
    look_lam,look_phi

USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength

IMPLICIT NONE

! -------------------------------
! Parameters not read in by namelist.
! -------------------------------
! Set up indices for arrays holding interpolation options,
! eg. high_order_scheme(moist_SL) contains the chosen option
!     for high-order interpolation of moisture.
! See /atmosphere/dynamics_advection/highos_mod.F90.
INTEGER, PARAMETER :: Number_SL_choices = 5
INTEGER, PARAMETER :: Theta_SL = 1
INTEGER, PARAMETER :: moist_SL = 2
INTEGER, PARAMETER :: Wind_SL = 3
INTEGER, PARAMETER :: rho_SL = 4
INTEGER, PARAMETER :: tracer_SL = 5

INTEGER, PARAMETER :: unselected = 0
INTEGER, PARAMETER :: standard = 1
INTEGER, PARAMETER :: accurate = 2


! ===============================================
LOGICAL :: L_2d_sl_geometry = .FALSE.         ! (To be retired)
LOGICAL :: L_sl_halo_reprod = .FALSE.         ! (To be retired)
INTEGER :: Instability_diagnostics   = imdi    ! (To be retired)
! ===============================================


! -------------------------------
! Items not read in by namelist.
! -------------------------------

LOGICAL :: L_conserv(Number_SL_choices) = .FALSE.
                                    ! T if conservation to be enforced.

LOGICAL :: L_mono(Number_SL_choices) = .FALSE.
                                    ! T if monotonicity to be enforced.
LOGICAL :: L_high(Number_SL_choices) = .FALSE.
                                    ! T if high order interpolation scheme
                                          ! to be used.
LOGICAL :: L_Ritchie_high = .FALSE. ! T if high order scheme to be
                                    ! used in Ritchie routine.
LOGICAL :: L_Ritchie_mono = .FALSE. ! T if monotone scheme to be used
                                    ! in Ritchie routine.

INTEGER :: thmono_levels = 0   ! levels for monotone fully-interp theta

INTEGER :: Depart_scheme = 1   ! which departure point scheme to use (Ritchie)

INTEGER , PARAMETER ::  quasi_cubic = 2 
                        !ECMWF monotone quasi-cubic interpolation

! pmf_identifier is the integer number associated with
! using the PMF scheme for monotonicity
!                       
INTEGER , PARAMETER ::   pmf_identifier = 3
INTEGER , PARAMETER ::  spmf_identifier = 4
                         

!-----------------------------------------------------------
! run_sl namelist items
! ----------------------------------------------------------

LOGICAL :: L_conserve_tracers = .FALSE. ! Run tracer advection code with
                                        ! conservation on

! select what type of moisture conservation is required.
INTEGER :: moisture_conservation = imdi

! a code saying which high monotone scheme to use.
INTEGER :: high_order_scheme(Number_SL_choices) = imdi

! a code saying which  monotone scheme to use.
INTEGER :: monotone_scheme(Number_SL_choices) = imdi

INTEGER :: Ritchie_high_order_scheme = imdi ! which high order scheme to use
INTEGER :: Ritchie_monotone_scheme   = imdi ! which monotone scheme to use
INTEGER :: Depart_order              = imdi ! how many iterations/term
                                            ! to use for Ritchie scheme
INTEGER :: interp_vertical_search_tol= imdi ! number of levels either side
                                            ! of default level to search.

REAL :: thmono_height = rmdi             ! top for monotone fully-interp theta

! Logicals to control use of Priestley conservation (correction) scheme
! for moisture fields and Classic tracers
LOGICAL :: L_priestley_correct_moist = .FALSE.
LOGICAL :: L_priestley_correct_tracers = .FALSE.

INTEGER :: tr_priestley_opt     = imdi    ! Version of Priestley Schmeme
INTEGER :: moist_priestley_opt  = imdi    ! 1= original, 2= optimised

LOGICAL :: L_tr_src_in_conserve    = .FALSE.  ! Whether 'Fast' physics changes
LOGICAL :: L_moist_src_in_conserve = .FALSE.  ! (sources) are considered

! Logical to control application of Priestley to thetav
LOGICAL :: L_priestley_correct_thetav = .FALSE.

! ----------------------
! Namelist
! ----------------------

NAMELIST/run_sl/ thmono_height,                                  &
 high_order_scheme, monotone_scheme,                             &
 Ritchie_high_order_scheme,Ritchie_monotone_scheme,              &
 Depart_order,interp_vertical_search_tol,                        &
 L_conserve_tracers,                                             &
 moisture_conservation,                                          &
 L_priestley_correct_moist, L_priestley_correct_tracers,         &
 tr_priestley_opt, moist_priestley_opt,                          &
 L_tr_src_in_conserve, L_moist_src_in_conserve,                  &
 L_priestley_correct_thetav

!DrHook-related parameters
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER, PRIVATE :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SL_INPUT_MOD'

CONTAINS

SUBROUTINE check_run_sl()

! Description:
!   Subroutine to apply logic checks and setup variables based on the
!   options selected in the run_sl namelist.

USE highos_mod, ONLY:       &
    hLag3_vHerm3_d2, hLag3_vHerm3_d4, ECMWF_quasiCubic, &
    ECMWF_mono_quasiCubic, hQuasiCubic_vQuintic

IMPLICIT NONE

INTEGER :: errcode
CHARACTER (LEN=errormessagelength) :: cmessage

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='CHECK_RUN_SL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------
! SL input setup
!-------------------------------------------------------
! moisture conservation setup
!-------------------------------------------------------
IF (moisture_conservation==accurate) THEN
  L_conserv(moist_SL) = .TRUE.
ELSE IF (moisture_conservation==standard) THEN
  L_conserv(moist_SL) = .TRUE.
END IF

IF ( L_priestley_correct_moist .AND.                               &
     (moist_priestley_opt < 1 .OR. moist_priestley_opt > 2) ) THEN
  errcode = 1
  WRITE(cmessage,'(A)') 'RUN_SL: Priestley scheme option '//     &
     'not correctly set (moist_priestley_opt should be 1 or 2)'
  CALL ereport ('CHECK_RUN_SL',errcode,cmessage)
END IF

!-------------------------------------------------------
! monotone scheme setup
!-------------------------------------------------------
IF (monotone_scheme(Theta_SL)/=unselected) THEN
  L_mono(Theta_SL) = .TRUE.
END IF
IF (monotone_scheme(moist_SL)/=unselected) THEN
  L_mono(moist_SL) = .TRUE.
END IF
IF (monotone_scheme(Wind_SL)/=unselected) THEN
  L_mono(Wind_SL) = .TRUE.
END IF
IF (monotone_scheme(rho_SL)/=unselected) THEN
  L_mono(rho_SL) = .TRUE.
END IF
IF (monotone_scheme(tracer_SL)/=unselected) THEN
  L_mono(tracer_SL) = .TRUE.         
END IF


IF (ANY (high_order_scheme == ECMWF_quasiCubic) .OR.       &
    ANY (high_order_scheme == ECMWF_mono_quasiCubic) .OR.  &
    ANY (high_order_scheme == hQuasiCubic_vQuintic) ) THEN
  cmessage =                                                 &
  'An invalid high order scheme has been selected for ENDGame' 
  errcode = 26
  CALL ereport('CHECK_RUN_SL', errcode, cmessage)
END IF

! ECMWF monotone quasi-cubic interpolation is only valid for New Dynamics.
IF (ANY (monotone_scheme == quasi_cubic) ) THEN  
  cmessage =                                                 &
  'An invalid monotone scheme has been selected for ENDGame' 
  errcode = 27
  CALL ereport('CHECK_RUN_SL', errcode, cmessage)
END IF

!-------------------------------------------------------
! high order scheme setup
!-------------------------------------------------------
IF (high_order_scheme(Theta_SL)/=unselected) THEN
  L_high(Theta_SL) = .TRUE.
END IF
IF (high_order_scheme(moist_SL)/=unselected) THEN
  L_high(moist_SL) = .TRUE.
END IF
IF (high_order_scheme(Wind_SL)/=unselected) THEN
  L_high(Wind_SL) = .TRUE.
END IF
IF (high_order_scheme(rho_SL)/=unselected) THEN
  L_high(rho_SL) = .TRUE.
END IF
IF (high_order_scheme(tracer_SL)/=unselected) THEN
  L_high(tracer_SL) = .TRUE.
END IF
!-------------------------------------------------------
! Ritchie scheme setup
!-------------------------------------------------------
IF (Ritchie_high_order_scheme/=unselected) THEN
  L_Ritchie_high = .TRUE.
END IF
IF (Ritchie_monotone_scheme/=unselected) THEN
  L_Ritchie_mono = .TRUE.
END IF

!-------------------------------------------------------
! Tracer - Priestley conservation scheme setup
!-------------------------------------------------------
L_conserv(tracer_SL) = l_conserve_tracers

IF ( L_priestley_correct_tracers .AND.                             &
     (tr_priestley_opt < 1 .OR. tr_priestley_opt > 2) ) THEN
  errcode = 2
  WRITE(cmessage,'(A)') 'RUN_SL: Priestley scheme option '//     &
     'not correctly set (tr_priestley_opt should be 1 or 2)'
  CALL ereport ('CHECK_RUN_SL',errcode,cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE check_run_sl


SUBROUTINE print_nlist_run_sl()
USE umPrintMgr, ONLY: umPrint
IMPLICIT NONE
CHARACTER(LEN=50000) :: lineBuffer
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_NLIST_RUN_SL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL umPrint('Contents of namelist run_sl', &
    src='sl_input_mod')

WRITE(lineBuffer,*)' L_conserv = ',L_conserv
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_mono = ',L_mono
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_high = ',L_high
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_Ritchie_high = ',L_Ritchie_high
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_Ritchie_mono = ',L_Ritchie_mono
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' thmono_height = ',thmono_height
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' high_order_scheme = ',high_order_scheme
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' monotone_scheme = ',monotone_scheme
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' Ritchie_high_order_scheme = ',Ritchie_high_order_scheme
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' Ritchie_monotone_scheme = ',Ritchie_monotone_scheme
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' Depart_scheme = ',Depart_scheme
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' Depart_order = ',Depart_order
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' interp_vertical_search_tol = ',interp_vertical_search_tol
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_conserve_tracers = ',L_conserve_tracers
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_priestley_correct_moist = ',                 &
                                         L_priestley_correct_moist
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_priestley_correct_tracers=',                 &
                                       L_priestley_correct_tracers
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' tr_priestley_opt = ',tr_priestley_opt
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' moist_priestley_opt = ',moist_priestley_opt
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_tr_src_in_conserve = ',L_tr_src_in_conserve
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,*)' L_moist_src_in_conserve = ',L_moist_src_in_conserve
CALL umPrint(lineBuffer,src='sl_input_mod')
WRITE(lineBuffer,'(A,L1)') ' L_priestley_correct_thetav =',             &
                                              L_priestley_correct_thetav
CALL umPrint(lineBuffer,src='sl_input_mod')

CALL umPrint('- - - - - - end of namelist - - - - - -', &
    src='sl_input_mod')

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE print_nlist_run_sl

#if !defined(LFRIC)
SUBROUTINE read_nml_run_sl(unit_in)

USE um_parcore, ONLY: mype

USE check_iostat_mod, ONLY: check_iostat

USE setup_namelist, ONLY: setup_nml_type

IMPLICIT NONE

INTEGER, INTENT(IN) ::unit_in
INTEGER :: my_comm
INTEGER :: mpl_nml_type
INTEGER :: ErrorStatus
INTEGER :: icode
CHARACTER(LEN=errormessagelength) :: iomessage
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='READ_NML_RUN_SL'

! set number of each type of variable in my_namelist type
INTEGER, PARAMETER :: no_of_types = 3
INTEGER, PARAMETER :: n_int = 7 + 2 * Number_SL_choices
INTEGER, PARAMETER :: n_real = 1
INTEGER, PARAMETER :: n_log = 6

TYPE my_namelist
  SEQUENCE
  INTEGER :: high_order_scheme(Number_SL_choices)
  INTEGER :: monotone_scheme(Number_SL_choices)
  INTEGER :: Ritchie_high_order_scheme
  INTEGER :: Ritchie_monotone_scheme
  INTEGER :: Depart_order
  INTEGER :: interp_vertical_search_tol
  INTEGER :: moisture_conservation
  INTEGER :: tr_priestley_opt
  INTEGER :: moist_priestley_opt
  REAL :: thmono_height
  LOGICAL :: L_conserve_tracers
  LOGICAL :: L_priestley_correct_moist
  LOGICAL :: L_moist_src_in_conserve
  LOGICAL :: L_priestley_correct_tracers
  LOGICAL :: L_tr_src_in_conserve
  LOGICAL :: L_priestley_correct_thetav
END TYPE my_namelist

TYPE (my_namelist) :: my_nml

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

CALL gc_get_communicator(my_comm, icode)

CALL setup_nml_type(no_of_types, mpl_nml_type, n_int_in=n_int,     &
                    n_real_in=n_real, n_log_in=n_log)

IF (mype == 0) THEN

  READ(UNIT=unit_in, NML=run_sl, IOSTAT=ErrorStatus, IOMSG=iomessage)
  CALL check_iostat(errorstatus, "namelist RUN_SL", iomessage)

  my_nml % high_order_scheme          = high_order_scheme
  my_nml % monotone_scheme            = monotone_scheme
  my_nml % Ritchie_high_order_scheme  = Ritchie_high_order_scheme
  my_nml % Ritchie_monotone_scheme    = Ritchie_monotone_scheme
  my_nml % Depart_order               = Depart_order
  my_nml % interp_vertical_search_tol = interp_vertical_search_tol
  my_nml % moisture_conservation      = moisture_conservation
  my_nml % tr_priestley_opt           = tr_priestley_opt
  my_nml % moist_priestley_opt        = moist_priestley_opt
  ! end of integers
  my_nml % thmono_height = thmono_height
  ! end of reals
  my_nml % L_conserve_tracers = L_conserve_tracers
  my_nml % L_priestley_correct_moist   = L_priestley_correct_moist
  my_nml % L_moist_src_in_conserve     = L_moist_src_in_conserve
  my_nml % L_priestley_correct_tracers = L_priestley_correct_tracers
  my_nml % L_tr_src_in_conserve        = L_tr_src_in_conserve
  my_nml % L_priestley_correct_thetav  = L_priestley_correct_thetav

END IF

CALL mpl_bcast(my_nml,1,mpl_nml_type,0,my_comm,icode)

IF (mype /= 0) THEN

  high_order_scheme          = my_nml % high_order_scheme
  monotone_scheme            = my_nml % monotone_scheme
  Ritchie_high_order_scheme  = my_nml % Ritchie_high_order_scheme
  Ritchie_monotone_scheme    = my_nml % Ritchie_monotone_scheme
  Depart_order               = my_nml % Depart_order
  interp_vertical_search_tol = my_nml % interp_vertical_search_tol
  moisture_conservation      = my_nml % moisture_conservation
  tr_priestley_opt           = my_nml % tr_priestley_opt
  moist_priestley_opt        = my_nml % moist_priestley_opt
  ! end of integers
  thmono_height = my_nml % thmono_height
  ! end of reals
  L_conserve_tracers = my_nml % L_conserve_tracers
  L_priestley_correct_moist = my_nml % L_priestley_correct_moist
  L_moist_src_in_conserve   = my_nml % L_moist_src_in_conserve
  L_priestley_correct_tracers = my_nml % L_priestley_correct_tracers
  L_tr_src_in_conserve   = my_nml % L_tr_src_in_conserve
  L_priestley_correct_thetav = my_nml % L_priestley_correct_thetav

END IF

CALL mpl_type_free(mpl_nml_type,icode)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE read_nml_run_sl
#endif

END MODULE  sl_input_mod
