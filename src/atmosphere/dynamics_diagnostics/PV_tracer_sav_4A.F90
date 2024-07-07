! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          Initializes, saves and computes important variables
!          for the PV-tracer, such as initial density, pv_0 or pv_1
!

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: the variable flag controls where in atm_step is called
!                   and what functions it do.
!
MODULE PV_tracer_sav_4A_mod

USE free_tracers_inputs_mod, ONLY:                                      &
    l_pv_tracer, l_pv_dyn, num_dpv, l_pv_split_rad, l_pv_adv_only

! error handling
USE ereport_mod, ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE umPrintMgr, ONLY: newline

! DrHook parameters
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

REAL, ALLOCATABLE ::                                                    &
  rho_0(:,:,:)                                                          &
, pv_0(:,:,:)                                                           &
, pv_1(:,:,:)

REAL, PARAMETER :: si_to_pvunits = 1.0e6

CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'PV_TRACER_SAV_4A_MOD'

CONTAINS

SUBROUTINE PV_tracer_sav_4A ()

! Use Bounds for arrays
USE atm_fields_bounds_mod,       ONLY: tdims, pdims_s,                  & 
                                       vdims_s, udims_s

USE atm_step_local,              ONLY: first_atmstep_call

USE cloud_inputs_mod,            ONLY: i_cld_vn

USE pc2_constants_mod,           ONLY: i_cld_pc2

!Call routines to set num_dpv
USE stochastic_physics_run_mod,  ONLY: l_skeb2, l_spt
USE IAU_mod,                     ONLY: l_iau
USE nudging_input_mod,           ONLY: l_nudging

USE cv_run_mod,                  ONLY: l_param_conv

USE horiz_grid_mod, ONLY: cartesian_grid

USE model_domain_mod, ONLY: l_regular

USE physics_tendencies_mod,      ONLY:                                  &
    l_retain_slow_tendencies,                                           &
    l_retain_rad_tendencies, l_retain_mic_tendencies,                   &
    l_retain_gwd_tendencies, l_retain_ph1_tendencies,                   &
    l_retain_conv_tendencies, l_retain_conv_all_tendencies,             &
    l_retain_conv_mom_tendencies, l_retain_bl_tendencies,               &
    l_retain_stph_tendencies, l_retain_cld_tendencies,                  &
    l_retain_nud_tendencies

IMPLICIT NONE

!Local variables
INTEGER :: i,j,k
           ! indexing for DO loops

! Error handling
CHARACTER(LEN=errormessagelength) :: cmessage      ! out error message
CHARACTER(LEN=*), PARAMETER  :: routinename='PV_TRACER_SAV_4A'
INTEGER :: icode    ! error code

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Initialize and compute num_dpv
IF (first_atmstep_call) THEN

  icode = 0
  IF (.NOT. l_regular) THEN
     WRITE(cmessage,'(A)')                                           &
    'PV-tracer code is not compatible with'//               newline//&
    'variable resolution models. '//                        newline//&
    'Please set l_pv_tracer = false in' //                  newline//&
    'free_tracers_inputs namelist (sect. 33).'
    icode   = 1
    CALL ereport(routinename, icode, cmessage)
  END IF

  IF (cartesian_grid) THEN
     WRITE(cmessage,'(A)')                                           &
    'PV-tracer code is not compatible with'//               newline//&
    'coordinate systems other than '//                      newline//&
    'spherical coordinates'
    icode   = 1
    CALL ereport(routinename, icode, cmessage)
  END IF

  
! Activate SI-SL advection of tracers:
! variable num_dpv sets the number to tracers to be advected.

! Hardwired value to 7: 
! dPV_rad, dPV_mic, dPV_gwd, dPV_ph1, dPV_bl, dPV_tot, PV_0
  num_dpv = 7

! add dPV_conv if param is switched on 
  IF (l_param_conv)            num_dpv = num_dpv + 1
! add dPV_cld if PC2 is on
  IF (i_cld_vn == i_cld_pc2)   num_dpv = num_dpv + 1

! Add dPV_stph if SKEB2 or SPT are on
  IF (l_skeb2 .OR. l_spt) num_dpv = num_dpv + 1
! Add dPV_iau if IAU is on  
  IF (l_iau)              num_dpv = num_dpv + 1
! Add dPV_nud if nudging scheme is on    
  IF (l_nudging)          num_dpv = num_dpv + 1
! If dynamics PV are requested, add
! dPV_adv, dPV_sol, dPV_mass, PV_1
  IF (l_pv_dyn)           num_dpv = num_dpv + 4 
! If split rad requested, replace dPV_rad by dPV_sw and dPV_lw (add one more)
  IF (l_pv_split_rad)     num_dpv = num_dpv + 1
! If adv only passive tracer req add adv_only_PV
  IF (l_pv_adv_only)      num_dpv = num_dpv + 1 !adv only tracer

  ! Set up values for retain_tendencies flags
  l_retain_slow_tendencies = .TRUE.
  l_retain_rad_tendencies  = .TRUE.
  l_retain_mic_tendencies  = .TRUE.
  l_retain_gwd_tendencies  = .TRUE.
  l_retain_ph1_tendencies  = .TRUE.
  IF (l_param_conv) THEN
    l_retain_conv_tendencies = .TRUE.
    l_retain_conv_mom_tendencies = .TRUE.
    l_retain_conv_all_tendencies = .TRUE.
  END IF
  l_retain_bl_tendencies   = .TRUE.
  IF (i_cld_vn == i_cld_pc2) l_retain_cld_tendencies  = .TRUE.

  IF (l_skeb2 .OR. l_spt) l_retain_stph_tendencies =.TRUE.
  IF (l_nudging) l_retain_nud_tendencies =.TRUE.

END IF !end if over first_atmstep_call

! save rho_0 array
IF (.NOT. ALLOCATED(rho_0))                                           &
  ALLOCATE(rho_0(pdims_s%i_start:pdims_s%i_end,                       &
                 pdims_s%j_start:pdims_s%j_end,                       &
                 pdims_s%k_start:pdims_s%k_end))
 rho_0 = 0.0
 ! save pv_0 array
IF (.NOT. ALLOCATED(pv_0))                                            &
  ALLOCATE(pv_0(tdims%i_start:tdims%i_end,                            &
                tdims%j_start:tdims%j_end,                            &
                tdims%k_start:tdims%k_end))
 pv_0 = 0.0

IF (l_pv_dyn .OR. l_iau) THEN
  ! save pv_1 array
  IF (.NOT. ALLOCATED(pv_1))                                          &
    ALLOCATE(pv_1(tdims%i_start:tdims%i_end,                          &
                  tdims%j_start:tdims%j_end,                          &
                  tdims%k_start:tdims%k_end))
    pv_1 = 0.0
END IF

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE PV_tracer_sav_4A

! ------------------------------------------------------------------------

SUBROUTINE PV_tracer_dealloc()

USE IAU_mod,                     ONLY: l_iau

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER  :: routinename='PV_TRACER_DEALLOC'
!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ( (l_pv_dyn  .OR. l_iau) .AND. ALLOCATED(pv_1))  DEALLOCATE(pv_1)

IF (ALLOCATED(pv_0))   DEALLOCATE(pv_0)
IF (ALLOCATED(rho_0))   DEALLOCATE(rho_0)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE PV_tracer_dealloc

END MODULE PV_tracer_sav_4A_mod
