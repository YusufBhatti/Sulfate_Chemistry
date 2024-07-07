! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          This module contains the routines for the diabatic 
!          tracer scheme. diab_tracer_sav request the tendencies to be 
!          saved, diab_tracer_slow and diab_tracer_fast adds the new
!          tendencies slow and physics tendencies into the tracers
!         

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: the variable flag controls where in atm_step is called
!                   and what functions it do.
!
MODULE diab_tracer_mod

USE atm_fields_bounds_mod,   ONLY: tdims, tdims_s

USE atm_step_local,          ONLY: first_atmstep_call

USE free_tracers_inputs_mod, ONLY: l_diab_tr_bl,  l_diab_tr_rad,        &
                                   num_dtheta

USE cv_run_mod,              ONLY: l_param_conv

! DrHook parameters
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE

CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'DIAB_TRACER_MOD'

! For loop iteration
INTEGER :: i,j,k

CONTAINS

SUBROUTINE diab_tracer_sav ()

USE cloud_inputs_mod,            ONLY: i_cld_vn
USE pc2_constants_mod,           ONLY: i_cld_pc2

USE physics_tendencies_mod,  ONLY:                                      &
    l_retain_rad_tendencies, l_retain_mic_tendencies,                   &
    l_retain_conv_tendencies, l_retain_bl_tendencies,                   &
    l_retain_q_cl_bl_tendencies, l_retain_cld_tendencies,               &
    l_retain_slow_tendencies,l_retain_ph1_tendencies


IMPLICIT NONE
!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER  :: routinename='DIAB_TRACER_SAV'
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Initialize and compute num_dpv
IF (first_atmstep_call) THEN
  
! Activate SI-SL advection of tracers:
! variable num_dtheta sets the number to tracers to be advected.

! Hardwired value to 5:
! dtheta_0, dtheta_bl, dtheta_rad, dtheta_mic, dtheta_slow
  num_dtheta = 5

! add dtheta_conv if param is switched on 
  IF (l_param_conv)  num_dtheta = num_dtheta + 1
! add dtheta_cld if PC2 is on
  IF (i_cld_vn == i_cld_pc2)   num_dtheta = num_dtheta + 1
! Add one more tracer if rad split is on (SW + LW  to rad)
  IF (l_diab_tr_rad) num_dtheta = num_dtheta + 1
! Add one more tracer if BL split is on (mic + LH to bl)
  IF (l_diab_tr_bl)  num_dtheta = num_dtheta + 1
  
  ! Set up values for retain_tendencies flags
  l_retain_rad_tendencies  = .TRUE.
  l_retain_mic_tendencies  = .TRUE.
  l_retain_slow_tendencies = .TRUE.
  l_retain_ph1_tendencies  = .TRUE.
  IF (l_param_conv)  l_retain_conv_tendencies = .TRUE.
  l_retain_bl_tendencies   = .TRUE.
  IF (l_diab_tr_bl)  l_retain_q_cl_bl_tendencies = .TRUE.
  IF (i_cld_vn == i_cld_pc2) l_retain_cld_tendencies  = .TRUE.
  
 
END IF !end if over first_atmstep_call

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE diab_tracer_sav

! ----------------------------------------------

SUBROUTINE diab_tracer_slow (exner_theta_levels, dtheta_0, dtheta_rad,      &
                             dtheta_SW, dtheta_LW, dtheta_mic, dtheta_slow)

USE atm_fields_mod,          ONLY: theta

USE physics_tendencies_mod,  ONLY: dt_sw, dt_lw, dt_mic, dtheta_ph1

IMPLICIT NONE

REAL, INTENT (IN) ::                                                    &
! Exner pressure on theta levels (to convert T -> theta)
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::                                                  &
! diab. Increments from radiation parametrization
  dtheta_0    (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! diab. Increments from radiation parametrization
, dtheta_rad  (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! diab. Increments from SW radiation parametrization
, dtheta_SW   (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! diab. Increments from LW radiation parametrization
, dtheta_LW   (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! diab. Increments from microphysics parametrization
, dtheta_mic  (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! diab. Increments from slow physiscs (atmos_physics1)
, dtheta_slow (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER  :: routinename='DIAB_TRACER_SLOW'
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Save Initial theta for first timestep
IF (first_atmstep_call) THEN
  DO k = tdims%k_start, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_0(i,j,k) = theta(i,j,k)
      END DO
    END DO
  END DO
END IF ! End first_atmstep_call

IF (l_diab_tr_rad) THEN
!Add SW increments to SW tracer
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_SW(i,j,k) = dtheta_SW(i,j,k)  +                          &
                           dt_sw(i,j,k)/exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  
  ! Set level 0 as level 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_SW(i,j,0)=dtheta_SW(i,j,1)
    END DO
  END DO
  
!Add LW increments to LW tracer
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_LW(i,j,k) = dtheta_LW(i,j,k)  +                          &
                           dt_lw(i,j,k)/exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  
  ! Set level 0 as level 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_LW(i,j,0)=dtheta_LW(i,j,1)
    END DO
  END DO
  

ELSE 
!Add increments (SW + LW) to rad tracer
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_rad(i,j,k) = dtheta_rad(i,j,k)  +                       &
                           ( dt_sw(i,j,k) + dt_lw(i,j,k) ) /           &    
                            exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
 ! Set level 0 as level 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_rad(i,j,0)=dtheta_rad(i,j,1)
    END DO
  END DO
  
END IF


! Add increments to dtheta_mic
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_mic(i,j,k) = dtheta_mic(i,j,k) +                           &
                          dt_mic(i,j,k)/exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
! Set level 0 as level 1
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    dtheta_mic(i,j,0)=dtheta_mic(i,j,1)
  END DO
END DO

!Add slow physics increments to slow physics tracer
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_slow(i,j,k) = dtheta_slow(i,j,k)  +                        &
                           dtheta_ph1(i,j,k)
    END DO
  END DO
END DO
! Set level 0 as level 1
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    dtheta_slow(i,j,0)=dtheta_slow(i,j,1)
  END DO
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE diab_tracer_slow

! ----------------------------------------------
SUBROUTINE diab_tracer_fast (exner_theta_levels,                        &
                             dtheta_bl, dtheta_bl_mix, dtheta_bl_LH,    &
                             dtheta_conv)

USE physics_tendencies_mod,  ONLY: dt_conv, dt_bl, dq_cl_bl

! Factors to compute the BL LH contribution if requested
USE planet_constants_mod, ONLY: cp
USE water_constants_mod, ONLY: lc

IMPLICIT NONE

REAL, INTENT (IN) ::                                                    &
! Exner pressure on theta levels (to convert T -> theta)
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::                                                  &
! diab. Increments from BL parametrization
  dtheta_bl    (tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)                              &
! diab. Increments from BL mixing parametrization
, dtheta_bl_mix(tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)                              & 
! diab. Increments from LH BL parametrization
, dtheta_bl_LH (tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)                              &
! diab. Increments from convection parametrization
, dtheta_conv  (tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)

REAL, ALLOCATABLE ::                                                    &
  inc_theta_bl_LH(:,:,:)

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER  :: routinename='DIAB_TRACER_FAST'
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Add increments from convection
IF (l_param_conv) THEN
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_conv(i,j,k) = dtheta_conv(i,j,k) +                       &
                             dt_conv(i,j,k)/exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  ! Set level 0 as level 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_conv(i,j,0)=dtheta_conv(i,j,1)
    END DO
  END DO

END IF ! End first_atmstep_call

! Add increments from BL
IF (l_diab_tr_bl) THEN
  !Compute LH contribution, first allocate
  IF (.NOT. ALLOCATED(inc_theta_bl_LH))                                 &
    ALLOCATE ( inc_theta_bl_LH(tdims%i_start:tdims%i_end,               &
                               tdims%j_start:tdims%j_end,               &
                               1:tdims%k_end))
  
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        inc_theta_bl_LH(i,j,k) = lc/cp*dq_cl_bl(i,j,k) /                &
                                 exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
    
  ! Add Mixing contribution to tracer
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_bl_mix(i,j,k) = dtheta_bl_mix(i,j,k) +                   &
                               (dt_bl(i,j,k)/exner_theta_levels(i,j,k)- &
                                inc_theta_bl_LH(i,j,k) ) 
      END DO
    END DO
  END DO
  ! Set level 0 as level 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_bl_mix(i,j,0)=dtheta_bl_mix(i,j,1)
    END DO
  END DO
  
! LH contribution
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_bl_LH(i,j,k) = dtheta_bl_LH(i,j,k) +                     &
                              inc_theta_bl_LH(i,j,k)
      END DO
    END DO
  END DO
  ! Set level 0 as level 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_bl_LH(i,j,0)=dtheta_bl_LH(i,j,1)
    END DO
  END DO

ELSE 
  !Allocate a dummy array for inc_theta_bl_LH
  IF (.NOT. ALLOCATED(inc_theta_bl_LH))                                 &
    ALLOCATE ( inc_theta_bl_LH(1,1,1))
  inc_theta_bl_LH(:,:,:) = 0.0
  
  ! All contribution to BL combined
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dtheta_bl(i,j,k) = dtheta_bl(i,j,k) +                             &
                           dt_bl(i,j,k)/exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  ! Set level 0 as level 1
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_bl(i,j,0)=dtheta_bl(i,j,1)
    END DO
  END DO

END IF 

! deallocate inc_theta_bl_LH
IF (ALLOCATED(inc_theta_bl_LH)) DEALLOCATE(inc_theta_bl_LH)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE diab_tracer_fast

! ----------------------------------------------
SUBROUTINE diab_tracer_np1 (exner_theta_levels, dtheta_cld)

USE physics_tendencies_mod,  ONLY: dt_cld

IMPLICIT NONE

REAL, INTENT (IN) ::                                                    &
! Exner pressure on theta levels (to convert T -> theta)
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)


REAL, INTENT(INOUT) ::                                                  &
! diab. Increments from cloud rebalancing
  dtheta_cld   (tdims%i_start:tdims%i_end,                              &
                tdims%j_start:tdims%j_end,                              &
                tdims%k_start:tdims%k_end)
!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER  :: routinename='DIAB_TRACER_NP1'
! End of header
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Add increments to dtheta_cld
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_cld(i,j,k) = dtheta_cld(i,j,k) +                           &
                          dt_cld(i,j,k) / exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
! Set level 0 as level 1
DO j = tdims%j_start, tdims%j_end
  DO i = tdims%i_start, tdims%i_end
    dtheta_cld(i,j,0)=dtheta_cld(i,j,1)
  END DO
END DO

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE diab_tracer_np1

END MODULE diab_tracer_mod
