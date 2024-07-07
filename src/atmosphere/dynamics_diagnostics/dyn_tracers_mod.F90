! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          This module contains the routines to merge the diabatic
!          tracer and PV tracer schemes.
!

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: the variable flag controls where in atm_step is called
!                   and what functions it do.


MODULE dyn_tracers_mod

!dims
USE atm_fields_bounds_mod,  ONLY: tdims, tdims_s, pdims_s, udims_s,     &
                                  vdims_s
! PV-tracer
USE PV_tracer_sav_4A_mod,   ONLY: PV_tracer_sav_4A
USE PV_tracer_slow_mod,     ONLY: PV_tracer_slow
USE eg_sl_PV_trac_mod,      ONLY: eg_sl_PV_trac
USE PV_tracer_fast_mod,     ONLY: PV_tracer_fast
USE PV_tracer_np1_mod,      ONLY: PV_tracer_np1
USE PV_tracer_end_mod,      ONLY: PV_tracer_end

!Diabatic tracer
USE eg_sl_diab_trac_mod,    ONLY: eg_sl_diab_trac
USE diab_tracer_mod,        ONLY: diab_tracer_sav, diab_tracer_slow,    &
                                  diab_tracer_fast,diab_tracer_np1

USE free_tracers_inputs_mod,ONLY: l_pv_tracer, l_diab_tracer

! DrHook parameters
USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim


IMPLICIT NONE

CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'DYN_TRACERS_MOD'

! For loop iteration
INTEGER :: i,j,k

CONTAINS

SUBROUTINE dyn_tr_sav()

IMPLICIT NONE
CHARACTER(LEN=*), PARAMETER  :: routinename='DYN_TR_SAV'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_pv_tracer) CALL PV_tracer_sav_4A()
IF (l_diab_tracer) CALL diab_tracer_sav()

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dyn_tr_sav

!-----------------------------------------------------------------

SUBROUTINE dyn_tr_slow( rho, exner_theta_levels,                        &
                        ! PV-tracers
                        dPV_rad, dPV_sw, dPV_lw, dPV_mic, dPV_gwd,      &
                        dPV_ph1, adv_only_PV,                           &
                        ! diab-tracers
                        dtheta_0, dtheta_rad, dtheta_SW, dtheta_LW,     &
                        dtheta_mic,dtheta_slow)

IMPLICIT NONE

REAL, INTENT (IN) ::                                                    &
! main prog. Rho field
  rho          (pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end)                          &

! Exner pressure on theta levels (to convert T -> theta)
, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::                                                  &
! PV Increments from radiation parametrization
  dPV_rad     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from SW radiation parametrization
, dPV_sw      (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from LW radiation parametrization
, dPV_lw      (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from microphysics parametrization
, dPV_mic     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from Gravity Wave Drag parametrization
, dPV_gwd     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from All slow physics
, dPV_ph1     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV at timestep 0 advected
, adv_only_PV (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! diabatic advection only tracer, a.k.a. theta at T+0 advected
, dtheta_0    (tdims%i_start:tdims%i_end,                               &
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
! diab. Increments from all slow physics (atmos_physics1)
, dtheta_slow (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)

CHARACTER(LEN=*), PARAMETER  :: routinename='DYN_TR_SLOW'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!End of header


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! PV TRACER:  Compute slow dPV increments
IF (l_pv_tracer) CALL PV_tracer_slow( rho, exner_theta_levels,          &
                                      dPV_rad, dPV_sw, dPV_lw,          &
                                      dPV_mic, dPV_gwd, dPV_ph1,        &
                                      adv_only_PV)

! diabatic tracers for slow physics
IF (l_diab_tracer) CALL diab_tracer_slow( exner_theta_levels,           &
                                          dtheta_0, dtheta_rad,         &
                                          dtheta_SW, dtheta_LW,         &
                                          dtheta_mic, dtheta_slow)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dyn_tr_slow

!-----------------------------------------------------------------
SUBROUTINE eg_sl_dyn_tr(                                              &
             row_length, rows, n_rows, model_levels,                  &
             g_i_pe, high_order_scheme, monotone_scheme,              &
             l_high, l_mono,                                          &
             dPV_rad, dPV_sw, dPV_lw,                                 &
             dPV_mic, dPV_gwd, dPV_ph1,                               &
             dPV_conv, dPV_bl, dPV_stph, dPV_cld,                     &
             dPV_iau, dPV_nud, dPV_tot,                               &
             dPV_adv, dPV_sol, dPV_mass, adv_only_PV,                 &
             dtheta_0, dtheta_bl, dtheta_bl_mix, dtheta_bl_LH,        &
             dtheta_conv, dtheta_mic, dtheta_rad,                     &
             dtheta_SW, dtheta_LW, dtheta_slow, dtheta_cld,           &
             errorstatus)

USE um_parvars,           ONLY: halo_i
USE nlsizes_namelist_mod, ONLY: global_row_length


IMPLICIT NONE
! Arguments with Intent IN. ie: Input variables.
INTEGER, INTENT(IN) ::                                                &
! number of points on a row
  row_length                                                          &
! number of rows.
, rows                                                                &
! number of v-rows.
, n_rows                                                              &
! number of model levels.
, model_levels                                                        &
! processor on my processor-row holding a given value in i direction
, g_i_pe(1-halo_i:global_row_length+halo_i)                           &
! a code saying which high order scheme to use for interpolation
, high_order_scheme                                                   &
! a code saying which monotone scheme to use for interpolation
, monotone_scheme

LOGICAL, INTENT(IN) ::                                                &
! True, if high order interpolation required
  l_high                                                              &
 ! True, if interpolation required to be monotone
, l_mono

! Prognostic fields for the diabatic-tracer scheme
REAL, INTENT(INOUT) ::                                                &
    dPV_rad        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_sw         (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_lw         (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_mic        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_gwd        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_ph1        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_conv       (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_bl         (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_stph       (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_cld        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_iau        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_nud        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_tot        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_adv        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_sol        (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dPV_mass       (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   adv_only_PV    (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_0       (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_bl      (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_bl_mix  (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_bl_LH   (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_conv    (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_mic     (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_rad     (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_SW      (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_LW      (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &
                    
,   dtheta_slow    (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)           &

,   dtheta_cld     (tdims%i_start:tdims%i_end,           &
                    tdims%j_start:tdims%j_end,           &
                    tdims%k_start:tdims%k_end)


! Non-zero on exit if error detected.
INTEGER, INTENT(INOUT) :: errorstatus

CHARACTER(LEN=*), PARAMETER  :: routinename='EG_SL_DYN_TR'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_pv_tracer) CALL eg_sl_PV_trac(                                  &
             row_length, rows, n_rows, model_levels,                  &
             g_i_pe, high_order_scheme, monotone_scheme,              &
             l_high, l_mono,                                          &
             dPV_rad, dPV_sw, dPV_lw,                                 &
             dPV_mic, dPV_gwd, dPV_ph1,                               &
             dPV_conv, dPV_bl, dPV_stph, dPV_cld,                     &
             dPV_iau, dPV_nud, dPV_tot,                               &
             dPV_adv, dPV_sol, dPV_mass, adv_only_PV,                 &
             errorstatus)

IF (l_diab_tracer)   CALL eg_sl_diab_trac(                            &
             row_length, rows, n_rows, model_levels,                  &
             g_i_pe, high_order_scheme, monotone_scheme,              &
             l_high, l_mono,                                          &
             dtheta_0, dtheta_bl, dtheta_bl_mix, dtheta_bl_LH,        &
             dtheta_conv, dtheta_mic, dtheta_rad,                     &
             dtheta_SW, dtheta_LW, dtheta_slow, dtheta_cld,           &
             errorstatus)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE eg_sl_dyn_tr

!-----------------------------------------------------------------
SUBROUTINE dyn_tr_fast( r_u, r_v, theta,                                &
                        exner_theta_levels,                             &
                        ! PV-tracers
                        dPV_adv, dPV_conv, dPV_bl, dPV_stph,            &
                        ! theta tracers
                        dtheta_bl, dtheta_bl_mix, dtheta_bl_LH,         &
                        dtheta_conv)

USE physics_tendencies_mod,  ONLY: init_convection_tendencies
USE cv_run_mod,              ONLY: l_param_conv

IMPLICIT NONE
REAL, INTENT (IN) ::                                                    &
! main prog. U field after advection
  r_u         (udims_s%i_start:udims_s%i_end,                           &
               udims_s%j_start:udims_s%j_end,                           &
               udims_s%k_start:udims_s%k_end)                           &
! main prog. V field after advection
, r_v         (vdims_s%i_start:vdims_s%i_end,                           &
               vdims_s%j_start:vdims_s%j_end,                           &
               vdims_s%k_start:vdims_s%k_end)                           &
! main prog. theta field after advection
, theta       (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! Exner pressure on theta levels (to convert T -> theta)
, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::                                                  &
! PV Increments from advection
  dPV_adv        (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)                            &
! PV Increments from convection parametrization
, dPV_conv       (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)                            &
! PV Increments from Boundary Layer parametrization
, dPV_bl         (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)                            &
! PV Increments from stochastic-physics
, dPV_stph       (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)                            &
! diab. Increments from BL parametrization
, dtheta_bl      (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)                            &
! diab. Increments from BL mixing parametrization
, dtheta_bl_mix  (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)                            &
! diab. Increments from LH BL parametrization
, dtheta_bl_LH   (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)                            &
! diab. Increments from convection parametrization
, dtheta_conv    (tdims%i_start:tdims%i_end,                            &
                  tdims%j_start:tdims%j_end,                            &
                  tdims%k_start:tdims%k_end)

CHARACTER(LEN=*), PARAMETER  :: routinename='DYN_TR_FAST'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! If convection is not switched on, initialise tendency arrays
IF (.NOT. l_param_conv) THEN
  
  CALL init_convection_tendencies
  
END IF

! PV TRACER:  Compute fast + stph dPV increments
IF (l_pv_tracer)  CALL PV_tracer_fast( r_u, r_v, theta,                 &
                                       exner_theta_levels, dPV_adv,     &
                                       dPV_conv, dPV_bl, dPV_stph)

! Theta tracer:  Compute fast + stph dPV increments
IF (l_diab_tracer ) CALL diab_tracer_fast(exner_theta_levels,           &
                                          dtheta_bl, dtheta_bl_mix,     &
                                          dtheta_bl_LH, dtheta_conv)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dyn_tr_fast

!-----------------------------------------------------------------
SUBROUTINE dyn_tr_np1(rho, exner_theta_levels,                          &
                      dPV_sol, dPV_mass, dPV_cld,dtheta_cld)

USE cloud_inputs_mod,        ONLY: i_cld_vn
USE pc2_constants_mod,       ONLY: i_cld_pc2

IMPLICIT NONE

REAL, INTENT (IN) ::                                                    &
  rho        (pdims_s%i_start:pdims_s%i_end,                            &
              pdims_s%j_start:pdims_s%j_end,                            &
              pdims_s%k_start:pdims_s%k_end)                            &
           ! main prog. Rho field

, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)
           ! Exner pressure on theta levels (to convert T -> theta)

REAL, INTENT(INOUT) ::                                                  &
  dPV_sol     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
           ! PV Increments from Helmholtz solver
, dPV_mass    (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
           ! PV Increments from mass update
, dPV_cld     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
           ! PV Increments from cloud rebalancing
, dtheta_cld  (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)
           ! theta increments from cloud rebalancing
CHARACTER(LEN=*), PARAMETER  :: routinename='DYN_TR_NP1'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_pv_tracer) CALL PV_tracer_np1( rho, exner_theta_levels,           &
                                     dPV_sol, dPV_mass, dPV_cld)



IF (l_diab_tracer .AND. (i_cld_vn == i_cld_pc2))                        &
   CALL diab_tracer_np1(exner_theta_levels, dtheta_cld)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dyn_tr_np1

!-----------------------------------------------------------------
SUBROUTINE dyn_tr_end ( rho, exner_theta_levels,                        &
                        dPV_iau, dPV_nud, dPV_tot)

IMPLICIT NONE

REAL, INTENT (IN) ::                                                    &
! main prog. Rho field
   rho         (pdims_s%i_start:pdims_s%i_end,                          &
                pdims_s%j_start:pdims_s%j_end,                          &
                pdims_s%k_start:pdims_s%k_end)                          &

! Exner pressure on theta levels (to convert T -> theta)
, exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                     &
                     tdims_s%j_start:tdims_s%j_end,                     &
                     tdims_s%k_start:tdims_s%k_end)

REAL, INTENT(INOUT) ::                                                  &
! PV Increments from IAU
  dPV_iau     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from nudging
, dPV_nud     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments for all timestep
, dPV_tot     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)

CHARACTER(LEN=*), PARAMETER  :: routinename='DYN_TR_END'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
!End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (l_pv_tracer) CALL PV_tracer_end( rho, exner_theta_levels,           &
                                     dPV_iau, dPV_nud, dPV_tot)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE dyn_tr_end

END MODULE dyn_tracers_mod
