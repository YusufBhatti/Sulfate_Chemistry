 ! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          Compute PV increments from slow physics parametrizations
!          (radiation, microphysics and GWD) and add them in to the
!          dPV. an additional array holds all inc from slow physics
!          (which holds energy correction plus pc2_turbulence)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: Computes the diabatic component of the radiation
!                   and microphysics, and the diabatic and friction of
!                   GWD. Also dPV_ph1 is computed and PV_0 and PV_1 for
!                   as the PV value at the beginning of the timestep and
!                   right before advection respectively

MODULE PV_tracer_slow_mod

IMPLICIT NONE
CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'PV_TRACER_SLOW_MOD'
CONTAINS

SUBROUTINE PV_tracer_slow( rho, exner_theta_levels,                     &
                           dPV_rad, dPV_sw, dPV_lw,                     &
                           dPV_mic, dPV_gwd, dPV_ph1,                   &
                           adv_only_PV)


USE halo_exchange,           ONLY: swap_bounds
USE mpp_conf_mod,            ONLY: swap_field_is_vector, swap_field_is_scalar

USE atm_fields_bounds_mod,   ONLY: tdims, pdims, udims , vdims,         &
                                   tdims_s, pdims_s, udims_s, vdims_s

USE free_tracers_inputs_mod, ONLY: l_pv_dyn, l_pv_split_rad,            &
                                   l_pv_adv_only, l_calc_pv_full


USE level_heights_mod,       ONLY: r_rho_levels

USE Field_Types,             ONLY: fld_type_u, fld_type_v, fld_type_p

USE atm_step_local,          ONLY: first_atmstep_call

USE atm_fields_mod,          ONLY: u, v, theta

USE physics_tendencies_mod,  ONLY: dt_lw, dt_sw, dt_gwd, du_gwd,        &
                                   dv_gwd, dt_mic, dtheta_ph1

USE PV_tracer_sav_4A_mod,    ONLY: rho_0, pv_0, pv_1, si_to_pvunits

USE calc_PV_full_mod,        ONLY: calc_pv_new_or_old

! DrHook parameters
USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim

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
! PV Increments from SW radiation parametrization
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
               tdims%k_start:tdims%k_end)

! Internal variables
INTEGER :: i,j,k
           ! intergers for do loops

LOGICAL :: l_coriolis
! Variable to set coriolis terms to zero (F) in PV computation

REAL ::                                                                 &
! Updated u to compute PV
  unow        (udims_s%i_start:udims_s%i_end,                           &
               udims_s%j_start:udims_s%j_end,                           &
               udims_s%k_start:udims_s%k_end)                           &
! Updated v to compute PV
, vnow        (vdims_s%i_start:vdims_s%i_end,                           &
               vdims_s%j_start:vdims_s%j_end,                           &
               vdims_s%k_start:vdims_s%k_end)                           &
! Updated theta to compute PV
, thetanow    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)                           &
! Temporal field of PV for physics
, pv_temp     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_end)

CHARACTER(LEN=*), PARAMETER  :: routinename='PV_TRACER_SLOW'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ++++++++++++++++++++++++++++++++++++
! Initialization of internal variables
! ++++++++++++++++++++++++++++++++++++
IF (l_calc_pv_full) THEN
  DO k = pdims_s%k_start, pdims_s%k_end
    DO j = pdims_s%j_start, pdims_s%j_end
      DO i = pdims_s%i_start, pdims_s%i_end
        rho_0(i,j,k) = rho(i,j,k)
      END DO
    END DO
  END DO
ELSE
  ! Note: The default PV computation routine Calc_PV_at_theta
  ! divides rho by a r square factor (old legacy from wet rho).
  ! Hence it needs to be multiplied by r square.
  DO k = pdims%k_start, pdims%k_end
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        rho_0(i,j,k) = rho(i,j,k)*r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
      END DO
    END DO
  END DO

  ! Swap bounds
  CALL swap_bounds(rho_0,                                               &
                   pdims%i_len, pdims%j_len, pdims%k_len,               &
                   pdims_s%halo_i, pdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)

END IF ! end if over l_calc_pv_full

! ++++++++++++++++++++++++++++++++++++
! Compute initial PV_0 for dPV_tot
! ++++++++++++++++++++++++++++++++++++

! Compute theta
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetanow(i,j,k) = theta(i,j,k)
    END DO
  END DO
END DO
! Set level zero equal to level 1
thetanow(:,:,0) = thetanow(:,:,1)

CALL swap_bounds(thetanow,                                              &
                 tdims%i_len, tdims%j_len, tdims%k_len,                 &
                 tdims_s%halo_i, tdims_s%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)


l_coriolis=.TRUE.
CALL calc_pv_new_or_old (l_coriolis, u, v, thetanow, rho_0, pv_temp)

! Compute pv_0
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      pv_0(i,j,k) = pv_temp(i,j,k)
    END DO
  END DO
END DO


! Initialise advection only PV tracer to initial PV
IF (first_atmstep_call .AND. l_pv_adv_only) THEN
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        adv_only_PV(i,j,k) = pv_0(i,j,k)*si_to_pvunits
      END DO
    END DO
  END DO
END IF


! ++++++++++++++++++++++++++++++++++++
! Compute radiation dPVs
! ++++++++++++++++++++++++++++++++++++

IF (l_pv_split_rad) THEN
! ++ Compute SW first
  ! Update theta
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = dt_sw(i,j,k) / exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  ! Do swap_bounds for 'new' variables
  CALL swap_bounds(thetanow,                                              &
                   tdims%i_len, tdims%j_len, tdims%k_len,                 &
                   tdims_s%halo_i, tdims_s%halo_j,                        &
                   fld_type_p,swap_field_is_scalar)
  
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, u, v, thetanow, rho_0, pv_temp)

  ! Update dPV_sw (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_sw(i,j,k) = dPV_sw(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

! ++ Then Compute LW
  ! Update theta
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = dt_lw(i,j,k) / exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  ! Do swap_bounds for 'new' variables
  CALL swap_bounds(thetanow,                                              &
                   tdims%i_len, tdims%j_len, tdims%k_len,                 &
                   tdims_s%halo_i, tdims_s%halo_j,                        &
                   fld_type_p,swap_field_is_scalar)

  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, u, v, thetanow, rho_0, pv_temp)

  ! Update dPV_lw (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_lw(i,j,k) = dPV_lw(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

ELSE
!++ Compute SW and LW together

  ! Update theta
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = ( dt_sw(i,j,k) + dt_lw(i,j,k) )                 &
                          / exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  ! Do swap_bounds for 'new' variables
  CALL swap_bounds(thetanow,                                              &
                   tdims%i_len, tdims%j_len, tdims%k_len,                 &
                   tdims_s%halo_i, tdims_s%halo_j,                        &
                   fld_type_p,swap_field_is_scalar)

  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, u, v, thetanow, rho_0, pv_temp)

  ! Update dPV_rad (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_rad(i,j,k) = dPV_rad(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

END IF ! end if over l_pv_split_rad

! ++++++++++++++++++++++++++++++++++++
! Compute Microphysics dPV
! ++++++++++++++++++++++++++++++++++++

! Update theta
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetanow(i,j,k) = dt_mic(i,j,k) / exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
! Set level zero equal to level 1
thetanow(:,:,0) = thetanow(:,:,1)

! Do swap_bounds for 'new' variables
CALL swap_bounds(thetanow,                                              &
                 tdims%i_len, tdims%j_len, tdims%k_len,                 &
                 tdims_s%halo_i, tdims_s%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)

l_coriolis=.TRUE.
CALL calc_pv_new_or_old (l_coriolis, u, v, thetanow, rho_0, pv_temp)

! Update dPV_mic (convert to PVU = 1e(-6) m**2 K / (Kg s) )
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_mic(i,j,k) = dPV_mic(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
    END DO
  END DO
END DO

! ++++++++++++++++++++++++++++++++++++
! Compute GWD dPV
! ++++++++++++++++++++++++++++++++++++

! ---  Compute diabatic term first --

! Update theta
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetanow(i,j,k) = dt_gwd(i,j,k) / exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
! Set level zero equal to level 1
thetanow(:,:,0) = thetanow(:,:,1)

! Do swap_bounds for 'new' variables
CALL swap_bounds(thetanow,                                              &
                 tdims%i_len, tdims%j_len, tdims%k_len,                 &
                 tdims_s%halo_i, tdims_s%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)

l_coriolis=.TRUE.
CALL calc_pv_new_or_old (l_coriolis, u, v, thetanow, rho_0, pv_temp)

! Update dPV_gwd (convert to PVU = 1e(-6) m**2 K / (Kg s) )
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_gwd(i,j,k) = dPV_gwd(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
    END DO
  END DO
END DO

! ---  Add in the friction term to dPV_gwd --

DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetanow(i,j,k) = theta(i,j,k)
    END DO
  END DO
END DO
! Set level zero equal to level 1
thetanow(:,:,0) = thetanow(:,:,1)


DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      unow(i,j,k) = du_gwd(i,j,k)
    END DO
  END DO
END DO

DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      vnow(i,j,k) = dv_gwd(i,j,k)
    END DO
  END DO
END DO

! Do swap_bounds for 'new' variables
CALL swap_bounds(thetanow,                                              &
                 tdims%i_len, tdims%j_len, tdims%k_len,                 &
                 tdims_s%halo_i, tdims_s%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)

CALL swap_bounds(unow,                                                  &
                 udims%i_len, udims%j_len, udims%k_len,                 &
                 udims_s%halo_i, udims_s%halo_j,                        &
                 fld_type_u,swap_field_is_vector)

CALL swap_bounds(vnow,                                                  &
                 vdims%i_len, vdims%j_len, vdims%k_len,                 &
                 vdims_s%halo_i, vdims_s%halo_j,                        &
                 fld_type_v,swap_field_is_vector)

l_coriolis=.FALSE.
CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

! Add friction incr to dPV_gwd
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_gwd(i,j,k) = dPV_gwd(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
    END DO
  END DO
END DO

! ++++++++++++++++++++++++++++++++++++
! Compute All slow physics dPV
! ++++++++++++++++++++++++++++++++++++

! ---  Add in the friction term --
! The friction sources are only GWD, so we can take the previous pv_temp

! Update dPV_ph1 with GWD friction pv increment
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_ph1(i,j,k) = dPV_ph1(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
    END DO
  END DO
END DO

! ---  Compute diabatic term for dPV_ph1 --

! Update theta
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetanow(i,j,k) = dtheta_ph1(i,j,k)
    END DO
  END DO
END DO
! Set level zero equal to level 1
thetanow(:,:,0) = thetanow(:,:,1)

! Do swap_bounds for 'new' variables
CALL swap_bounds(thetanow,                                              &
                 tdims%i_len, tdims%j_len, tdims%k_len,                 &
                 tdims_s%halo_i, tdims_s%halo_j,                        &
                 fld_type_p,swap_field_is_scalar)

l_coriolis=.TRUE.
CALL calc_pv_new_or_old (l_coriolis, u, v, thetanow, rho_0, pv_temp)

! Update dPV_ph1 with diabatic sources
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_ph1(i,j,k) = dPV_ph1(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
    END DO
  END DO
END DO

! ++++++++++++++++++++++++++++++++++++
! Compute PV_1 if requested
! ++++++++++++++++++++++++++++++++++++

! Add it in PV_1 to compute dPV_adv if requested
IF (l_pv_dyn) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = theta(i,j,k) + dtheta_ph1(i,j,k)
      END DO
    END DO
  END DO
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = u(i,j,k) + du_gwd(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = v(i,j,k) + dv_gwd(i,j,k)
      END DO
    END DO
  END DO

  ! Do swap_bounds for 'new' variables
  CALL swap_bounds(thetanow,                                            &
                   tdims%i_len, tdims%j_len, tdims%k_len,               &
                   tdims_s%halo_i, tdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)

  CALL swap_bounds(unow,                                                &
                   udims%i_len, udims%j_len, udims%k_len,               &
                   udims_s%halo_i, udims_s%halo_j,                      &
                   fld_type_u,swap_field_is_vector)

  CALL swap_bounds(vnow,                                                &
                   vdims%i_len, vdims%j_len, vdims%k_len,               &
                   vdims_s%halo_i, vdims_s%halo_j,                      &
                   fld_type_v,swap_field_is_vector)

 l_coriolis=.TRUE.
 CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Compute PV_1
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        pv_1(i,j,k) = pv_temp(i,j,k)
      END DO
    END DO
  END DO

END IF !end if over l_pv_dyn

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE PV_tracer_slow

END MODULE PV_tracer_slow_mod
