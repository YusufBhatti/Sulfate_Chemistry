! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          Compute PV increments from fast physics parametrizations
!          (convection and boundary layer) and add them in to the
!          respective dPV arrays. An extra option can compute the
!          dPV of stochastic physics schemes if any of them is
!          active.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: Computes the diabatic and friction term for
!                   convection. Updates fields and computes the dPV
!                   for BL. If requested, fields are updated to compute
!                   the stochastic physics dPV and store values for the
!                   dynamical solver. It does compute the advection dPV
!                   If requested
!

MODULE PV_tracer_fast_mod

IMPLICIT NONE
CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'PV_TRACER_FAST_MOD'
CONTAINS

SUBROUTINE PV_tracer_fast(r_u, r_v, theta, exner_theta_levels,          &
                          dPV_adv, dPV_conv, dPV_bl, dPV_stph)

USE halo_exchange,              ONLY: swap_bounds
USE mpp_conf_mod,               ONLY: swap_field_is_vector,             &
                                      swap_field_is_scalar

USE atm_fields_bounds_mod,      ONLY: tdims, pdims, udims, vdims,       &
                                      tdims_s, pdims_s, udims_s,        &
                                      vdims_s

USE free_tracers_inputs_mod,    ONLY: l_pv_dyn

USE Field_Types,                ONLY: fld_type_u, fld_type_v,           &
                                      fld_type_p

USE cv_run_mod,                 ONLY: l_param_conv

USE stochastic_physics_run_mod, ONLY: l_spt, l_skeb2

USE atm_fields_mod,             ONLY: u,v

USE physics_tendencies_mod,     ONLY: dt_conv, du_conv, dv_conv,        &
                                      dt_bl, du_bl, dv_bl,              &
                                      dtheta_stph, du_stph, dv_stph

USE PV_tracer_sav_4A_mod,       ONLY: rho_0, pv_1, si_to_pvunits

USE calc_PV_full_mod,           ONLY: calc_pv_new_or_old

! DrHook parameters
USE yomhook,                    ONLY: lhook, dr_hook
USE parkind1,                   ONLY: jprb, jpim

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
  dPV_adv     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from convection parametrization
, dPV_conv    (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from Boundary Layer parametrization
, dPV_bl      (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)                               &
! PV Increments from stochastic-physics
, dPV_stph    (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_start:tdims%k_end)

! Internal variables
! intergers for do loops
INTEGER :: i,j,k

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
! Temporal value of PV
, pv_temp     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_end)

CHARACTER(LEN=*), PARAMETER  :: routinename='PV_TRACER_FAST'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ++++++++++++++++++++++++++++++++++++
! Compute Initial variables
! ++++++++++++++++++++++++++++++++++++

DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      unow(i,j,k) = u(i,j,k) + r_u(i,j,k)
    END DO
  END DO
END DO

DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      vnow(i,j,k) = v(i,j,k) + r_v(i,j,k)
    END DO
  END DO
END DO

CALL swap_bounds(unow,                                                  &
                 udims%i_len, udims%j_len, udims%k_len,                 &
                 udims_s%halo_i, udims_s%halo_j,                        &
                 fld_type_u,swap_field_is_vector)

CALL swap_bounds(vnow,                                                  &
                 vdims%i_len, vdims%j_len, vdims%k_len,                 &
                 vdims_s%halo_i, vdims_s%halo_j,                        &
                 fld_type_v,swap_field_is_vector)

! ++++++++++++++++++++++++++++++++++++
! Compute Advection dPV
! ++++++++++++++++++++++++++++++++++++

IF (l_pv_dyn) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = theta(i,j,k)
      END DO
    END DO
  END DO
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  CALL swap_bounds(thetanow,                                            &
                   tdims%i_len, tdims%j_len, tdims%k_len,               &
                   tdims_s%halo_i, tdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)
  
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Update dPV_adv
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_adv(i,j,k) = dPV_adv(i,j,k) +                               &
                         ( pv_temp(i,j,k) - pv_1(i,j,k) ) *si_to_pvunits
      END DO
    END DO
  END DO

END IF ! end if l_pv_dyn

! ++++++++++++++++++++++++++++++++++++
! Compute Convection dPV
! ++++++++++++++++++++++++++++++++++++

! If convection is not switched on, create arrays
IF (l_param_conv) THEN

  ! ---  Compute diabatic term first --

  ! Update theta
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = dt_conv(i,j,k) /exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  CALL swap_bounds(thetanow,                                            &
                   tdims%i_len, tdims%j_len, tdims%k_len,               &
                   tdims_s%halo_i, tdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)

  ! U and V  are already computed in the initialization
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Update dPV_conv
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_conv(i,j,k) = dPV_conv(i,j,k) + pv_temp(i,j,k)*si_to_pvunits
      END DO
    END DO
  END DO

  ! ---  Add in the friction term --

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
        unow(i,j,k) = du_conv(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = dv_conv(i,j,k)
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

  l_coriolis=.FALSE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Add friction source term of PV to dPV array
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_conv(i,j,k) = dPV_conv(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

END IF !end if over l_param_conv

! ++++++++++++++++++++++++++++++++++++
! Compute Boundary Layer dPV
! ++++++++++++++++++++++++++++++++++++

! Update theta
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetanow(i,j,k) = dt_bl(i,j,k) /exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
! Set level zero equal to level 1
thetanow(:,:,0) = thetanow(:,:,1)

! Update values for horizontal wind
DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      unow(i,j,k) = u(i,j,k) + r_u(i,j,k) + du_conv(i,j,k)
    END DO
  END DO
END DO

DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      vnow(i,j,k) = v(i,j,k) + r_v(i,j,k) + dv_conv(i,j,k)
   END DO
  END DO
END DO

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

l_coriolis=.TRUE.
CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

! Add PV inc to dPV_bl
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_bl(i,j,k) = dPV_bl(i,j,k) +  pv_temp(i,j,k)*si_to_pvunits
    END DO
  END DO
END DO

! ---  Add in the friction term --
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      thetanow(i,j,k) = theta(i,j,k) +                                &
                        dt_conv(i,j,k) /exner_theta_levels(i,j,k)
    END DO
  END DO
END DO
! Set level zero equal to level 1
thetanow(:,:,0) = thetanow(:,:,1)

DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      unow(i,j,k) = du_bl(i,j,k)
    END DO
  END DO
END DO

DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      vnow(i,j,k) = dv_bl(i,j,k)
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

! Add friction source term of PV to dPV array
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_bl(i,j,k) = dPV_bl(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
    END DO
  END DO
END DO

! ++++++++++++++++++++++++++++++++++++
! Compute Stochastic physics dPV_stph
! ++++++++++++++++++++++++++++++++++++

! ---  Add in the friction term --
! Only if L_SPT or L_SKEB2

IF (l_spt .OR. l_skeb2) THEN

  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = theta(i,j,k) + (dt_conv(i,j,k) + dt_bl(i,j,k))  &
                                         /exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = du_stph(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = dv_stph(i,j,k)
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
  
  l_coriolis=.FALSE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Update dPV_conv
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_stph(i,j,k) = dPV_stph(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

END IF  !end if over l_skeb2 or l_spt

! ---  Compute diabatic term (only if SPT is active) --

IF (l_spt) THEN
  ! Update theta
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = dtheta_stph(i,j,k)
      END DO
    END DO
  END DO
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  ! Update values for horizontal wind
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = u(i,j,k) + r_u(i,j,k) +                           &
                      du_conv(i,j,k) + du_bl(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) =  v(i,j,k) + r_v(i,j,k) +                          &
                       dv_conv(i,j,k) + dv_bl(i,j,k)
      END DO
    END DO
  END DO

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

  ! Update dPV_stph
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_stph(i,j,k) = dPV_stph(i,j,k) +  pv_temp(i,j,k)*si_to_pvunits
      END DO
    END DO
  END DO

END IF !end if over l_spt


! ++++++++++++++++++++++++++++++++++++
! Compute PV prior to the solver
! ++++++++++++++++++++++++++++++++++++
IF (l_pv_dyn) THEN

  ! Update unow variables with BL
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = theta(i,j,k) + (dt_conv(i,j,k) +              &
                                          dt_bl(i,j,k))                 &
                                         /exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = u(i,j,k) + r_u(i,j,k) +                           &
                      du_conv(i,j,k) + du_bl(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = v(i,j,k) + r_v(i,j,k) +                           &
                      dv_conv(i,j,k) + dv_bl(i,j,k)
      END DO
    END DO
  END DO

  IF (l_spt .OR. l_skeb2) THEN
    ! compute thetanow
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          thetanow(i,j,k) = thetanow(i,j,k) + dtheta_stph(i,j,k)
        END DO
      END DO
    END DO

    DO k = udims%k_start, udims%k_end
      DO j = udims%j_start, udims%j_end
        DO i = udims%i_start, udims%i_end
          unow(i,j,k) = unow(i,j,k) + du_stph(i,j,k)
        END DO
      END DO
    END DO

    DO k = vdims%k_start, vdims%k_end
      DO j = vdims%j_start, vdims%j_end
        DO i = vdims%i_start, vdims%i_end
          vnow(i,j,k) = vnow(i,j,k) + dv_stph(i,j,k)
        END DO
      END DO
    END DO

  END IF !end if over stoch physics

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

  ! Overwrite PV_1 with PV value prior to solver (and diff)
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Update pv_1
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        pv_1(i,j,k) = pv_temp(i,j,k)
      END DO
    END DO
  END DO

END IF ! end if l_pv_dyn

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE PV_tracer_fast

END MODULE PV_tracer_fast_mod
