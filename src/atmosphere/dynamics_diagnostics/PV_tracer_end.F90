! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          Compute PV increments for all timestep (dPV_tot),
!          if requested it also computes nudging and iau dPV
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: Computes the diabatic and friction terms for
!                   nudging and the PV difference for IAU.
!                   Computes final PV to obtain the dPV_tot.
!

MODULE PV_tracer_end_mod

IMPLICIT NONE
CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'PV_TRACER_END_MOD'
CONTAINS

SUBROUTINE PV_tracer_end( rho, exner_theta_levels,                      &
                          dPV_iau, dPV_nud, dPV_tot)


USE halo_exchange,           ONLY: swap_bounds
USE mpp_conf_mod,            ONLY: swap_field_is_vector, swap_field_is_scalar

USE atm_fields_bounds_mod,   ONLY: tdims, pdims, udims, vdims,          &
                                   tdims_s, pdims_s, udims_s, vdims_s

USE level_heights_mod,       ONLY: r_rho_levels

USE Field_Types,             ONLY: fld_type_u, fld_type_v, fld_type_p

USE nudging_input_mod,       ONLY: l_nudging
USE IAU_mod,                 ONLY: l_iau

USE free_tracers_inputs_mod, ONLY: l_calc_pv_full

USE atm_fields_mod,          ONLY: u,v, theta

USE PV_tracer_sav_4A_mod,    ONLY: PV_tracer_dealloc,                  &
                                   rho_0, pv_1, pv_0, si_to_pvunits

USE calc_PV_full_mod,        ONLY: calc_pv_new_or_old

USE physics_tendencies_mod,  ONLY: dtheta_nud, du_nud, dv_nud

! DrHook parameters
USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim

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

CHARACTER(LEN=*), PARAMETER  :: routinename='PV_TRACER_END'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ++++++++++++++++++++++++++++++++++++
! Initialize fields
! ++++++++++++++++++++++++++++++++++++

IF (l_nudging) THEN
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = u(i,j,k) - du_nud(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = v(i,j,k) - dv_nud(i,j,k)
      END DO
    END DO
  END DO
ELSE ! if not l_nudging
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = u(i,j,k)
      END DO
    END DO
  END DO
  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = v(i,j,k)
      END DO
    END DO
  END DO
END IF !end if over nudging

CALL swap_bounds(unow,                                                  &
                 udims%i_len, udims%j_len, udims%k_len,                 &
                 udims_s%halo_i, udims_s%halo_j,                        &
                 fld_type_u,swap_field_is_vector)

CALL swap_bounds(vnow,                                                  &
                 vdims%i_len, vdims%j_len, vdims%k_len,                 &
                 vdims_s%halo_i, vdims_s%halo_j,                        &
                 fld_type_v,swap_field_is_vector)

! ++++++++++++++++++++++++++++++++++++
! Compute IAU dPV
! ++++++++++++++++++++++++++++++++++++

IF (l_iau) THEN
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
    ! divides rho by a a r square factor (old legacy from wet rho).
    ! Hence it needs to be multiplied by r square.
    DO k = pdims%k_start, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          rho_0(i,j,k) = rho(i,j,k)*                                    &
                         r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
        END DO
      END DO
    END DO

    ! Swap bounds
    CALL swap_bounds(rho_0,                                             &
                     pdims%i_len, pdims%j_len, pdims%k_len,             &
                     pdims_s%halo_i, pdims_s%halo_j,                    &
                     fld_type_p,swap_field_is_scalar)

  END IF

  IF (l_nudging) THEN
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          thetanow(i,j,k) = theta(i,j,k) - dtheta_nud(i,j,k)
        END DO
      END DO
    END DO
  ELSE
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          thetanow(i,j,k) = theta(i,j,k)
        END DO
      END DO
    END DO
  END IF ! end if over l_nudging

  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)


  CALL swap_bounds(thetanow,                                            &
                   tdims%i_len, tdims%j_len, tdims%k_len,               &
                   tdims_s%halo_i, tdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)
  
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Update dPV_iau (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO  k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_iau(i,j,k) = dPV_iau(i,j,k) +                               &
                         ( pv_temp(i,j,k) - pv_1(i,j,k) ) *si_to_pvunits
      END DO
    END DO
  END DO

END IF !end if IAU

! ++++++++++++++++++++++++++++++++++++
! Compute Nudging dPV
! ++++++++++++++++++++++++++++++++++++

IF (l_nudging) THEN

  ! ---  Compute diabatic term first --
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = dtheta_nud(i,j,k)
      END DO
    END DO
  END DO
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)


  ! U and v has been computed for dPV_IAU

  CALL swap_bounds(thetanow,                                            &
                   tdims%i_len, tdims%j_len, tdims%k_len,               &
                   tdims_s%halo_i, tdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)
  
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Update dPV_nud (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO  k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_nud(i,j,k) = dPV_nud(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

  ! ---  Add in the friction term --
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = theta(i,j,k) - dtheta_nud(i,j,k)
      END DO
    END DO
  END DO
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = du_nud(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = dv_nud(i,j,k)
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

  ! Add friction source term of PV to dPV array
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_nud(i,j,k) = dPV_nud(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

END IF ! l_nudging

! ++++++++++++++++++++++++++++++++++++
! Compute total dPV
! ++++++++++++++++++++++++++++++++++++

! If nudging -> Get u,v; if not, they've been already computed
! on the initial section
IF (l_nudging) THEN
  DO k = udims%k_start, udims%k_end
    DO j = udims%j_start, udims%j_end
      DO i = udims%i_start, udims%i_end
        unow(i,j,k) = u(i,j,k)
      END DO
    END DO
  END DO

  DO k = vdims%k_start, vdims%k_end
    DO j = vdims%j_start, vdims%j_end
      DO i = vdims%i_start, vdims%i_end
        vnow(i,j,k) = v(i,j,k)
      END DO
    END DO
  END DO

  CALL swap_bounds(unow,                                                &
                   udims%i_len, udims%j_len, udims%k_len,               &
                   udims_s%halo_i, udims_s%halo_j,                      &
                   fld_type_u,swap_field_is_vector)

  CALL swap_bounds(vnow,                                                &
                   vdims%i_len, vdims%j_len, vdims%k_len,               &
                   vdims_s%halo_i, vdims_s%halo_j,                      &
                   fld_type_v,swap_field_is_vector)

END IF !end if over l_nudging

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
CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

! Update dPV_tot (convert to PVU = 1e(-6) m**2 K / (Kg s) )
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dPV_tot(i,j,k) = dPV_tot(i,j,k) +                                 &
                       ( pv_temp(i,j,k) - pv_0(i,j,k) ) * si_to_pvunits
    END DO
  END DO
END DO

! ++++++++++++++++++++++++++++++++++++
! Deallocate arrays
! ++++++++++++++++++++++++++++++++++++

! general variables
CALL PV_tracer_dealloc()

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE PV_tracer_end

END MODULE PV_tracer_end_mod
