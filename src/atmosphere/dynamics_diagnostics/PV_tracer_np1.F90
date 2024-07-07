! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose:
!          Compute PV increments from cloud rebalancing and,
!          if requested, dynamical dPV (solver + density update)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics diagnostics
!
! Code Description: Computes the diabatic term for cloud and if requested
!                   the differences in PV before and after the
!                   dynamical processes of sover and density update
!                   convection.
!

MODULE PV_tracer_np1_mod

IMPLICIT NONE
CHARACTER(LEN=*), PRIVATE, PARAMETER :: ModuleName = 'PV_TRACER_NP1_MOD'
CONTAINS

SUBROUTINE PV_tracer_np1( rho, exner_theta_levels,                      &
                          dPV_sol, dPV_mass, dPV_cld)

USE halo_exchange,           ONLY: swap_bounds
USE mpp_conf_mod,            ONLY: swap_field_is_vector, swap_field_is_scalar

USE atm_fields_bounds_mod,   ONLY: tdims, pdims, udims , vdims,         &
                                   tdims_s, pdims_s, udims_s, vdims_s

USE free_tracers_inputs_mod, ONLY: l_pv_dyn, l_calc_pv_full

USE Field_Types,             ONLY: fld_type_u, fld_type_v, fld_type_p

USE calc_PV_full_mod,        ONLY: calc_pv_new_or_old

USE IAU_mod,                 ONLY: l_iau

USE atm_fields_mod,          ONLY: u, v, theta

USE level_heights_mod,       ONLY: r_rho_levels

USE PV_tracer_sav_4A_mod,    ONLY: rho_0, pv_1, si_to_pvunits

USE physics_tendencies_mod,  ONLY: dt_cld

USE cloud_inputs_mod,        ONLY: i_cld_vn

USE pc2_constants_mod,       ONLY: i_cld_pc2

! DrHook parameters
USE yomhook,                 ONLY: lhook, dr_hook
USE parkind1,                ONLY: jprb, jpim

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
               tdims%k_start:tdims%k_end)
           ! PV Increments from cloud rebalancing

! Internal variables
INTEGER :: i,j,k
           ! intergers for do loops

LOGICAL :: l_coriolis
! Variable to set coriolis terms to zero (F) in PV computation

REAL ::                                                                 &
  unow        (udims_s%i_start:udims_s%i_end,                           &
               udims_s%j_start:udims_s%j_end,                           &
               udims_s%k_start:udims_s%k_end)                           &
           ! Updated u to compute PV
, vnow        (vdims_s%i_start:vdims_s%i_end,                           &
               vdims_s%j_start:vdims_s%j_end,                           &
               vdims_s%k_start:vdims_s%k_end)                           &
           ! Updated v to compute PV
, thetanow    (tdims_s%i_start:tdims_s%i_end,                           &
               tdims_s%j_start:tdims_s%j_end,                           &
               tdims_s%k_start:tdims_s%k_end)                           &
           ! Updated theta to compute PV
, pv_temp     (tdims%i_start:tdims%i_end,                               &
               tdims%j_start:tdims%j_end,                               &
               tdims%k_end)
           ! Temporal value of PV

REAL, ALLOCATABLE ::                                                    &
  pv_np1      (:,:,:)
           ! PV value at n+1 if requested by l_pv_dyn

CHARACTER(LEN=*), PARAMETER  :: RoutineName='PV_TRACER_NP1'

!DrHook variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
! End of header

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ++++++++++++++++++++++++++++++++++++
! Initialization of internal variables
! ++++++++++++++++++++++++++++++++++++

! Update winds
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

! Do swap_bounds for 'new' variables
CALL swap_bounds(unow,                                                  &
                 udims%i_len, udims%j_len, udims%k_len,                 &
                 udims_s%halo_i, udims_s%halo_j,                        &
                 fld_type_u,swap_field_is_vector)

CALL swap_bounds(vnow,                                                  &
                 vdims%i_len, vdims%j_len, vdims%k_len,                 &
                 vdims_s%halo_i, vdims_s%halo_j,                        &
                 fld_type_v,swap_field_is_vector)

! ++++++++++++++++++++++++++++++++++++++
! Compute dPV Solver + dPV mass update
! ++++++++++++++++++++++++++++++++++++++

IF (l_pv_dyn) THEN

  IF (i_cld_vn == i_cld_pc2) THEN
    DO k = 1, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          thetanow(i,j,k) = theta(i,j,k) -                              &
                            dt_cld(i,j,k) / exner_theta_levels(i,j,k)
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
  END IF
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  CALL swap_bounds(thetanow,                                            &
                   tdims%i_len, tdims%j_len, tdims%k_len,               &
                   tdims_s%halo_i, tdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)

  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

! allocace pv_np1
  IF (.NOT. ALLOCATED(pv_np1))                                          &
    ALLOCATE(pv_np1(tdims%i_start:tdims%i_end,                          &
                    tdims%j_start:tdims%j_end,                          &
                    tdims%k_start:tdims%k_end))
           ! Temporal value of PV

 ! Copy np1
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        pv_np1(i,j,k) = pv_temp(i,j,k)
      END DO
    END DO
  END DO



 ! Update dPV_sol (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_sol(i,j,k) = dPV_sol(i,j,k) +                               &
                        (pv_np1(i,j,k) - pv_1(i,j,k) ) *si_to_pvunits
      END DO
    END DO
  END DO

END IF !end if l_pv_dyn

! Update rho
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
        rho_0(i,j,k) = rho(i,j,k)*r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
      END DO
    END DO
  END DO

  ! Swap bounds
  CALL swap_bounds(rho_0,                                               &
                   pdims%i_len, pdims%j_len, pdims%k_len,               &
                   pdims_s%halo_i, pdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)

END IF

IF (l_pv_dyn) THEN
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

 ! Update dPV_mass (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_mass(i,j,k) = dPV_mass(i,j,k) +                             &
                          (pv_temp(i,j,k) - pv_np1(i,j,k) ) *si_to_pvunits
      END DO
    END DO
  END DO

  IF (ALLOCATED(pv_np1)) DEALLOCATE(pv_np1)

END IF
! ++++++++++++++++++++++++++++++++++++
! Compute Cloud dPV
! ++++++++++++++++++++++++++++++++++++

IF (i_cld_vn == i_cld_pc2) THEN

  ! Update theta
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        thetanow(i,j,k) = dt_cld(i,j,k) /exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO
  ! Set level zero equal to level 1
  thetanow(:,:,0) = thetanow(:,:,1)

  ! Do swap_bounds for 'new' variables
  CALL swap_bounds(thetanow,                                            &
                   tdims%i_len, tdims%j_len, tdims%k_len,               &
                   tdims_s%halo_i, tdims_s%halo_j,                      &
                   fld_type_p,swap_field_is_scalar)
  
  l_coriolis=.TRUE.
  CALL calc_pv_new_or_old (l_coriolis, unow, vnow, thetanow, rho_0, pv_temp)

  ! Update dPV_cld (convert to PVU = 1e(-6) m**2 K / (Kg s) )
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        dPV_cld(i,j,k) = dPV_cld(i,j,k) + pv_temp(i,j,k) *si_to_pvunits
      END DO
    END DO
  END DO

END IF ! end if over i_cld_vn == i_cld_pc2

! ++++++++++++++++++++++++++++++++++++
! Initialize IAU PV
! ++++++++++++++++++++++++++++++++++++
! Get PV before theta
IF (l_iau) THEN
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

 ! Update pv_1
  DO k = 1, tdims%k_end
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        pv_1(i,j,k) = pv_temp(i,j,k)
      END DO
    END DO
  END DO

END IF ! end if IAU

! End of routine
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE PV_tracer_np1

END MODULE PV_tracer_np1_mod
