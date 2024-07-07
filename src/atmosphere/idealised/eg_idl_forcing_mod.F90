! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_idl_forcing_mod

IMPLICIT NONE
  ! Description: Module to assign, and enact forcing regimes
  ! both temperature and parameterised boundary layer drag
  !  
  !
  ! Method: Temperature forcing uses Newtonian Relaxation of
  ! potential temperature, and the boudnary layer drag is 
  ! simple Rayleigh friction
  !  
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_IDL_FORCING_MOD'

CONTAINS
SUBROUTINE eg_idl_forcing(                                           &
  u, v, theta,  exner_theta_levels, exner, p_theta_levels,           &
  theta_star, q, q_star, qcl_star, qcf_star, qcf2_star, qrain_star,  &
  qgraup_star, rho, r_u, r_v, wetrho_r_sq_n,                         &
  t_surface, problem_number,                                         &
! error information
  error_code  )

USE idl_col_int_diag_mod,      ONLY: idl_col_int_diag
USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod,     ONLY: tdims, tdims_s, pdims, pdims_s, &
                                     udims, udims_s, vdims, vdims_s, &
                                     tdims_l, wdims
USE held_suarez_mod,           ONLY: eg_held_suarez
USE mphys_inputs_mod,          ONLY: l_mcr_qcf2,l_mcr_qrain,l_mcr_qgraup
USE idealise_run_mod,          ONLY: l_HeldSuarez, l_heldsuarez1_drag, &
                                     l_geo_for, p_surface
USE local_heat_mod,            ONLY: local_heat_option, omit

USE geostrophic_forcing_mod,                                       &
                               ONLY: geostrophic_forcing
USE timestep_mod,              ONLY: timestep
USE planet_constants_mod,      ONLY: r, cp, pref, kappa, recip_kappa
USE idealised_diag_mod,        ONLY: dt_inc_ideal_um, dq_inc_ideal_um, &
                                     du_inc_ideal_um, dv_inc_ideal_um, &
                                     dtheta_inc_ideal_um
USE eg_helmholtz_mod,          ONLY: ec_vol
USE eg_idl_theta_forcing_mod,  ONLY: eg_idl_theta_forcing
USE eg_idl_friction_mod,       ONLY: eg_idl_friction
USE tforce_mod,                ONLY: tf_none
USE trelax_mod,                ONLY: tr_none
USE planet_suite_mod,          ONLY: tforce_number, trelax_number
USE bl_option_mod,             ONLY: fric_number, fric_none
USE problem_mod,               ONLY: idealised_planet, dynamical_core,         &
                                     idealised_problem
USE dynamics_testing_mod,      ONLY: l_idealised_data, L_prognostic_level0
USE gen_phys_inputs_mod,       ONLY: l_mr_physics
USE field_types,               ONLY: fld_type_u, fld_type_v, fld_type_w
USE profiles_mod,              ONLY: theta_dry, temperature, u_relax, v_relax, &
                                     theta_relax, theta_inc, mv_inc, w_subs,   &
                                     l_p_profile, p_prof_theta, p_prof_rho
USE p_profile_mod,             ONLY: p_profile
USE profile_relax_mod,         ONLY: profile_relax
USE profile_increment_mod,     ONLY: profile_increment
USE profile_w_force_mod,       ONLY: profile_w_force
USE apply_w_forcing_mod,       ONLY: apply_w_forcing

USE idl_local_heat_mod,        ONLY: idl_local_heat

USE Ereport_mod,               ONLY: Ereport
USE errormessagelength_mod,    ONLY: errormessagelength

IMPLICIT NONE

! Description: Initialises star state for ENDGame when there is no
!              physics. Optionally applies Held-Suarez test of
!              dynamical core (BAMS 75, 1825-1830).


! Data arrays
REAL, INTENT(IN) :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,        &
                                   tdims_s%j_start:tdims_s%j_end,        &
                                   tdims_s%k_start:tdims_s%k_end),       &
  wetrho_r_sq_n(     pdims_s%i_start:pdims_s%i_end,                      &
                     pdims_s%j_start:pdims_s%j_end,                      &
                     pdims_s%k_start:pdims_s%k_end+1)
REAL, INTENT (INOUT) ::                                               &
  theta(             tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end )

REAL, INTENT (INOUT) ::                                               &
  theta_star(        tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  q(                 tdims_l%i_start:tdims_l%i_end,                   &
                     tdims_l%j_start:tdims_l%j_end,                   &
                     tdims_l%k_start:tdims_l%k_end),                  &
  rho(               pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end),                  &
  q_star(            tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  qcl_star(          tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  qcf_star(          tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  qcf2_star(         tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  qrain_star       ( tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  qgraup_star(       tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  r_u(               udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  r_v(               vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end,                   &
                     vdims_s%k_start:vdims_s%k_end),                  &
  u(                 udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  v(                 vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end ,                  &
                     vdims_s%k_start:vdims_s%k_end),                  &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  exner(             pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end+1)

REAL, INTENT(IN) :: t_surface

INTEGER :: problem_number
INTEGER :: error_code
INTEGER :: i, j, k, first_physics_level


! - temporary array holding increment forcing
REAL    :: temp_inc(tdims%i_start:tdims%i_end,                                &
                    tdims%j_start:tdims%j_end,                                &
                    tdims%k_start:tdims%k_end)
REAL    :: w_force(wdims%k_start:wdims%k_end)   ! w profile for forcing (m/s)
                    
LOGICAL :: l_mixing_ratio = .FALSE.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

INTEGER                           :: ErrorStatus
CHARACTER(LEN=errormessagelength) :: Cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='EG_IDL_FORCING'

! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_code = 0

! initial NULL change atmos_physics1 creates tendencies. So if
! nothing happens all tendencies are zero.
theta_star  = 0.0
q_star      = 0.0
qcl_star    = 0.0
qcf_star    = 0.0
IF (l_mcr_qcf2)   qcf2_star   = 0.0
IF (l_mcr_qrain)  qrain_star  = 0.0
IF (l_mcr_qgraup) qgraup_star = 0.0

IF (problem_number == idealised_planet) THEN
   ! Only apply friction if required
   IF (fric_number > fric_none)                                      &
   ! Implement friction or drag (i.e. rayleigh friction)
   CALL eg_idl_friction(                                             &
        ! In data fields
        u, v, exner_theta_levels, exner,                             &
        theta_star, q_star, p_theta_levels,                          &
        ! In/out data fields
        r_u, r_v,                                                    &
        ! Error information
        error_code, l_mixing_ratio)
   ! Again only apply forcing if required
   IF (tforce_number /= tf_none .AND.                                &
        trelax_number /= tr_none)                                    &
   ! Force the potential temperature
   CALL eg_idl_theta_forcing(                                        &
        ! In data fields
        theta, exner_theta_levels, exner,                            &
        ! In/out data fields
        theta_star,                                                  & 
        ! Profile switches
        tforce_number, trelax_number,                                &
        t_surface,                                                   &
        ! Error information
        error_code)   
ELSE IF (problem_number == dynamical_core) THEN
   IF (L_HeldSuarez)                                                 &
   CALL eg_held_suarez(                                              &
        ! in data fields.
        u, v, theta,  exner_theta_levels,                            &
        exner,                                                       &
        ! in/out
        theta_star, r_u, r_v,                                        &
        L_HeldSuarez,                                                &
        L_HeldSuarez1_drag,                                          &
        ! error information
        error_code  )

ELSE IF (problem_number == idealised_problem) THEN

! If any profiles are specified against pressure, then compute
! averaged hydrostatic pressure for current state
  IF (l_p_profile) THEN

    IF (.NOT. ALLOCATED(p_prof_theta)) THEN
      ALLOCATE(p_prof_theta(tdims%k_start:tdims%k_end))
    END IF

    IF (.NOT. ALLOCATED(p_prof_rho)) THEN
      ALLOCATE(p_prof_rho(pdims%k_start:pdims%k_end))
    END IF

    CALL p_profile(p_surface, rho, p_prof_theta, p_prof_rho)

  END IF

  IF (l_geo_for) CALL geostrophic_forcing(r_u, r_v)

! Relaxation -
! - relax fields to their current reference states
  IF (u_relax % n_times > 0) THEN
    CALL profile_relax(u_relax, udims, udims_s, udims_s, fld_type_u, u, r_u)
  END IF
  IF (v_relax % n_times > 0) THEN
    CALL profile_relax(v_relax, vdims, vdims_s, vdims_s, fld_type_v, v, r_v)
  END IF
  IF (theta_relax % n_times > 0) THEN
    CALL profile_relax(theta_relax, tdims, tdims_s, tdims, fld_type_w,       &
                       theta, theta_star)
  END IF

  ! w forcing apply increments to theta and q

  IF (w_subs % n_times > 0) THEN
    ! first work out w profile for forcing at this time
    CALL profile_w_force(w_subs, wdims, fld_type_w, w_force)

    CALL apply_w_forcing(tdims_s,pdims_s,tdims,w_force, rho, theta, temp_inc)
    ! - apply theta increment
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          theta_star(i,j,k) = theta_star(i,j,k) + temp_inc(i,j,k)
        END DO
      END DO
    END DO

    CALL apply_w_forcing(tdims_l,pdims_s,tdims,w_force, rho, q, temp_inc)
    ! - apply moisture increment
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_star(i,j,k) = q_star(i,j,k) + temp_inc(i,j,k)
        END DO
      END DO
    END DO
 
  END IF

! Increment forcing
! - temperature increment
  IF (theta_inc % n_times > 0) THEN
    CALL profile_increment(theta_inc, fld_type_w, tdims, temp_inc)

! - only apply forcing at lower boundary if level 0 is prognostic
    IF (l_prognostic_level0) THEN
      first_physics_level = tdims%k_start
    ELSE
      first_physics_level = 1
    END IF

! - apply temperature increment according to specified field type.
    SELECT CASE (theta_inc % field_type)

    CASE(theta_dry)
      DO k = first_physics_level, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            theta_star(i,j,k) = theta_star(i,j,k) + temp_inc(i,j,k)
          END DO
        END DO
      END DO

    CASE(temperature)
      DO k = first_physics_level, tdims%k_end
        DO j = tdims%j_start, tdims%j_end
          DO i = tdims%i_start, tdims%i_end
            theta_star(i,j,k) = theta_star(i,j,k) + temp_inc(i,j,k) /          &
                                                    exner_theta_levels(i,j,k)
          END DO
        END DO
      END DO

    CASE DEFAULT
      ErrorStatus = 1
      WRITE(Cmessage, '(''Unknown theta_inc_field_type = '', I0)')             &
                      theta_inc % field_type
      CALL Ereport(RoutineName, ErrorStatus, Cmessage)
    END SELECT

  END IF

! - moisture increment
  IF (mv_inc % n_times > 0) THEN
    CALL profile_increment(mv_inc, fld_type_w, tdims, temp_inc)

! - apply moisture increment
    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_star(i,j,k) = q_star(i,j,k) + temp_inc(i,j,k)
        END DO
      END DO
    END DO

  END IF

  ! Special case of localised heating
  IF (local_heat_option > omit) THEN
    
    CALL idl_local_heat(temp_inc)

    DO k = tdims%k_start, tdims%k_end
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          theta_star(i,j,k) = theta_star(i,j,k) + temp_inc(i,j,k)
        END DO
      END DO
    END DO
  END IF   ! test on local heating

END IF

IF (l_idealised_data) THEN
!-------------------------------------------------------------------------
! Work out model increments using the stored values from the beginning of 
! the routine. Note *_star values zeroed at the beginning of this routine so
! they will contain just the increments from this routine
DO k =             1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dtheta_inc_ideal_um(i,j,k) = theta_star(i,j,k)
      dt_inc_ideal_um(i,j,k) = dtheta_inc_ideal_um(i,j,k)              &
                                         * exner_theta_levels(i,j,k)
    END DO
  END DO
END DO

DO k =             1, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      du_inc_ideal_um(i,j,k) = R_u(i,j,k)
    END DO
  END DO
END DO
DO k =             1, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      dv_inc_ideal_um(i,j,k) = R_v(i,j,k)
    END DO
  END DO
END DO

DO k =             1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      dq_inc_ideal_um(i,j,k) = q_star(i,j,k)
    END DO
  END DO
END DO

  ! Calculate column integrals from forcing increments
  CALL idl_col_int_diag(l_mr_physics, u, v, rho, wetrho_r_sq_n)


END IF ! l_idealised_data


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_forcing

!======================================================================

SUBROUTINE eg_idl_forcing2(                                           &
  u, v,exner_theta_levels, exner, theta, q, rho, wetrho_r_sq_n,       &
  r_u, r_v, q_star, theta_star, qlimit,                               &
! error information
  error_code  )

USE idl_col_int_diag_mod,      ONLY: idl_col_int_diag
USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE atm_fields_bounds_mod,     ONLY: tdims, tdims_s, pdims, pdims_s, &
                                     udims, udims_s, vdims, vdims_s, &
                                     tdims_l
USE field_types
USE halo_exchange,             ONLY: swap_bounds

USE gen_phys_inputs_mod,       ONLY: l_mr_physics
USE timestep_mod,              ONLY: timestep
USE planet_constants_mod,      ONLY: r, cp
USE idealised_diag_mod,        ONLY: dt_inc_ideal_um, dq_inc_ideal_um, &
                                     du_inc_ideal_um, dv_inc_ideal_um, &
                                     dtheta_inc_ideal_um
USE eg_helmholtz_mod,          ONLY: ec_vol
USE horiz_grid_mod,            ONLY: intw_w2rho, cell_area_surface, intw_p2u,  &
                                     intw_p2v

IMPLICIT NONE

! Description: Initialises star state for ENDGame when there is no
!              physics. Optionally applies Held-Suarez test of
!              dynamical core (BAMS 75, 1825-1830).

REAL, INTENT (IN) ::                                                  &
  u(                 udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  v(                 vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end ,                  &
                     vdims_s%k_start:vdims_s%k_end),                  &
  exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
  exner(             pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end+1),                &
  theta(             tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end),                  &
       q(            tdims_l%i_start:tdims_l%i_end,                   &
                     tdims_l%j_start:tdims_l%j_end,                   &
                     tdims_l%k_start:tdims_l%k_end),                  &
       rho(          pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end),                  &
  wetrho_r_sq_n(     pdims_s%i_start:pdims_s%i_end,                   &
                     pdims_s%j_start:pdims_s%j_end,                   &
                     pdims_s%k_start:pdims_s%k_end+1)
REAL, INTENT (INOUT) ::                                               &
  r_u(               udims_s%i_start:udims_s%i_end,                   &
                     udims_s%j_start:udims_s%j_end,                   &
                     udims_s%k_start:udims_s%k_end),                  &
  r_v(               vdims_s%i_start:vdims_s%i_end,                   &
                     vdims_s%j_start:vdims_s%j_end,                   &
                     vdims_s%k_start:vdims_s%k_end),                  &
  q_star(            tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end),                      &
  theta_star(        tdims%i_start:tdims%i_end,                       &
                     tdims%j_start:tdims%j_end,                       &
                     tdims%k_start:tdims%k_end)
REAL ::                                                               &
    theta_in(        tdims_s%i_start:tdims_s%i_end,                   &
                     tdims_s%j_start:tdims_s%j_end,                   &
                     tdims_s%k_start:tdims_s%k_end)

INTEGER :: error_code


INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_IDL_FORCING2'

INTEGER :: i, j, k

REAL :: q_in (tdims_l%i_start:tdims_l%i_end,     &
           tdims_l%j_start:tdims_l%j_end,     &
           tdims_l%k_start:tdims_l%k_end)

REAL :: qlimit

REAL, ALLOCATABLE  ::                                 &
  th_temp(:,:,:), q_temp(:,:,:), u_temp(:,:,:), v_temp(:,:,:)

! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Copy star value before increments so that it is possible to work out
! total increments to fields from the idealised forcing. Note as 2 different
! calls to eg_idl_forcing it is possible that increments may not be zero
! so cannot use diagnostic arrays to hold copies 

! Allocate temporaray arrays
ALLOCATE(th_temp(tdims%i_start: tdims%i_end,    &
                 tdims%j_start: tdims%j_end,    &
                 tdims%k_start: tdims%k_end ) )
ALLOCATE( q_temp(tdims%i_start: tdims%i_end,    &
                 tdims%j_start: tdims%j_end,    &
                 tdims%k_start: tdims%k_end ) )
ALLOCATE( u_temp(udims%i_start: udims%i_end,    &
                 udims%j_start: udims%j_end,    &
                 udims%k_start: udims%k_end ) )
ALLOCATE( v_temp(vdims%i_start: vdims%i_end,    &
                 vdims%j_start: vdims%j_end,    &
                 vdims%k_start: vdims%k_end ) )

! copy star value before increments so that it is possible to work out
! total increments to fields from the idealised forcing

DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      th_temp(i,j,k) = theta_star(i,j,k)
    END DO
  END DO
END DO
DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      q_temp(i,j,k) = q_star(i,j,k)
    END DO
  END DO
END DO
DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      u_temp(i,j,k) = R_u(i,j,k)
    END DO
  END DO
END DO
DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      v_temp(i,j,k) = R_v(i,j,k)
    END DO
  END DO
END DO


!-------------------------------------------------------------------------
! Work out model increments using the stored values from the beginning of 
! the routine.
DO k = tdims%k_start , tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      th_temp(i,j,k)   = theta_star(i,j,k) - th_temp(i,j,k)
      dtheta_inc_ideal_um(i,j,k) = dtheta_inc_ideal_um(i,j,k)+           &
                                              th_temp(i,j,k)
      ! Convert increment to dT
      th_temp(i,j,k)   = th_temp(i,j,k)  *exner_theta_levels(i,j,k)
      dt_inc_ideal_um(i,j,k) = dt_inc_ideal_um(i,j,k)  +           &
                                      th_temp(i,j,k)
    END DO
  END DO
END DO

DO k = udims%k_start, udims%k_end
  DO j = udims%j_start, udims%j_end
    DO i = udims%i_start, udims%i_end
      du_inc_ideal_um(i,j,k) = du_inc_ideal_um(i,j,k)+           &
                               R_u(i,j,k) - u_temp(i,j,k)
    END DO
  END DO
END DO
DO k = vdims%k_start, vdims%k_end
  DO j = vdims%j_start, vdims%j_end
    DO i = vdims%i_start, vdims%i_end
      dv_inc_ideal_um(i,j,k) = dv_inc_ideal_um(i,j,k) +           &
                               R_v(i,j,k) - v_temp(i,j,k)
    END DO
  END DO
END DO

DO k = tdims%k_start, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      q_temp(i,j,k) = q_star(i,j,k) - q_temp(i,j,k)
      dq_inc_ideal_um(i,j,k) = dq_inc_ideal_um(i,j,k)+             &
                               q_temp(i,j,k)
    END DO
  END DO
END DO


DEALLOCATE( v_temp )
DEALLOCATE( u_temp )
DEALLOCATE( q_temp )
DEALLOCATE( th_temp )

! Calculate column integrals from forcing increments
CALL idl_col_int_diag(l_mr_physics, u, v, rho, wetrho_r_sq_n)

!----------------------------------------------------------------------------

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE eg_idl_forcing2

END MODULE eg_idl_forcing_mod
