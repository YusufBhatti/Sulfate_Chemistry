! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE eg_idl_friction_mod

! Modules for setups
USE held_suarez_forcing_mod, ONLY: held_suarez_friction

IMPLICIT NONE

  ! Description: 
  !  Module containing subroutines to assign, and employ
  !  selected Rayleigh friction (i.e. paramterised boundary
  !  layer drag)
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EG_IDL_FRICTION_MOD'

CONTAINS

  SUBROUTINE eg_idl_friction(u, v, exner_theta_levels, exner,           &
                             theta_star, q_star, p_theta_levels,        &
                             r_u, r_v, error_code, l_mixing_ratio)

    ! Description: This applies parameterised boundary layer
    !              drag/friction to the horizontal winds
    !              for example (rayleigh friction).
    !              fric_number controls the scheme stated in fric_mod

    USE jules_sea_seaice_mod,      ONLY: beta_evap
    USE bl_option_mod,             ONLY: fric_number, fric_none,             &
                                         fric_HeldSuarez_1, fric_HeldSuarez_2
    USE level_heights_mod,         ONLY: r_theta_levels
    USE timestep_mod,              ONLY: timestep
    USE atm_fields_bounds_mod,     ONLY: udims_s, vdims_s, tdims_s, pdims_s, &
                                         tdims, pdims

    !Redirect routine names to avoid clash with existing qsat routines
    USE qsat_mod, ONLY: qsat_new     => qsat,                                 &
                        qsat_mix_new => qsat_mix,                             &
                        l_new_qsat_ideal !Currently defaults to FALSE

    USE umPrintMgr,                ONLY: umPrint, umMessage
    USE ereport_mod,               ONLY: ereport
    USE parkind1,                  ONLY: jpim, jprb       !DrHook
    USE yomhook,                   ONLY: lhook, dr_hook   !DrHook

    IMPLICIT NONE

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    CHARACTER(LEN=*), PARAMETER :: RoutineName='EG_IDL_FRICTION'

    ! Subroutine arguments
    LOGICAL, INTENT(IN) :: l_mixing_ratio
    REAL, INTENT(IN) :: theta_star(tdims%i_start:tdims%i_end,                &
                                   tdims%j_start:tdims%j_end,                &
                                   tdims%k_start:tdims%k_end)
    REAL, INTENT(IN) :: p_theta_levels(tdims_s%i_start:tdims_s%i_end,        &
                                       tdims_s%j_start:tdims_s%j_end,        &
                                       tdims_s%k_start:tdims_s%k_end)
    REAL, INTENT (INOUT) ::                                                  &
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
                            vdims_s%j_start:vdims_s%j_end,                   &
                            vdims_s%k_start:vdims_s%k_end),                  &
         exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                   &
                            tdims_s%j_start:tdims_s%j_end,                   &
                            tdims_s%k_start:tdims_s%k_end),                  &
         exner(             pdims_s%i_start:pdims_s%i_end,                   &
                            pdims_s%j_start:pdims_s%j_end,                   &
                            pdims_s%k_start:pdims_s%k_end+1)
    REAL, INTENT(INOUT) :: q_star    (tdims%i_start:tdims%i_end,             &
                                      tdims%j_start:tdims%j_end,             &
                                      tdims%k_start:tdims%k_end)
    INTEGER, INTENT(INOUT) :: error_code

    REAL :: t_lev1(tdims%i_start:tdims%i_end,                                &
                   tdims%j_start:tdims%j_end)     ! t at level 1
    REAL :: qsat_lev1(tdims%i_start:tdims%i_end,                             &
                      tdims%j_start:tdims%j_end)  ! qsat at level 1
    REAL :: q_source(tdims%i_start:tdims%i_end,                              &
                     tdims%j_start:tdims%j_end)  ! q source
    REAL :: ch_u      ! ch*u parameter for moisture source

    INTEGER :: i,j, err_code

    ! 1.0 Start of subroutine code: perform the calculation.
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! ----------------------------------------------------------------------
    ! Section 1. Select the Scheme
    ! ----------------------------------------------------------------------
    IF (fric_number == fric_none) THEN
       ! No drag applied
    ELSE IF (fric_number== fric_HeldSuarez_1 .OR.                            &
         fric_number== fric_HeldSuarez_2) THEN
       CALL held_suarez_friction(u, v, exner_theta_levels, exner,            &
                                 r_u, r_v, fric_number)
    ELSE
       WRITE(umMessage,'(A,I3)') '-Friction/Drag Scheme not supported:',     &
            fric_number
       CALL umPrint(umMessage,src='eg_idl_friction_mod.F90')
       err_code = 1
       CALL ereport("eg_idl_friction_mod", err_code,                         &
            "Friction Scheme not supported" )
    END IF

    IF (beta_evap > 0.0) THEN
    ! add a simple source of moisture
      ! calculate level 1 temperature
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          t_lev1(i,j) = theta_star(i,j,1) * exner_theta_levels(i,j,1)
        END DO
      END DO
      ! calculate level 1 saturation specific humidity

      IF (l_new_qsat_ideal) THEN
        IF (l_mixing_ratio) THEN
          CALL qsat_mix_new(qsat_lev1,t_lev1,                                 &
                            p_theta_levels(pdims%i_start:pdims%i_end,         &
                                           pdims%j_start:pdims%j_end,1),      &
                            pdims%i_len,pdims%j_len)
        ELSE
          CALL qsat_new(qsat_lev1,t_lev1,                                     &
                            p_theta_levels(pdims%i_start:pdims%i_end,         &
                                           pdims%j_start:pdims%j_end,1),      &
                            pdims%i_len,pdims%j_len)
        END IF
      ELSE
        ! DEPENDS ON: qsat_mix
        CALL qsat_mix(qsat_lev1,t_lev1,                                       &
        p_theta_levels(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,1),&
                      pdims%i_len*pdims%j_len,l_mixing_ratio)
      END IF

      ! flux = beta * ch_u * (qsat - q)
      ! hence source dq = beta * ch_u * (qsat - q) * dt / dz
      ! take ch = 0.001, u = 10 m/s
      ch_u = 0.001 * 10
      ! beta_evap is moisture availability parameter, 1=sea, 0=none
      ! calculate the moisture source - do this implicitly otherwise it will
      ! blow up
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_source(i,j) = ( q_star(i,j,1) + beta_evap * ch_u *             &
                            qsat_lev1(i,j) * timestep /                    &
                   ( r_theta_levels(i,j,1) - r_theta_levels(i,j,0) ) ) /   &
                          ( 1 + beta_evap * ch_u * timestep /              &
                   ( r_theta_levels(i,j,1) - r_theta_levels(i,j,0) ) ) -   &
                   q_star(i,j,1)
        END DO
      END DO
      ! add it to the lowest level q
      DO j = tdims%j_start, tdims%j_end
        DO i = tdims%i_start, tdims%i_end
          q_star(i,j,0) = q_star(i,j,0) + q_source(i,j)
          q_star(i,j,1) = q_star(i,j,1) + q_source(i,j)
        END DO
      END DO
    END IF


    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE eg_idl_friction
END MODULE eg_idl_friction_mod
