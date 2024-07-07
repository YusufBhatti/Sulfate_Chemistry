! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE external_force_2_mod

IMPLICIT NONE
  ! Description:
  !
  ! Method: 
  !  
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: IDEALISED
  !
  ! Code Description:
  !   Language: FORTRAN 90
  !   This code is written to UMDP3 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='EXTERNAL_FORCE_2_MOD'

CONTAINS
SUBROUTINE external_force_2(l_spec_z0, flux_e, flux_h, ustar_in,               &
                            z0m_scm, z0h_scm)

USE parkind1,                  ONLY: jpim, jprb       !DrHook
USE yomhook,                   ONLY: lhook, dr_hook   !DrHook
USE Ereport_mod,               ONLY: Ereport
USE umPrintMgr,                ONLY: umPrint
USE errormessagelength_mod,    ONLY: errormessagelength
USE chk_opts_mod,              ONLY: chk_var, def_src
USE conversions_mod,           ONLY: pi
USE atm_fields_bounds_mod,     ONLY: pdims
USE nlsizes_namelist_mod,      ONLY: global_row_length, global_rows
USE horiz_grid_mod,            ONLY: glob_xi1_p, glob_xi2_p, xi1_p, xi2_p
USE model_time_mod,            ONLY: i_hour, i_minute, i_second
USE timestep_mod,              ONLY: timestep, timestep_number
USE missing_data_mod,          ONLY: rmdi
USE idealise_run_mod,          ONLY: roughlen_z0m, roughlen_z0h
USE surface_flux_mod,          ONLY: IdlSurfFluxSeaOption,                     &
                                     zero_flux, diurnal_flux, constant_flux,   &
                                     hot_spot, time_varying,                   &
                                     IdlSurfFluxSeaParams,                     &
                                     num_surface_flux_times, surface_flux_time,&
                                     sh_flux, lh_flux
USE bl_option_mod,             ONLY: flux_bc_opt, interactive_fluxes,          &
                                     specified_fluxes_only,                    &
                                     specified_fluxes_tstar,                   &
                                     specified_fluxes_cd

IMPLICIT NONE

! Description: 

! Switch indicating whether roughness length is specified
LOGICAL, INTENT(IN)  :: l_spec_z0

! Surface latent heat flux
REAL,    INTENT(OUT) :: flux_e(pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end)
! Surface sensible heat flux
REAL,    INTENT(OUT) :: flux_h(pdims%i_start:pdims%i_end,                      &
                               pdims%j_start:pdims%j_end)
! Surface friction velocity
REAL,    INTENT(OUT) :: ustar_in(pdims%i_start:pdims%i_end,                    &
                                 pdims%j_start:pdims%j_end)
! Roughness length m
REAL,    INTENT(OUT) :: z0m_scm(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end)
! Roughness length h
REAL,    INTENT(OUT) :: z0h_scm(pdims%i_start:pdims%i_end,                     &
                                pdims%j_start:pdims%j_end)


! Local variables
INTEGER :: i, j, k
REAL    :: x_centre, y_centre, x_length, y_length, xi, eta
REAL    :: hrl, xfact
REAL    ::   &
  t          & ! current time in seconds
 ,tau        & ! tau
 ,sh_now     & ! sensible heat this timestep
 ,lh_now       ! latent heat this timestep

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


INTEGER                           :: ErrorStatus
CHARACTER(LEN=errormessagelength) :: Cmessage
CHARACTER(LEN=*), PARAMETER       :: RoutineName='EXTERNAL_FORCE_2'

! 1.0 Start of subroutine code: perform the calculation.
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!--------------------------------------------------------------
!
!                 Idealised Surface Fluxes
!
!--------------------------------------------------------------
! Initialise surface fluxes to missing data
DO j = pdims%j_start, pdims%j_end
  DO i = pdims%i_start, pdims%i_end
    ustar_in(i,j) = rmdi
    flux_h(i,j)   = rmdi
    flux_e(i,j)   = rmdi
  END DO
END DO

SELECT CASE(flux_bc_opt)

CASE(specified_fluxes_only, specified_fluxes_tstar)
  ! Note that if specified_fluxes_tstar is selected, the surface 
  ! temperature, tstar, from the input dump will be used to generate 
  ! the surface upward LW flux.  If specified_fluxes_only, a value of 
  ! surface temperature will be calculated to be approximately 
  ! consistent with the specified flux and local atmospheric state.

  ! Overwrite sea surface scalar fluxes with idealised fluxes
  SELECT CASE(IdlSurfFluxSeaOption)

  !---------------------------------
  ! Option 1: Zero surface fluxes
  !---------------------------------
  CASE(zero_flux)

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        flux_h(i,j) = 0.0
        flux_e(i,j) = 0.0
      END DO
    END DO

  !---------------------------------
  ! Option 2: Diurnal cycle
  ! (positive surface fluxes during the day, zero at night)
  !---------------------------------
  CASE(diurnal_flux)

    ! IdlSurfFluxSeaParams(1) = max sensible heat flux
    ! IdlSurfFluxSeaParams(2) = max latent heat flux
    ! IdlSurfFluxSeaParams(3) = Time (UTC) of max flux (hours)
    ! IdlSurfFluxSeaParams(4) = Length of day (hours)

    ! Calculate current time (UTC) in hours
    hrl   = ((i_hour*60.0 + i_minute)*60.0 + i_second)/3600.0

    ! Set up diurnally varying function
    xfact  = COS( 0.5*pi*(IdlSurfFluxSeaParams(3) - hrl) /                     &
                         (IdlSurfFluxSeaParams(4) / 2.0)  )

    ! Limit fluxes to being positive (upward)
    IF (xfact <= 0.0) xfact=0.0

    ! Set diurnally varying fluxes
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        flux_h(i,j) = IdlSurfFluxSeaParams(1)*xfact**1.5
        flux_e(i,j) = IdlSurfFluxSeaParams(2)*xfact**1.3
      END DO
    END DO

  !---------------------------------
  ! Option 3: Fixed surface fluxes
  !---------------------------------
  CASE(constant_flux)

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        flux_h(i,j) = IdlSurfFluxSeaParams(1)
        flux_e(i,j) = IdlSurfFluxSeaParams(2)
      END DO
    END DO

  CASE(hot_spot)

    x_centre = glob_xi1_p(1) + 0.5 * ( glob_xi1_p(global_row_length) -         &
                                       glob_xi1_p(1))
    y_centre = glob_xi2_p(1) + 0.5 * ( glob_xi2_p(global_rows) -               &
                                       glob_xi2_p(1) )

    x_length = xi1_p(2) - xi1_p(1)
    y_length = xi2_p(2) - xi2_p(1)

    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        xi  = (xi1_p(i) - x_centre) / x_length
        eta = (xi2_p(j) - y_centre) / y_length
        IF ( ABS(xi) <= 2.1 .AND. ABS(eta) <= 2.1) THEN
          flux_h(i,j) = IdlSurfFluxSeaParams(3)
          flux_e(i,j) = IdlSurfFluxSeaParams(4)
        ELSE
          flux_h(i,j) = IdlSurfFluxSeaParams(1)
          flux_e(i,j) = IdlSurfFluxSeaParams(2)
        END IF
      END DO
    END DO

  CASE(time_varying)
    ! Sensible and latent heat vary with time but the same value everywhere

    ! Calculate current time from start of simulation in seconds
    t = timestep * timestep_number
    IF ( t < surface_flux_time(1) ) THEN
      ! If before first time use flux at first time
      sh_now = sh_flux(1)  
      lh_now = lh_flux(1)  
    ELSE IF (t >= surface_flux_time(num_surface_flux_times)) THEN
      ! If after last time use  flux at last time
      sh_now = sh_flux(num_surface_flux_times)  
      lh_now = lh_flux(num_surface_flux_times)  
    ELSE 
      ! interpolate to current time assumes times are in ascending order
      DO i = 1, num_surface_flux_times - 1
        IF (surface_flux_time(i+1) > t .AND. t >= surface_flux_time(i)) THEN
          tau = (t - surface_flux_time(i)) /                                   &
                (surface_flux_time(i+1) - surface_flux_time(i))
          sh_now = (1.0 - tau) * sh_flux(i) + tau * sh_flux(i+1)
          lh_now = (1.0 - tau) * lh_flux(i) + tau * lh_flux(i+1)
        END IF
      END DO
    END IF

    ! Set varying fluxes
    DO j = pdims%j_start, pdims%j_end
      DO i = pdims%i_start, pdims%i_end
        flux_h(i,j) = sh_now
        flux_e(i,j) = lh_now
      END DO
    END DO

  CASE DEFAULT
    ErrorStatus = 1
    WRITE(Cmessage, '(''Unknown IdlSurfFluxSeaOption = '', I0)')               &
                    IdlSurfFluxSeaOption
    CALL Ereport(RoutineName, ErrorStatus, Cmessage)

  END SELECT

CASE (specified_fluxes_cd)

  ErrorStatus = 1
  WRITE(Cmessage,'(A)') 'Code not yet written to specify ustar in '//          &
                        'the idealised UM'
  CALL Ereport(RoutineName, ErrorStatus, Cmessage)

END SELECT

!--------------------------------------------------------------
!
!                 Idealised Roughness Length
!
!--------------------------------------------------------------

IF (l_spec_z0) THEN
  DO j = pdims%j_start, pdims%j_end
    DO i = pdims%i_start, pdims%i_end
      z0m_scm(i,j) = roughlen_z0m
      z0h_scm(i,j) = roughlen_z0h
    END DO
  END DO
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,                  &
                        zhook_out,zhook_handle)

RETURN
END SUBROUTINE external_force_2

END MODULE external_force_2_mod
