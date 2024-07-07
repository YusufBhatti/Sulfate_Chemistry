! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Set thermodynamic properties
!
! Purpose:
!   Layer masses, densities, heat capacities, and temperatures at the
!   edges of layers are set.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Radiation Control
!
!------------------------------------------------------------------------------
MODULE set_thermodynamic_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SET_THERMODYNAMIC_MOD'

CONTAINS

SUBROUTINE set_thermodynamic(row_length, rows, off_x, off_y, t_layer, rho_r2,  &
  p_layer_centres, q, qcl, qcf, qcf2, qrain, qgraup,                           &
  l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,                                       &
  p_layer_boundaries, t_layer_boundaries, p_extra_layer, t_extra_layer,        &
  d_mass, density, layer_heat_capacity)

USE nlsizes_namelist_mod,   ONLY: model_levels
USE planet_constants_mod,   ONLY: l_planet_g, g, cp, r, rv, c_virtual
USE gravity_mod,            ONLY: g_theta
USE water_constants_mod,    ONLY: hcapv, hcapw, hcapi
USE rad_input_mod,          ONLY: l_hydrostatic_mass, l_moist_heat_capacity, &
                                  l_extra_top
USE level_heights_mod,      ONLY: r_theta_levels, r_rho_levels
USE science_fixes_mod,      ONLY: l_fail_p_layers_inconsis
USE gen_phys_inputs_mod,    ONLY: l_mr_physics
USE ereport_mod,            ONLY: ereport
USE errormessagelength_mod, ONLY: errormessagelength
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE

INTEGER, INTENT(IN) :: row_length, rows, off_x, off_y
!   Model grid

REAL, INTENT(IN) :: t_layer(row_length, rows, model_levels)
!   Temperature of layers
REAL, INTENT(IN) :: &
  rho_r2(1-off_x:row_length+off_x, 1-off_y:rows+off_y, model_levels)
!   Air density*radius**2 / kg m-1
REAL, INTENT(IN) :: p_layer_centres(row_length, rows, 0:model_levels)
!   Pressure at layer centres
REAL, INTENT(IN) :: q(row_length, rows, model_levels)
!   Vapour content / kg kg-1
REAL, INTENT(IN) :: qcl(row_length, rows, model_levels)
!   Liquid water content / kg kg-1
REAL, INTENT(IN) :: qcf(row_length, rows, model_levels)
!   Ice content / kg kg-1
REAL, INTENT(IN) :: qcf2(row_length, rows, model_levels)
!   Second ice content / kg kg-1
REAL, INTENT(IN) :: qrain(row_length, rows, model_levels)
!   Rain water content / kg kg-1
REAL, INTENT(IN) :: qgraup(row_length, rows, model_levels)
!   Graupel content / kg kg-1

LOGICAL, INTENT(IN) :: l_mcr_qcf2
!   Use second prognostic ice
LOGICAL, INTENT(IN) :: l_mcr_qrain
!   Use prognostic rain
LOGICAL, INTENT(IN) :: l_mcr_qgraup
!   Use prognostic graupel

REAL, INTENT(INOUT) :: p_layer_boundaries(row_length, rows, 0:model_levels)
!   Pressure at layer boundaries

REAL, INTENT(OUT) :: t_layer_boundaries(row_length, rows, 0:model_levels)
!   Temperature at layer boundaries
REAL, INTENT(OUT) :: p_extra_layer(row_length, rows)
!   Pressure at centre of extra top layer
REAL, INTENT(OUT) :: t_extra_layer(row_length, rows)
!   Temperature at centre of extra top layer
REAL, INTENT(OUT) :: d_mass(row_length, rows, model_levels+1)
!   Mass of layer (kg m-2)
REAL, INTENT(OUT) :: density(row_length, rows, model_levels+1)
!   Density of layer (kg m-3)
REAL, INTENT(OUT) :: layer_heat_capacity(row_length, rows, model_levels)
!   Heat capacity of layer

INTEGER :: i, j, k
!   Loop variables
REAL :: rho(row_length, rows, model_levels)
!   Air density on rho levels / kg m-3
REAL :: q_total(row_length, rows, model_levels)
!   Total moisture / kg kg-1
REAL :: weight_lower, weight_upper
!   Weights for lower and upper layers
REAL :: rm
!   Gas constant for moist air

CHARACTER (LEN=*), PARAMETER       :: RoutineName = 'SET_THERMODYNAMIC'
CHARACTER (LEN=errormessagelength) :: cmessage
INTEGER                            :: n_points_negative_mass
INTEGER                            :: ErrorStatus

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Calculate the layer masses and densities
! ----------------------------------------------------------------------

IF (l_hydrostatic_mass) THEN
  ! ----------------------------------------------------------------------  
  ! Assuming hydrostatic equilibrium
  ! ----------------------------------------------------------------------  
  n_points_negative_mass = 0
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP& SHARED( model_levels, rows, row_length, l_planet_g, g_theta, d_mass, &
!$OMP&         p_layer_boundaries, g, l_mr_physics, &
!$OMP&         r, q, density, p_layer_centres, t_layer, c_virtual ) &
!$OMP& PRIVATE(i,j,k,rm) &
!$OMP& REDUCTION(+:n_points_negative_mass)
  DO k = 1, model_levels - 1
    DO j = 1, rows
      DO i = 1, row_length
        IF (l_planet_g .AND. ALLOCATED(g_theta)) THEN
          d_mass(i, j, k) = &
            ( p_layer_boundaries(i,j,k-1) - p_layer_boundaries(i,j,k) ) / &
            g_theta(i,j,k)
        ELSE
          d_mass(i, j, k) = &
            ( p_layer_boundaries(i,j,k-1) - p_layer_boundaries(i,j,k) ) / g
        END IF
        IF (d_mass(i, j, k) <= 0.0) THEN
          n_points_negative_mass = n_points_negative_mass + 1
          d_mass(i, j, k) = ABS(d_mass(i, j, k))
        END IF
        IF (l_mr_physics) THEN
          rm = ( r + rv * q(i,j,k) ) / ( 1.0 + q(i,j,k) )
          density(i, j, k) = p_layer_centres(i, j, k) / (rm * t_layer(i, j, k))
        ELSE
          density(i, j, k) = p_layer_centres(i, j, k) &
            / (r * t_layer(i, j, k) * ( 1.0 + c_virtual * q(i,j,k) ))
        END IF
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  IF (n_points_negative_mass > 0) THEN
    ! Prepare a meaningful error message/warning:
    WRITE(cmessage, '(A,I5,A)') 'A total of ', n_points_negative_mass, &
      ' points had negative mass in set_thermodynamic.'         // &
      ' This indicates the pressure fields are inconsistent'    // &
      ' between different levels and the model is about to fail.'
  
    ! What happens here depends on l_fail_p_layers_inconsis
    IF (l_fail_p_layers_inconsis) THEN
      ! This is a full error, need to fail.
      ErrorStatus = 100
    ELSE
      ! raise a warning (the UM will probably fail soon anyway!)
      ErrorStatus = -100
    END IF
    CALL ereport('set_thermodynamic', ErrorStatus, cmessage)
  END IF

ELSE
  ! ----------------------------------------------------------------------
  ! Use air densities interpolated from rho levels
  ! ----------------------------------------------------------------------
  ! We should note that this formulation, although better than the
  ! hydrostatic formulation, is still not entirely conservative.

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        ! Calculate densities by removing the r**2 term from rho_r2.
        rho(i,j,k) = rho_r2(i,j,k) &
          / ( r_rho_levels(i,j,k) * r_rho_levels(i,j,k) )
      END DO
    END DO
  END DO
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP& SHARED( model_levels, rows, row_length, d_mass, rho, r_theta_levels, &
!$OMP&         r_rho_levels, density ) &
!$OMP& PRIVATE(i,j,k)
  DO k = 1, model_levels - 1
    DO j = 1, rows
      DO i = 1, row_length
  
        ! Calculate the average value of rho across the layer
        ! multiplied by the layer thickness, and then the average density.
        d_mass(i,j,k) = &
          rho(i,j,k+1) * ( r_theta_levels(i,j,k) - r_rho_levels(i,j,k) ) + &
          rho(i,j,k)   * ( r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k) )
        density(i,j,k) = d_mass(i,j,k) &
          / ( r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k) )
  
        IF (k == 1) THEN
          ! For the lowest layer we need to extend the lower
          ! boundary from the first rho level to the surface.
          ! The surface is the 0'th theta level.
          d_mass(i,j,1) = density(i,j,1) &
            * ( r_rho_levels(i,j,2) - r_theta_levels(i,j,0) )
        END IF  ! k  ==  1
  
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

END IF


! ----------------------------------------------------------------------
! Adjust the top layer for an extra layer added above
! ----------------------------------------------------------------------

IF (l_extra_top) THEN
  ! Here, an extra level is inserted at the top of the atmosphere by replacing
  ! the top layer used in the rest of the model with two layers. The masses are
  ! set by assuming the hydrostatic approximation.

  ! Note that we will continue to use this calculation for the mass of the
  ! artificial top layer when the more accurate calculation (non-hydostatic)
  ! is used for the lower levels because we do not have the model information
  ! to make a better estimate.
  k = model_levels
  DO j = 1, rows
    DO i = 1, row_length
      ! Extrapolate pressure
      p_layer_boundaries(i,j,k) = 10**( &
        LOG10(p_layer_centres(i,j,k)) - LOG10(p_layer_boundaries(i,j,k-1)) + &
        LOG10(p_layer_centres(i,j,k)) )
      p_extra_layer(i,j) = 0.5 * p_layer_boundaries(i,j,k)
      ! Extrapolate temperature in log10(pressure)
      t_extra_layer(i,j) = (t_layer(i,j,k) - t_layer(i,j,k-1))/ &
        (LOG10(p_layer_centres(i,j,k)) - LOG10(p_layer_centres(i,j,k-1)))
      t_layer_boundaries(i,j,k) = &
        t_extra_layer(i,j) * ( LOG10(p_layer_boundaries(i,j,k)) &
        - LOG10(p_layer_centres(i,j,k-1)) ) + t_layer(i,j,k-1)
      t_extra_layer(i,j) = &
        t_extra_layer(i,j) * ( LOG10(p_extra_layer(i,j)) &
        - LOG10(p_layer_centres(i,j,k-1)) ) + t_layer(i,j,k-1)
      IF (l_planet_g .AND. ALLOCATED(g_theta)) THEN
        d_mass(i, j, k) = &
          (p_layer_boundaries(i,j,k-1) - p_layer_boundaries(i,j,k)) / &
          g_theta(i,j,k)
        d_mass(i, j, k+1) = p_layer_boundaries(i,j,k) / &
          g_theta(i,j,k)
      ELSE
        d_mass(i, j, k) = &
          (p_layer_boundaries(i,j,k-1) - p_layer_boundaries(i,j,k)) / g
        d_mass(i, j, k+1) = p_layer_boundaries(i,j,k) / g
      END IF
      IF (l_mr_physics) THEN
        rm = ( r + rv * q(i,j,k) ) / ( 1.0 + q(i,j,k) )
        density(i, j, k) = p_layer_centres(i, j, k) / (rm * t_layer(i, j, k))
        density(i, j, k+1) = p_extra_layer(i, j) / (rm * t_extra_layer(i, j))
      ELSE
        density(i, j, k) = p_layer_centres(i, j, k) &
          / (r * t_layer(i, j, k) * ( 1.0 + c_virtual * q(i,j,k) ))
        density(i, j, k+1) = p_extra_layer(i, j) &
          / (r * t_extra_layer(i, j) * ( 1.0 + c_virtual * q(i,j,k) ))
      END IF
    END DO
  END DO
ELSE
  k = model_levels
  DO j = 1, rows
    DO i = 1, row_length
      ! Top layer mass represents the total mass above.
      ! Use the hydostatic approximation as other information not available.
      IF (l_planet_g .AND. ALLOCATED(g_theta)) THEN
        d_mass(i, j, k) = p_layer_boundaries(i,j,k-1) / g_theta(i,j,k)
      ELSE
        d_mass(i, j, k) = p_layer_boundaries(i,j,k-1) / g
      END IF
      IF (l_mr_physics) THEN
        rm = ( r + rv * q(i,j,k) ) / ( 1.0 + q(i,j,k) )
        density(i, j, k) = p_layer_centres(i, j, k) / (rm * t_layer(i, j, k))
      ELSE
        density(i, j, k) = p_layer_centres(i, j, k) &
          / (r * t_layer(i, j, k) * ( 1.0 + c_virtual * q(i,j,k) ))
      END IF
    END DO
  END DO
END IF


! ----------------------------------------------------------------------
! Set the layer heat capacities.
! Convert masses and densities to dry values when using mixing ratios.
! ----------------------------------------------------------------------

IF (l_moist_heat_capacity .OR. l_mr_physics) THEN
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP& SHARED( model_levels, rows, row_length, q_total, q, qcl, qcf, &
!$OMP&         layer_heat_capacity, l_mcr_qcf2, qcf2, l_mcr_qrain, qrain, &
!$OMP&         l_mcr_qgraup, qgraup, l_mr_physics, l_moist_heat_capacity, &
!$OMP&         d_mass, cp, density ) &
!$OMP& PRIVATE(i,j,k)
  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
    
        ! Calculate total moisture
        q_total(i,j,k) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)
  
        ! Calculate the specific heat capacity of the moist air
        layer_heat_capacity(i,j,k) = q(i,j,k)*hcapv +  &
          qcl(i,j,k)*hcapw + qcf(i,j,k)*hcapi
  
        ! Add on contributions from optional prognostics
        IF (l_mcr_qcf2) THEN
          q_total(i,j,k) = q_total(i,j,k) + qcf2(i,j,k)
          layer_heat_capacity(i,j,k) = layer_heat_capacity(i,j,k) &
            + qcf2(i,j,k)*hcapi
        END IF
  
        IF (l_mcr_qrain) THEN
          q_total(i,j,k) = q_total(i,j,k) + qrain(i,j,k)
          layer_heat_capacity(i,j,k) = layer_heat_capacity(i,j,k) &
            + qrain(i,j,k)*hcapw
        END IF
  
        IF (l_mcr_qgraup) THEN
          q_total(i,j,k) = q_total(i,j,k) + qgraup(i,j,k)
          layer_heat_capacity(i,j,k) = layer_heat_capacity(i,j,k) &
            + qgraup(i,j,k)*hcapi
        END IF
  
  
        ! Calculate heat capacity of the layer.
        IF (l_mr_physics) THEN
          IF (.NOT. l_moist_heat_capacity) THEN
            layer_heat_capacity(i, j, k) = d_mass(i, j, k)*cp
          END IF
          ! d_mass and density represent dry values
          d_mass(i,j,k)  = d_mass(i,j,k)  / (1.0 + q_total(i,j,k))
          density(i,j,k) = density(i,j,k) / (1.0 + q_total(i,j,k))
          IF (l_moist_heat_capacity) THEN
            layer_heat_capacity(i,j,k) = &
              d_mass(i,j,k) * ( cp + layer_heat_capacity(i,j,k) )
          END IF
        ELSE
          ! d_mass and density represents moist values
          layer_heat_capacity(i,j,k) = d_mass(i,j,k) * &
            ( (1.0 - q_total(i,j,k))*cp + layer_heat_capacity(i,j,k) )
        END IF
  
      END DO
    END DO
  END DO
!$OMP END PARALLEL DO

  IF (l_extra_top) THEN
    IF (l_mr_physics) THEN
      k = model_levels
      DO j = 1, rows
        DO i = 1, row_length
          d_mass(i,j,k+1)  = d_mass(i,j,k+1)  / (1.0 + q_total(i,j,k))
          density(i,j,k+1) = density(i,j,k+1) / (1.0 + q_total(i,j,k))
        END DO
      END DO
    END IF
  END IF

ELSE

  DO k = 1, model_levels
    DO j = 1, rows
      DO i = 1, row_length
        layer_heat_capacity(i, j, k) = d_mass(i, j, k)*cp
      END DO
    END DO
  END DO

END IF


! ----------------------------------------------------------------------
! Set the boundary temperatures
! ----------------------------------------------------------------------

! Isothermal extrapolation at the surface
DO j = 1, rows
  DO i = 1, row_length
    t_layer_boundaries(i, j, 0) = t_layer(i, j, 1)
  END DO
END DO
IF (.NOT. l_extra_top) THEN
  DO j = 1, rows
    DO i = 1, row_length
      ! Isothermal extrapolation at top-of-atmosphere
      t_layer_boundaries(i, j, model_levels) = t_layer(i, j, model_levels)
    END DO
  END DO
END IF

! Now set the temperatures at interior boundaries of the profile supplied
! using interpolation in height.
!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(NONE) &
!$OMP& SHARED( model_levels, rows, row_length, r_rho_levels, r_theta_levels, &
!$OMP&         t_layer_boundaries, t_layer ) &
!$OMP& PRIVATE(i,j,k,weight_upper,weight_lower)
DO k = 1,  model_levels-1
  DO j = 1, rows
    DO i = 1, row_length
      weight_upper = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k)
      weight_lower = r_theta_levels(i,j,k+1) - r_rho_levels(i,j,k+1)
      t_layer_boundaries(i, j, k) = &
        ( weight_upper*t_layer(i,j,k+1) + weight_lower*t_layer(i,j,k) ) &
        / (weight_lower + weight_upper)
    END DO
  END DO
END DO
!$OMP END PARALLEL DO


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE set_thermodynamic

END MODULE set_thermodynamic_mod
