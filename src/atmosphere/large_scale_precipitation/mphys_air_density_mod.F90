! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Microphysics precipitation scheme. Air density calculation

MODULE mphys_air_density_mod
! Description:
! Calculates air density in one place for use within the microphysics

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large Scale Precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='MPHYS_AIR_DENSITY_MOD'

CONTAINS

! Subroutine Interface:
SUBROUTINE mphys_air_density( rho_r2, q, qcl, qcf, qcf2, qrain, qgraup,       &
                              rhodz_dry, rhodz_moist, deltaz )

USE atm_fields_bounds_mod, ONLY: tdims, pdims_s

USE level_heights_mod,     ONLY: r_rho_levels, r_theta_levels

USE mphys_inputs_mod,      ONLY: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup

USE gen_phys_inputs_mod,   ONLY: l_mr_physics

! Dr Hook Modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

REAL, INTENT(IN)  :: rho_r2( pdims_s%i_start : pdims_s%i_end,                 &
                             pdims_s%j_start : pdims_s%j_end,                 &
                             pdims_s%k_start : pdims_s%k_end )
                             ! IN Air density * earth radius**2

REAL, INTENT(IN) :: q( tdims%i_start:tdims%i_end,                             &
                       tdims%j_start:tdims%j_end,                             &
                                   1:tdims%k_end )

REAL, INTENT(IN) :: qcf( tdims%i_start:tdims%i_end,                           &
      tdims%j_start:tdims%j_end,                                              &
      1:tdims%k_end )
                                       ! Cloud ice (kg per kg air).
REAL, INTENT(IN) :: qcl( tdims%i_start:tdims%i_end,                           &
      tdims%j_start:tdims%j_end,                                              &
      1:tdims%k_end )
                                       ! Cloud liquid water (kg per
                                       ! kg air)
REAL, INTENT(IN) ::  qcf2( tdims%i_start:tdims%i_end,                         &
       tdims%j_start:tdims%j_end,                                             &
       1:tdims%k_end )
                                       ! Ice (kg per kg air)
REAL, INTENT(IN) ::  qrain( tdims%i_start:tdims%i_end,                        &
        tdims%j_start:tdims%j_end,                                            &
        1:tdims%k_end )
                                       ! Rain (kg per kg air)
REAL, INTENT(IN) ::  qgraup( tdims%i_start:tdims%i_end,                       &
         tdims%j_start:tdims%j_end,                                           &
         1:tdims%k_end )

REAL, INTENT(OUT) ::  rhodz_dry( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                 &
                                             1 : tdims%k_end )
                                       ! Dry density

REAL, INTENT(OUT) ::  rhodz_moist( tdims%i_start : tdims%i_end,               &
                                   tdims%j_start : tdims%j_end,               &
                                               1 : tdims%k_end )
                                       ! Moist density

REAL, INTENT(OUT) ::  deltaz( tdims%i_start : tdims%i_end,                    &
                              tdims%j_start : tdims%j_end,                    &
                                          1 : tdims%k_end )

                                       ! Layer thickness

!------------------------------------------------------------------------------
! Purpose:
!   Calculates air density in a single place, for use in either the 3D or
!   CASIM microphysics
!   Documentation: UMDP 26.

!------------------------------------------------------------------------------
! Subroutine Arguments


!------------------------------------------------------------------------------
! Local Variables

REAL :: rho1( tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end )
REAL :: rho2( tdims%i_start : tdims%i_end,                                    &
              tdims%j_start : tdims%j_end )

REAL :: q_total( tdims%i_start : tdims%i_end,                                 &
                 tdims%j_start : tdims%j_end )

INTEGER :: i, j, k ! Loop counters

! Declarations for Dr Hook:
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle
CHARACTER(LEN=*), PARAMETER   :: RoutineName='MPHYS_AIR_DENSITY'

!==============================================================================
! Start of calculations
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Calculate the (non-hydrostatic) layer thicknesses (deltaz) and air
! densities multiplied by deltaz (rhodz_moist and rhodz_dry).
! ----------------------------------------------------------------------
! We should note that this formulation, although better than the
! hydrostatic formulation, is still not entirely conservative. To ensure
! conservation we would need to rewrite the large-scale precipitation
! scheme to consider masses in terms of rho<q>, and
! not the current <rho>q formulation.

! We only need to calculate averages for the moist levels

!$OMP  PARALLEL DO SCHEDULE(STATIC) DEFAULT(SHARED) PRIVATE(i,j,k,      &
!$OMP& q_total,rho1,rho2)
DO k = 1, tdims%k_end
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end

            ! Calculate densities at the boundaries of the layer
            ! by removing the r**2 term from rho_r2.
            ! Rho1 is the density at the lower boundary.
      rho1(i,j)= rho_r2(i,j,k)/( r_rho_levels(i,j,k) *                        &
                             r_rho_levels(i,j,k) )

            ! Check whether there is a rho level above the current
            ! moist level.
      IF ( k  <   tdims%k_end ) THEN
              ! Rho2 is the density at the upper boundary.
        rho2(i,j)= rho_r2(i,j,k+1)/( r_rho_levels(i,j,k+1) *                  &
                               r_rho_levels(i,j,k+1) )

              ! Calculate the average value of rho across the layer
              ! multiplied by the layer thickness and the layer
              ! thickness.
        rhodz_moist(i,j,k) =                                                  &
                      rho2(i,j) * ( r_theta_levels(i,j,k) -                   &
                                         r_rho_levels(i,j,k) )                &
                   +  rho1(i,j) * ( r_rho_levels(i,j,k+1) -                   &
                                    r_theta_levels(i,j,k) )
        deltaz(i,j,k) = r_rho_levels(i,j,k+1)                                 &
                          - r_rho_levels(i,j,k)

        IF (k  ==  1) THEN
          ! For the lowest layer we need to extend the lower
          ! boundary from the first rho level to the surface.
          ! The surface is the 0'th theta level.
          deltaz(i,j,1) = r_rho_levels(i,j,2)                                 &
                            - r_theta_levels(i,j,0)
          rhodz_moist(i,j,1) = rhodz_moist(i,j,1)*deltaz(i,j,1)               &
                   / (r_rho_levels(i,j,2)-r_rho_levels(i,j,1))
        END IF  ! k  ==  1

      ELSE
        ! For a top layer higher than the highest rho level
        ! we can calculate a pseudo rho level. We will assume
        ! it has a similar density to the rho level below
        ! and that the intervening theta level is in the centre
        ! of the layer.
        deltaz(i,j,k) = 2.0*(r_theta_levels(i,j,k)                            &
                              -r_rho_levels(i,j,k))
        rhodz_moist(i,j,k) = rho1(i,j) * deltaz(i,j,k)

      END IF  ! k  <   tdims%k_end

      ! Calculate total moisture
      q_total(i,j) = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)

      IF (l_mcr_qcf2) THEN
        q_total(i,j) = q_total(i,j) + qcf2(i,j,k)
      END IF  ! l_mcr_qcf2

      IF (l_mcr_qrain) THEN
        q_total(i,j) = q_total(i,j) + qrain(i,j,k)
      END IF  ! l_mcr_qrain

      IF (l_mcr_qgraup) THEN
        q_total(i,j) = q_total(i,j) + qgraup(i,j,k)
      END IF  ! l_mcr_qgraup


      ! Rho_r2 uses the moist density of air. If the mixing
      ! ratio framework is in place then we need to also know
      ! the dry density of air.
      IF (l_mr_physics) THEN
        rhodz_dry(i,j,k) = rhodz_moist(i,j,k)                                 &
             / (1.0 + q_total(i,j))
      ELSE
        rhodz_dry(i,j,k) = rhodz_moist(i,j,k)                                 &
             * (1.0 - q_total(i,j))
      END IF  ! l_mr_physics

    END DO  ! i
  END DO  ! j
END DO  ! k
!$OMP END PARALLEL DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE mphys_air_density
END MODULE mphys_air_density_mod
