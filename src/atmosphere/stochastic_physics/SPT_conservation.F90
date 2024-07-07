! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Description: This routine calculates the vertical integrated q
!     before and after SPT, computes a ratio to modulate the
!     SPT increments (delta_q) so they conserve q in the column.
!
!     Code Owner: Please refer to the UM file CodeOwners.txt
!     This file belongs in section: Stochastic Physics

MODULE SPT_conservation_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SPT_CONSERVATION_MOD'

CONTAINS

SUBROUTINE SPT_conservation(delta_q, delta_theta,                       &
                            q,exner_theta_levels, rho_in)

! Bounds for arrays
USE atm_fields_bounds_mod, ONLY:                                  &
    tdims, tdims_s,                                               &
    pdims, pdims_s

! Load atm_step to save numbers
USE atm_step_local,      ONLY: first_atmstep_call

! SPT UMUI settings passed in via NameList READ
USE stochastic_physics_run_mod, ONLY:                             &
    l_spt_qcons,l_spt_mse,spt_bot_tap_lev,spt_top_tap_lev

! Model level-height modules
USE level_heights_mod,     ONLY:                                  &
    r_theta_levels, r_rho_levels

! Factors to compute the Moist Static Energy
USE planet_constants_mod, ONLY: cp
USE water_constants_mod, ONLY: lc

!Use unscaled dry rho or r2 scaled wet rho
USE gen_phys_inputs_mod, ONLY: l_mr_physics

! DrHook modules
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER :: i,j,k
    ! For loop indexes
     
REAL, SAVE :: LdivCp
   ! Save L/Cp factor to compute MSE conserving increments      
     
REAL ::                                                                &
  
      delta_q(tdims%i_start:tdims%i_end,                               &
              tdims%j_start:tdims%j_end,                               &
              1:tdims%k_end)                                           &
    ! SPT increments for Q               
  ,   delta_theta(tdims%i_start:tdims%i_end,                           &
                  tdims%j_start:tdims%j_end,                           &
                  1:tdims%k_end)                                       &
    ! Matix holding the increments for theta
  ,   q  (tdims%i_start:tdims%i_end,                                   &
          tdims%j_start:tdims%j_end,                                   &
          tdims%k_start:tdims%k_end)                                   &
    ! Q value before SPT perturbations
  ,   exner_theta_levels(tdims_s%i_start:tdims_s%i_end,                &
                       tdims_s%j_start:tdims_s%j_end,                  &
                       tdims_s%k_start:tdims_s%k_end)                  &
    ! Exner pressure on theta levels (to convert T -> theta)
  ,   rho_in(pdims_s%i_start:pdims_s%i_end,                            &
             pdims_s%j_start:pdims_s%j_end,                            &
             pdims_s%k_start:pdims_s%k_end)
    ! density (depends on the l_mr_physics definition)
REAL ::                                                                &
        rho  (pdims%i_start:pdims%i_end,                               &
              pdims%j_start:pdims%j_end,                               &
              pdims%k_start:pdims%k_end)                               &
          ! Final unscaled Rho (density)
  ,     Q_bef(tdims%i_start:tdims%i_end,                               &
              tdims%j_start:tdims%j_end)                               &
          ! Q Before
  ,     Q_aft(tdims%i_start:tdims%i_end,                               &
              tdims%j_start:tdims%j_end)                               &
          ! Q after
  ,     alpha(tdims%i_start:tdims%i_end,                               &
              tdims%j_start:tdims%j_end)
          ! alpha: Q_bef/ Q)a

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='SPT_CONSERVATION'

! ------------------------------------------------------------------
!      END OF VARIABLE DECLARATIONS - START OF THE CODE
! ------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! +++++++++++++++++++ Apply conservation to the water column

IF (l_spt_qcons) THEN
    
  IF (l_mr_physics) THEN
  ! If q is mixing ratio, then rho is unscaled by the r**2 factor
    DO k = pdims%k_start, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          rho(i,j,k) = rho_in(i,j,k)
        END DO
      END DO
    END DO
  ELSE
  ! If q is specific humidity, then rho is scaled by the r**2 factor
    DO k = pdims%k_start, pdims%k_end
      DO j = pdims%j_start, pdims%j_end
        DO i = pdims%i_start, pdims%i_end
          rho(i,j,k) = rho_in(i,j,k)/(r_rho_levels(i,j,k) *            &
                       r_rho_levels(i,j,k))
        END DO
      END DO
    END DO
  END IF 

   ! Initialize Q_bef and Q_aft
  Q_bef=0.0
  Q_aft=0.0

  ! Compute total vertical tendency before SPT forcing
  ! Q_bef = integral ( q x rho x delta_Z )
  DO k= spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        Q_bef(i,j)= Q_bef(i,j) + 0.5*( rho(i,j,k) +  rho(i,j,k+1) )    &
                                     * q(i,j,k) *                      &
                    ( r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k) )

      END DO
    END DO
  END DO

  ! Compute total vertical tendency After SPT forcing
  DO k= spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        Q_aft(i,j)= Q_aft(i,j) + 0.5*( rho(i,j,k) +  rho(i,j,k+1) ) *  &
                        ( q(i,j,k) + delta_q(i,j,k) )  *               &
                    ( r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k) )
      END DO
    END DO
  END DO

  ! Get ratio between both vertical estimates
  DO j = tdims%j_start, tdims%j_end
    DO i = tdims%i_start, tdims%i_end
      alpha(i,j)= Q_bef(i,j) / Q_aft(i,j)
    END DO
  END DO

  ! Rescale perturbations
  DO k= spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end
        delta_q(i,j,k) =  delta_q(i,j,k) * alpha(i,j)+ q(i,j,k) *     &
                                          ( alpha(i,j)-1.0 )
      END DO
    END DO
  END DO

END IF

! +++++++++++++++++++ Apply conservation to the MSE
IF (l_spt_mse) THEN

  IF (first_atmstep_call) LdivCp= lc/cp

  DO k= spt_bot_tap_lev, spt_top_tap_lev
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        delta_theta(i,j,k) = (-1.0)* LdivCp* delta_q(i,j,k)              &
                             / exner_theta_levels(i,j,k)
      END DO
    END DO
  END DO

ELSE

  IF (first_atmstep_call) LdivCp=1.0

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE SPT_conservation

END MODULE SPT_conservation_mod
