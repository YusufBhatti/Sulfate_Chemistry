! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE fr_mccaul_mod

! Purpose: Calculates flash rate based on McCaul et al (2009), Weather and
! forecasting.

! Full paper reference:
! McCaul, E. W., Goodman, S. J., LaCasse, K. M., Cecil, D. J. 2009.
! Forecasting Lightning Threat Using Cloud-Resolving Model Simulations.
! Weather and Forecasting, Volume 24, pp 709–729.
! doi: http://dx.doi.org/10.1175/2008WAF2222152.1

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

USE atm_fields_bounds_mod,  ONLY: tdims, wdims_s
USE electric_constants_mod, ONLY: i, j, k, minus15, mccaul_r1, mccaul_r2,      &
                                  min2sec
USE electric_inputs_mod,    ONLY: k1, k2

! Dr Hook modules
USE yomhook,               ONLY: lhook, dr_hook
USE parkind1,              ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='FR_MCCAUL_MOD'

CONTAINS

SUBROUTINE fr_mccaul(storm_field, t, w, qgraup, gwp, tiwp, flash, fr1_mc,      &
                     fr2_mc )

IMPLICIT NONE

LOGICAL, INTENT(IN) :: storm_field( tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end )
! Defines where storms exist in model

REAL, INTENT(IN) :: w(              wdims_s%i_start : wdims_s%i_end,           &
                                    wdims_s%j_start : wdims_s%j_end,           &
                                    wdims_s%k_start : wdims_s%k_end )
! Vertical velocity (m/s)

REAL, INTENT(IN) :: t(              tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end,               &
                                                1 : tdims%k_end )
! Temperature (K)
REAL, INTENT(IN) :: qgraup(         tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end,               &
                                                1 : tdims%k_end )
! Graupel mixing ratio (kg/kg)

REAL, INTENT(IN) :: gwp(            tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end)
REAL, INTENT(IN) :: tiwp(           tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end)
! Graupel and total ice water paths (kg m-2)

REAL, INTENT(INOUT) :: flash( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end )

REAL, INTENT(INOUT) :: fr1_mc(  tdims%i_start : tdims%i_end,                   &
                              tdims%j_start : tdims%j_end )

REAL, INTENT(INOUT) :: fr2_mc(  tdims%i_start : tdims%i_end,                   &
                              tdims%j_start : tdims%j_end )

! Flash rates (flash 1 and flash2 for diagnostics). Units: s-1


! local variables

INTEGER :: m15l, upper_limit, lower_limit

REAL :: flash1
REAL :: flash2

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='FR_MCCAUL'

!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Define upper and lower limits for minus 15 Celius calculation
upper_limit = tdims%k_end   - 2
lower_limit = tdims%k_start + 1

DO j = tdims%j_start, tdims%j_end

  DO i = tdims%i_start, tdims%i_end

    ! Initialise flash1 and flash2

    flash1 = 0.0
    flash2 = 0.0

    ! Only perform calculation if this point is defined as a storm point
    IF (storm_field(i,j)) THEN

      ! First find the minus 15 level
      DO k = tdims%k_end, 1, -1

        IF (t(i,j,k) > minus15) THEN

          m15l = k+1 ! define minus 15 level
          EXIT ! leave the loop early to avoid overwriting

        END IF

      END DO ! k

      ! Determine flash rates

      IF (m15l > upper_limit .OR. m15l < lower_limit) THEN
        ! Set flash1 to be zero if you cannot define the
        ! -15C level in the model.

        flash1 = 0.0

      ELSE ! m15l

        flash1 = k1 * w(i,j, m15l) * qgraup(i,j, m15l)

      END IF ! m15l

      flash2 = k2 * ( gwp(i,j) + tiwp(i,j) )

      flash(i,j) = min2sec * ( (mccaul_r1 * flash1) + (mccaul_r2 * flash2) )

      ! Determine the value of each component for diagnostics
      ! (these will have been initialised to zero in electric_init).

      fr1_mc(i,j) = min2sec * flash1
      fr2_mc(i,j) = min2sec * flash2

    END IF ! storm_field

  END DO ! i

END DO ! j

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE fr_mccaul

END MODULE fr_mccaul_mod
