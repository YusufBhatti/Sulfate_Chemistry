! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE distribute_flash_mod

! Purpose: Calculates if and where cloud-to-ground lightning flashes
! will occur based on inputs of flash potential (passed from previous
! timesteps) and current flash rate diagnosed in this timestep. If the
! flash potential exceeds 1.0 in any grid column, a cloud-to-ground
! strike is deemed to occur. Any remaining flash potential is then assigned
! to have the same distribution as the graupel and then passed out
! as a prognostic variable to the next timestep.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

USE atm_fields_bounds_mod,  ONLY: tdims
USE timestep_mod,           ONLY: timestep
USE electric_constants_mod, ONLY: i,j,k

! Dr Hook modules
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='DISTRIBUTE_FLASH_MOD'

CONTAINS

SUBROUTINE distribute_flash( qgraup, qgtot, flash, flash_pot1d, flash_pot,     &
                             num_flashes )

IMPLICIT NONE
!---------------------
! Subroutine arguments
!----------------------

REAL, INTENT(IN)     :: qgraup(      tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                1 : tdims%k_end )
! Graupel mass mixing ratio [ ]
REAL, INTENT(IN)     :: qgtot(       tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
! Vertical sum of qgraup over the column [ ]
REAL, INTENT(IN)     :: flash(       tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
! Lightning flash rate [s-1]
REAL, INTENT(IN)     :: flash_pot1d( tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
! 2D flash potential [ ]
REAL, INTENT(INOUT)  :: flash_pot  ( tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
! 3D flash potential mixing ratio [ ]

REAL, INTENT(OUT) :: num_flashes(    tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
! Number of lightning strikes in a given grid box [ ]


!----------------------
! Local variables
!----------------------

REAL :: flash_pot_1 ! Flash potential at a single point
REAL :: k_factor    ! Graupel-charge distribution factor

REAL, PARAMETER :: flashmin = 1.0e-10 ! Minimum flash potential carried

INTEGER :: num_strikes ! number of lightning strikes at a single point

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='DISTRIBUTE_FLASH'

!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

num_flashes(:,:) = 0.0

DO j = tdims%j_start, tdims%j_end

  DO i = tdims%i_start, tdims%i_end

    flash_pot_1 = flash_pot1d(i,j) + ( flash(i,j) * timestep )

    IF (flash_pot_1 >= 1.0) THEN
      ! Need to output a flash here

      num_strikes = FLOOR(flash_pot_1)

      num_flashes(i,j) = REAL(num_strikes)

      flash_pot_1 = flash_pot_1 - REAL(num_strikes)

    END IF ! flash_pot_1 > 1.0

    ! Second loop now needed as flash_pot_1 could have been
    ! above 1.0 but will now could be between 0.0 and 1.0

    IF (flash_pot_1 >= flashmin) THEN

      ! ignore very small flash amounts- this means flashes will not
      ! be conserved exactly, but we probably don't care about this
      ! mostly because there will be some leakage of electricity
      ! and it is not a closed system in the UM yet.

      k_factor = flash_pot_1 / qgtot(i,j)

      ! Distribute the flash in the vertical for passing on to next
      ! timestep

      DO k = 1, tdims%k_end

        flash_pot(i,j,k) = k_factor * qgraup(i,j,k)

      END DO ! k

    END IF ! flash_pot_1 > flashmin

  END DO ! i

END DO ! j

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE distribute_flash

END MODULE distribute_flash_mod
