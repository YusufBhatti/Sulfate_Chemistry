! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Initialises all variables required in the electric
!          scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

MODULE electric_init_mod

USE gen_phys_inputs_mod,ONLY: l_mr_physics
USE mphys_inputs_mod,       ONLY: l_mcr_qcf2
USE electric_constants_mod, ONLY: i, j, k
USE atm_fields_bounds_mod,  ONLY: tdims
USE mphys_ice_mod,          ONLY: qcfmin
USE mphys_bypass_mod,       ONLY: qcf2_idims_start, qcf2_idims_end,            &
                                  qcf2_jdims_start, qcf2_jdims_end,            &
                                  qcf2_kdims_start, qcf2_kdims_end
! Dr Hook modules
USE yomhook,                ONLY: lhook, dr_hook
USE parkind1,               ONLY: jprb, jpim

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='ELECTRIC_INIT_MOD'

CONTAINS

SUBROUTINE electric_init( storm_field, qcf, qcf2, qgraup, qcf_tot, rhodz_dry,  &
                          rhodz_moist, flash_pot, rhodz, flash, qgtot,         &
                          flash_pot1d, fr1_mc, fr2_mc )

IMPLICIT NONE

REAL, INTENT(IN)     :: qcf(         tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
REAL, INTENT(IN)     :: qcf2(     qcf2_idims_start : qcf2_idims_end,           &
                                  qcf2_jdims_start : qcf2_jdims_end,           &
                                  qcf2_kdims_start : qcf2_kdims_end )
REAL, INTENT(IN)     :: qgraup(      tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
REAL, INTENT(IN)     :: rhodz_dry(   tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
REAL, INTENT(IN)     :: rhodz_moist( tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
REAL, INTENT(IN)     :: flash_pot(   tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
REAL, INTENT(OUT)    :: rhodz(       tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
LOGICAL, INTENT(OUT) :: storm_field( tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
REAL, INTENT(OUT)    :: qcf_tot(     tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 1 : tdims%k_end )
REAL, INTENT(OUT)    :: flash(       tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
REAL, INTENT(OUT)    :: qgtot(       tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
REAL, INTENT(OUT)    :: flash_pot1d( tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
REAL, INTENT(OUT)    :: fr1_mc(      tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )
REAL, INTENT(OUT)    :: fr2_mc(      tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end )

! Local variables
REAL :: qg_sum ! component of qgtot
REAL :: fp_sum ! component of flash_pot

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='ELECTRIC_INIT'

!==================================================================
! Start the subroutine
!==================================================================

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise storm variables to zero or false

!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,j,k) SHARED(tdims,qgtot,flash_pot1d, &
!$OMP flash,fr1_mc,fr2_mc,storm_field,l_mcr_qcf2,qcf_tot,qcf,qcf2,qgraup,   &
!$OMP flash_pot,l_mr_physics,rhodz,rhodz_dry,rhodz_moist)

!$OMP DO SCHEDULE(STATIC)
DO j = tdims%j_start , tdims%j_end
  DO i= tdims%i_start , tdims%i_end
    qgtot(i,j)       = 0.0
    flash_pot1d(i,j) = 0.0
    flash(i,j)       = 0.0
    fr1_mc(i,j)      = 0.0
    fr2_mc(i,j)      = 0.0
    storm_field(i,j) = .FALSE.
  END DO
END DO
!$OMP END DO NOWAIT

IF (l_mcr_qcf2) THEN

   ! Two cloud ice categories in use. Total ice is the sum of
   ! the two categories

  DO k = 1, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        qcf_tot(i,j,k) = qcf(i,j,k) + qcf2(i,j,k)
        qgtot(i,j) = qgtot(i,j) + qgraup(i,j,k)

        ! Only include flash potential where graupel is present
        IF (qgraup(i,j,k) >= qcfmin) flash_pot1d(i,j) =                      &
                                     flash_pot1d(i,j) + flash_pot(i,j,k)

      END DO ! i
    END DO ! j
!$OMP END DO NOWAIT
  END DO ! k

ELSE   ! not l_mcr_qcf2

   ! Only one cloud ice category in use. All ice present is in one
   ! category only, so qcf_tot = qcf

  DO k = 1, tdims%k_end
!$OMP DO SCHEDULE(STATIC)
    DO j = tdims%j_start, tdims%j_end
      DO i = tdims%i_start, tdims%i_end

        qcf_tot(i,j,k) = qcf(i,j,k)
        qgtot(i,j) = qgtot(i,j) + qgraup(i,j,k)

        ! Only include flash potential where graupel is present
        IF (qgraup(i,j,k) >= qcfmin) flash_pot1d(i,j) =                      &
                                     flash_pot1d(i,j) + flash_pot(i,j,k)

      END DO ! i
    END DO ! j
!$OMP END DO NOWAIT
  END DO ! k

END IF ! l_mcr_qcf2

! Set rhodz based on whether using specific humidity
! or mixing ratio formulation

IF (l_mr_physics) THEN

!$OMP DO SCHEDULE(STATIC)
  DO k=1 , tdims%k_end
    DO j = tdims%j_start , tdims%j_end
      DO i= tdims%i_start , tdims%i_end
        rhodz(i,j,k) = rhodz_dry(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO NOWAIT

ELSE ! l_mr_physics

!$OMP DO SCHEDULE(STATIC)
  DO k=1 , tdims%k_end
    DO j = tdims%j_start , tdims%j_end
      DO i= tdims%i_start , tdims%i_end
        rhodz(i,j,k) = rhodz_moist(i,j,k)
      END DO ! i
    END DO ! j
  END DO ! k
!$OMP END DO NOWAIT

END IF ! l_mr_physics

!$OMP END PARALLEL

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE electric_init

END MODULE electric_init_mod
