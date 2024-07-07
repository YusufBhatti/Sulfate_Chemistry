! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stochastic Physics
MODULE stph_writeentry_mod

USE stochastic_physics_run_mod, ONLY: stphseed_unit

IMPLICIT NONE

INTERFACE stph_writeentry
MODULE PROCEDURE stph_writeentry_int, stph_writeentry_real
END INTERFACE stph_writeentry

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='STPH_WRITEENTRY_MOD'

CONTAINS

! -------------------------------------------------------------------
! Subroutine in which stpharray is passed as an 1D integer array
! -------------------------------------------------------------------
SUBROUTINE stph_writeentry_int( stpharray, sizearray)

! This routine is passed a integer array and its size as arguments
! and writes unformatted to the pre-defined stoch phys
! output file unit=stphseed_unit

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stochastic_physics_run_mod, ONLY: stphseed_unit

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN) :: sizearray            ! The size of the array
INTEGER, INTENT(IN) :: stpharray(sizearray) ! The random seed

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STPH_WRITEENTRY_INT'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Write the value to the file stream
WRITE(stphseed_unit) stpharray

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stph_writeentry_int

! -------------------------------------------------------------------
! Subroutine in which stpharray is passed as a 2D real array
! -------------------------------------------------------------------
SUBROUTINE stph_writeentry_real( stpharray, sizearray)

! This routine is passed a real array and its size as arguments
! and writes unformatted to the pre-defined stoch phys
! output file unit=stphseed_unit

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE stochastic_physics_run_mod, ONLY: stphseed_unit

IMPLICIT NONE

! Arguments
REAL,    INTENT(IN) :: stpharray(:,:) ! The random field coefficients
INTEGER, INTENT(IN) :: sizearray      ! The size of the array

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='STPH_WRITEENTRY_REAL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Write the value to the file stream
WRITE(stphseed_unit) stpharray

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE stph_writeentry_real

END MODULE stph_writeentry_mod
