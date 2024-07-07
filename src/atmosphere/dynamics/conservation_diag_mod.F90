! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

MODULE conservation_diag_mod

USE um_parcore, ONLY: mype
USE filenamelength_mod, ONLY: filenamelength
USE file_manager, ONLY: assign_file_unit, release_file_unit
USE parkind1, ONLY: jpim, jprb       !DrHook
USE yomhook,  ONLY: lhook, dr_hook   !DrHook
USE timestep_mod,   ONLY: timestep_number

IMPLICIT NONE

LOGICAL, PARAMETER :: initial_timestep = .TRUE.
LOGICAL, PARAMETER :: not_initial_timestep = .FALSE.

PUBLIC :: print_conservation_diag, initial_timestep, not_initial_timestep
PRIVATE

CHARACTER(LEN=filenamelength), PARAMETER :: consv_file="conservation_diag.dat"

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CONSERVATION_DIAG_MOD'

CONTAINS

SUBROUTINE print_conservation_diag(total_rho, total_aam, total_ke, initial)

IMPLICIT NONE
!
! Description:
!
!
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Dynamics
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

REAL :: total_rho, total_aam, total_ke
LOGICAL :: initial
INTEGER :: istat
INTEGER :: eg_unit

REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='PRINT_CONSERVATION_DIAG'
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in, zhook_handle)

IF ( mype == 0 ) THEN
   IF (initial) THEN
      CALL assign_file_unit(consv_file, eg_unit, handler="fortran")
      OPEN(UNIT=eg_unit,FILE=consv_file)
      WRITE(UNIT=eg_unit,FMT='(I7,3E25.16)') 0, total_rho, total_aam, total_ke
   ELSE
      CALL assign_file_unit(consv_file, eg_unit, handler="fortran")
      OPEN(UNIT=eg_unit,FILE=consv_file,POSITION='append')
      WRITE(UNIT=eg_unit,FMT='(I9,3E25.16)') timestep_number,                  &
                                             total_rho, total_aam, total_ke
   END IF
   CLOSE(UNIT=eg_unit)
   CALL release_file_unit(eg_unit, handler="fortran")
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out, zhook_handle)
END SUBROUTINE
END MODULE
