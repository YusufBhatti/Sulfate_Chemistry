! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE file_mod

USE filenamelength_mod, ONLY: filenamelength
USE missing_data_mod,   ONLY: imdi
! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook
IMPLICIT NONE

! Description:
!   A base class for any file. 
!     Subclasses of this class represent different types of files, e.g.
!     datafile
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'FILE_MOD'

TYPE, PUBLIC :: file_type
  PRIVATE
  CHARACTER(LEN=filenamelength), PUBLIC :: filename
  
  ! unit number needs to be public so subclasses can access it.
  INTEGER, PUBLIC :: unit_num = imdi

  CONTAINS

  PROCEDURE, PASS :: open_file
  PROCEDURE, PASS :: close_file

END TYPE file_type

CONTAINS

! Simple open/close methods. Subclasses may well override these to do
! more sophisticated things.

SUBROUTINE open_file(this, filename)
USE file_manager, ONLY: assign_file_unit
IMPLICIT NONE
CLASS(file_type), INTENT(INOUT) :: this
CHARACTER(LEN=filenamelength), INTENT(IN) :: filename
CHARACTER(LEN=*), PARAMETER :: routinename = 'OPEN_FILE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
this%filename = filename
CALL assign_file_unit(this%filename, this%unit_num, handler="fortran")
OPEN(UNIT=this%unit_num, FILE=this%filename)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE open_file


SUBROUTINE close_file(this)
USE file_manager, ONLY: release_file_unit
IMPLICIT NONE
CLASS(file_type), INTENT(INOUT) :: this
CHARACTER(LEN=*), PARAMETER :: routinename = 'CLOSE_FILE'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
CLOSE(this%unit_num)
CALL release_file_unit(this%unit_num, handler="fortran")

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE close_file

END MODULE file_mod
