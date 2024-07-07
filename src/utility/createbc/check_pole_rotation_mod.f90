! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
MODULE check_pole_rotation_mod

! Description:
!   A routine to compare to determine if the rotation of the poles
!   is the same for two grids.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CreateBC
!
! Language: Fortran 2003.
!
! This code is written to UMDP3 standards.

! DrHook Modules
USE parkind1, ONLY: jprb, jpim
USE yomhook, ONLY: lhook, dr_hook

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName = 'CHECK_POLE_ROTATION_MOD'

CONTAINS

SUBROUTINE check_pole_rotation(source_pole_lat, source_pole_long, &
                               target_pole_lat, target_pole_long, &
                               l_same_rotation, l_source_rotated, &
                               l_target_rotated)

USE umPrintMgr, ONLY:  umPrint, umMessage

IMPLICIT NONE

REAL, INTENT(IN) :: source_pole_lat
REAL, INTENT(IN) :: source_pole_long
REAL, INTENT(IN) :: target_pole_lat
REAL, INTENT(IN) :: target_pole_long
LOGICAL, INTENT(INOUT) :: l_same_rotation
LOGICAL, INTENT(INOUT) :: l_source_rotated
LOGICAL, INTENT(INOUT) :: l_target_rotated
CHARACTER(LEN=*), PARAMETER :: routinename = 'CHECK_POLE_ROTATION'
! DrHook Variables
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF ((source_pole_lat == target_pole_lat) .AND. &
     (source_pole_long == target_pole_long)) THEN
  l_same_rotation = .TRUE.
  WRITE(umMessage, '(A,I5)') '[INFO] Input file grid and lbc output grid '// &
       'have the same rotation.'
  CALL umPrint(umMessage, src='createbc')
ELSE
  l_same_rotation = .FALSE.
  WRITE(umMessage, '(A,I5)') '[INFO] Input file grid and lbc output grid '// &
       'have a different rotation.'
  CALL umPrint(umMessage, src='createbc')
END IF

! The logical flag in the input file has already been assigned during the
! reading of the file header.  Here we just print out the rotation.
IF ( .NOT. l_source_rotated ) THEN
  WRITE(umMessage, '(A,I5)') '[INFO] Input file grid has a standard pole'
  CALL umPrint(umMessage, src='createbc')
ELSE
  WRITE(umMessage, '(A,I5)') '[INFO] Input file grid has a rotated pole'
  CALL umPrint(umMessage, src='createbc')
END IF

IF ( target_pole_lat == 90.0 .AND. target_pole_long == 0.0 ) THEN
  WRITE(umMessage, '(A,I5)') '[INFO] Output file grid has a standard pole'
  CALL umPrint(umMessage, src='createbc')
  l_target_rotated = .FALSE.
ELSE
  l_target_rotated = .TRUE.
  WRITE(umMessage, '(A,I5)') '[INFO] Output file grid has a rotated pole'
  CALL umPrint(umMessage, src='createbc')
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_pole_rotation

END MODULE check_pole_rotation_mod
