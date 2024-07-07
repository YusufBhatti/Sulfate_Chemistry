! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE check_stm_codes_mod

IMPLICIT NONE

! Description:
! Checks that STASHmaster codes are consistent.

! Method:
! We perform three consistency checks, detailed below.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='CHECK_STM_CODES_MOD'

CONTAINS
SUBROUTINE check_stm_codes(model, sectn, item, grid, datat, dumpp)

USE ereport_mod, ONLY: ereport
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE errormessagelength_mod, ONLY: errormessagelength
IMPLICIT NONE

INTEGER, INTENT(IN) :: model
INTEGER, INTENT(IN) :: sectn
INTEGER, INTENT(IN) :: item
INTEGER, INTENT(IN) :: grid
INTEGER, INTENT(IN) :: datat
INTEGER, INTENT(IN) :: dumpp

! Local
INTEGER :: n1, n2, n3
INTEGER :: errorstatus
LOGICAL :: dumpp_land_packed
LOGICAL :: grid_land_only
CHARACTER(LEN=errormessagelength) :: cmessage

! Parameters
CHARACTER(LEN=*), PARAMETER :: RoutineName = "CHECK_STM_CODES"

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

n1 = MOD(dumpp,10)
n2 = MOD(dumpp/10,10)
n3 = MOD(dumpp/100,10)

dumpp_land_packed = n2 == 2 .AND. n3 == 1
grid_land_only = grid == 21

! Check grid type and dump packing is consistent.
IF ( dumpp_land_packed .neqv. grid_land_only ) THEN
  WRITE(cmessage,'(A,I3,A,I3,A,I4,A)') &
    "In STASHmaster entry (model = ", &
    model, ", section = ", sectn, ", item = ", item, &
    ") land packing option and grid definition are inconsistent"
  errorstatus = 1
  CALL ereport(routinename, errorstatus, cmessage)
END IF

! Check non-real datatype is not 32-bit packed in the dump.
IF ( datat /= 1 .AND. n1 == 2) THEN
  WRITE(cmessage,'(A,I3,A,I3,A,I4,A)') &
    "In STASHmaster entry (model = ", &
    model, ", section = ", sectn, ", item = ", item, &
    ") 32-bit packing is inconsistent with non-real datatypes"
  errorstatus = 2
  CALL ereport(routinename, errorstatus, cmessage)
END IF

! Check that data type is valid no: +/-1 to +/-3
IF (( ABS(datat)  <  1) .OR.                &
    ( ABS(datat)  >  3)) THEN
  WRITE(cmessage,'(A,I3,A,I3,A,I4,A,I3)') &
    "In STASHmaster entry (model = ", &
    model, ", section = ", sectn, ", item = ", item, &
    ") Invalid data type:", datat
  errorstatus = 3
  CALL ereport(routinename, errorstatus, cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
END SUBROUTINE check_stm_codes

END MODULE check_stm_codes_mod
