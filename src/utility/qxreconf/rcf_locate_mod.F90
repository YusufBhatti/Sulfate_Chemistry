! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Find a field in a field list by stash item number

MODULE Rcf_Locate_mod

!  Subroutine Rcf_Locate - locates a field in a field array
!
! Description:
! Fortran 90 replacement for locate1, using structures rather
! than a whole bunch of arrays. Will return the position of the
! field in the array if it exists otherwise will die gracefully.
! Optional argument zero_ok allows return of 0 in certain circumstances
!
! Method:
!   Loop through the array until the correct field is found
!   in section 0.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_LOCATE_MOD'

CONTAINS

SUBROUTINE Rcf_Locate( section, item, fields, nfields, pos, zero_ok_arg )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Ereport_Mod, ONLY: &
    Ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER, INTENT(IN)                :: nfields
INTEGER, INTENT(IN)                :: section
INTEGER, INTENT(IN)                :: item
INTEGER, INTENT(OUT)               :: pos
TYPE (field_type), INTENT(IN)      :: fields( nfields )
LOGICAL, INTENT(IN), OPTIONAL      :: zero_ok_arg


! Local variables
INTEGER                            :: i
INTEGER                            :: ErrorStatus
LOGICAL                            :: zero_ok
CHARACTER (LEN=*), PARAMETER       :: RoutineName='RCF_LOCATE'
CHARACTER (LEN=errormessagelength)   :: Cmessage
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
zero_ok = .FALSE.
IF ( PRESENT( zero_ok_arg ) ) THEN
  zero_ok = zero_ok_arg
END IF

pos = 0

DO i = 1, nfields
  IF (fields(i) % stashmaster % item    == item .AND. &
      fields(i) % stashmaster % section == section ) THEN
    pos = i
    EXIT
  END IF
END DO

IF (pos == 0 .AND. .NOT. zero_ok) THEN
  ErrorStatus = 10
  WRITE(Cmessage, '(A, I3, I5)') 'Unable to Rcf_Locate field with STASHcode ', &
                                  section, item
  CALL Ereport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Locate
END MODULE Rcf_Locate_Mod
