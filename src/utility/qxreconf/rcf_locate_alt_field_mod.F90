! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Sets the data_source array corresponding to the fields array.

MODULE rcf_locate_alt_field_mod
!  Subroutine rcf_locate_alt_field - sets field to use alternative.
!
! Description:
!   When the field is not in the input dump sometimes we can fallback to using
!   an alternative.
!
! Method:
!   Using predefined list we can fallback to using other fields if the field
!   is not available in the input dump.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_LOCATE_ALT_FIELD_MOD'

CONTAINS
SUBROUTINE rcf_locate_alt_field( field, fields, fields_count, pos)

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_Field_Type_Mod, ONLY: &
    Field_type

USE Rcf_Locate_mod, ONLY: &
    Rcf_Locate

USE Submodel_Mod, ONLY:   &
    atmos_im

USE um_stashcode_mod

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE Rcf_Ppx_Info_Mod, ONLY: &
    STM_record_type

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (field_type)                         :: field
TYPE (field_type), POINTER                :: fields(:)
INTEGER                                   :: fields_count
INTEGER                                   :: pos
CHARACTER(LEN=errormessagelength)         :: cmessage

! Local variables
CHARACTER(LEN=*), PARAMETER :: RoutineName = 'RCF_LOCATE_ALT_FIELD'
INTEGER :: other_sec
INTEGER :: other_item
INTEGER :: errorstatus
TYPE(STM_record_type), POINTER :: other_stm
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

other_sec = -1
other_item = -1
NULLIFY(other_stm)

! Lets see what field we want
SELECT CASE( field % stashmaster % model )
CASE ( atmos_im )
  SELECT CASE( field % stashmaster % section )
  CASE ( stashcode_prog_sec )
    SELECT CASE( field % stashmaster % item )
    CASE ( stashcode_theta )
      ! We keep the STM as stashcode_t.
      other_sec  = stashcode_t/1000
      other_item = MOD(stashcode_t, 1000)
      IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
        WRITE(umMessage,*) "Using ", stashcode_t, " for ", field % stashmaster % NAME
        CALL umPrint(umMessage,src='rcf_locate_alt_field_mod')
      END IF

    END SELECT
  END SELECT
END SELECT

IF (other_sec == -1 .OR. other_item == -1) THEN
  WRITE(cmessage,*) "Alternative field not found for ", &
                    field % stashmaster % NAME
  errorstatus = 10
  CALL ereport(RoutineName, errorstatus, cmessage)
END IF

! Now lets find the field.
CALL rcf_locate(other_sec, other_item, fields, fields_count, pos)

! Reset the STASHmaster if defined.
IF (ASSOCIATED(other_stm)) THEN
  fields(pos) % stashmaster => other_stm
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_locate_alt_field
END MODULE rcf_locate_alt_field_mod
