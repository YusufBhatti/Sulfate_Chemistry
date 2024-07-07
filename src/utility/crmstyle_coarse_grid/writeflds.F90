! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routines to write fields to a UM fieldsfile

SUBROUTINE WriteFlds ( NumFlds,     &  ! in
                       MaxFlds,     &  ! in
                       MaxFldsOut,  &  ! in
                       PackType,    &  ! in
                       Source,      &  ! in
                       PackAcc,     &  ! in
                       Fields,      &  ! inout
                       UMHdr,       &  ! inout
                       ErrorStatus )   ! inout

! Description:
!
! Method:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: CRMstyle coarse grid
!
! Code Description:
!   Language:           Fortran 90
!   Software Standards: UMDP3 v6

USE IO_Mod, ONLY:         &
  PP_Field_type,          &
  PP_Header_type,         &
  UM_Header_type
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning,          &
  StatusFatal
USE ereport_mod, ONLY: ereport
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage
IMPLICIT NONE


! Subroutine Arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: MaxFldsOut
CHARACTER(LEN=*), INTENT(IN) :: PackType
INTEGER, INTENT(IN) :: Source(NumFlds)
REAL,    INTENT(IN) :: PackAcc(NumFlds)

TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)
TYPE(UM_Header_type), INTENT(INOUT) :: UMHdr
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "WriteFlds"

! Local Variables:
INTEGER :: i
INTEGER :: ifld
LOGICAL :: compressed = .FALSE.
TYPE(PP_Field_type) :: TempField

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

IF ( .NOT. ASSOCIATED(UMHdr % Lookup) ) THEN
  Errorstatus = Statusfatal
  CALL ereport ( Routinename, Errorstatus, &
                'Output file not been initialised.')
END IF

NULLIFY( TempField % RData )

! Check there is space in the lookup table
IF ( NumFlds + UMHdr % NumFlds > UMHdr % Len2Lookup ) THEN
  CALL umPrint( "Insufficient space allocated for Lookup Table", &
      src='writeflds')
  CALL umPrint( "Increase the value of MaxFldsOut in main prog to correct", &
      src='writeflds')
  ErrorStatus = StatusFatal

  CALL EReport( RoutineName, ErrorStatus, &
               "Insufficient space allocated for Lookup Table." )
END IF

DO i = 1,NumFlds
  ifld = Source(i)
  compressed = .FALSE.

  ! Packing if required
  IF ( TRIM(PackType) /= "NONE" ) THEN
    ! DEPENDS ON: pack_single
    CALL Pack_Single( PackAcc(i), PackType, Fields(ifld), TempField, &
                      ErrorStatus )
    IF ( ErrorStatus == StatusOK ) THEN
      compressed = .TRUE.
    ELSE IF ( ErrorStatus == StatusWarning ) THEN
      ! Since only warning report warning but we should be able to continue
      ! with unpacked.  Reset errorstatus.
      CALL umPrint( "Warning in Pack_Single - writing unpacked field.", &
          src='writeflds')
      ErrorStatus = StatusOK
    ELSE
      ! Report error message but do not reset errorstatus.
      CALL umPrint( "Error in Pack_Single - attempt to write unpacked field.", &
          src='writeflds')
    END IF
  END IF

  IF (.NOT. compressed) THEN
    ALLOCATE(TempField % &
        Rdata(Fields(ifld) % Hdr % NumRows * Fields(ifld) % Hdr % NumCols, 1))

    TempField % Rdata = RESHAPE(Fields(ifld) % Rdata, &
        [Fields(ifld) % Hdr % NumRows * Fields(ifld) % Hdr % NumCols, 1])
    TempField % Hdr = Fields(ifld) % Hdr
    TempField % LookupPos = Fields(ifld) % LookupPos
    TempField % ArrayPos = Fields(ifld) % ArrayPos
    TempField % Hdr % lblrec = &
            TempField % Hdr % NumRows * TempField % Hdr % NumCols
  END IF

  ! DEPENDS ON: fldout
  CALL FldOut( TempField, UMHdr, ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, "Error writing field." )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF

  DEALLOCATE( TempField % RData )  ! deallocate compressed field
  NULLIFY( TempField % RData )

END DO

9999 CONTINUE

END SUBROUTINE WriteFlds

