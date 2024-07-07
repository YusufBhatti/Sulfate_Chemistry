! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routines to read requested fields from a UM ancillary fieldsfile.

MODULE get_anc_flds_mod

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='GET_ANC_FLDS_MOD'

CONTAINS

SUBROUTINE get_anc_flds ( NumFlds,         &  ! in
                          MaxFlds,         &  ! in
                          STCode,          &  ! in
                          Store,           &  ! in
                          PPHdrMod,        &  ! in
                          UMHdr,           &  ! in
                          bzy,bzx,bdy,bdx, &  ! inout
                          Fields,          &  ! inout
                          ErrorStatus)        ! inout


USE IO_Mod, ONLY:         &
  PP_Header_type,         &
  PP_Field_type,          &
  UM_Header_type,         &
  LSMField
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY: ereport

USE umPrintMgr                              ! for writing output
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

USE errormessagelength_mod, ONLY: errormessagelength

USE Packing_codes_mod, ONLY: PC_BitMask_CompressType

IMPLICIT NONE
! ------------------------------------------------------------------------------
! Description:
!   This subroutine reads a given number of fields, matching given
!   criteria from an open UM fieldsfile.  The fields are then store in
!   given locations in the 'Fields' array.
!
! Method:
!   See online documentation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section:  Utility - crmstyle_coarse_grid
!
! Code Description:
!   Language:           Fortran 90
!   This code is written to UM programming standards version 8.3
! ------------------------------------------------------------------------------
! Subroutine arguments
!-------------------------------------------------------------------------------

! Subroutine arguments:
INTEGER, INTENT(IN) :: NumFlds
INTEGER, INTENT(IN) :: MaxFlds
INTEGER, INTENT(IN) :: STCode(NumFlds)
INTEGER, INTENT(IN) :: Store(NumFlds)

LOGICAL, INTENT(IN) :: PPHdrMod

TYPE(UM_Header_type), INTENT(IN) :: UMHdr

REAL, INTENT(INOUT) :: &
  bzy                  & ! pp grid details
 ,bzx                  &
 ,bdy                  &
 ,bdx

TYPE(PP_Field_type), INTENT(INOUT) :: Fields(MaxFlds)

INTEGER, INTENT(INOUT) :: ErrorStatus


!------------------------------------------------------------------------------
! Local Variables:
!------------------------------------------------------------------------------

INTEGER :: i                      ! local counter
INTEGER :: ifld
INTEGER :: ihdr
INTEGER :: NumCols           ! Number of columns
INTEGER :: NumRows           ! Number of rows
INTEGER :: DFldCount
INTEGER :: WFldCount
INTEGER :: PFldCount
INTEGER :: NumStore

LOGICAL :: LDecode = .TRUE.
LOGICAL :: Field_Not_Found

REAL, POINTER :: TempArray(:,:)

CHARACTER(LEN=*), PARAMETER :: RoutineName = "GET_ANC_FLDS"
CHARACTER(LEN=errormessagelength) :: ErrMessage

! Required by Dr hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!------------------------------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

IF ( .NOT. ASSOCIATED(UMHdr % Lookup) ) THEN
  WRITE(umMessage,'(A)') "Cannot get field since file uninitialised - skipping"
  CALL umPrint(umMessage,src=RoutineName)
  errorstatus = statuswarning
  IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
  RETURN
END IF

! Set defaults
Field_Not_Found = .FALSE.
DfldCount = 0
WFldCount = 0
PFldCount = 0
NumStore  = 0

!-------------------------------
! Loop through required fields
!-------------------------------
DO i = 1,NumFlds

  ifld = Store(i)     ! Array index where field should be stored

  ! If Store is less than 1 then assume we can just skip
  IF (ifld < 1) THEN
    CYCLE
  END IF
  ! Search to make sure Store is not set
  IF (COUNT(Store(i) == Store(1:i)) > 1) THEN
    DFldCount = DFldCount + 1
    WRITE(umMessage,'(A,I7,2A)') "Field ", ifld, " is already found - ",  &
                         "storing first match"
    CALL umPrint(umMessage,src=RoutineName)
    CYCLE
  END IF

  NumStore = NumStore + 1

  ! Clear space for field
  IF ( ASSOCIATED( Fields(ifld) % RData ) ) THEN
    DEALLOCATE( Fields(ifld) % RData )
    NULLIFY( Fields(ifld) % RData )
  END IF

  !----------------------------------------
  !  Search lookup for the required FIELD
  !----------------------------------------
  ! Search on LBTYP, LBLEV, LBFT - similar check in copyflds
  Fields(ifld) % LookupPos = 0


  DO ihdr = 1,UMHdr % NumFlds
    ! Need to use ABS since a negative is used by getb_flds to signify
    ! field has been interpolated to B grid.
    IF ( (ABS(STCode(i)) == UMHdr % Lookup(ihdr) % STCode) ) THEN

      Fields(ifld) % LookupPos = ihdr
      WFldCount = WFldCount + 1
      ! Take first field it finds (LBProc is allowed to be -1 here...)
      EXIT

    END IF
  END DO

  IF ( Fields(ifld) % LookupPos == 0 ) THEN
    WRITE( ErrMessage, '(A55,4I5)' )                                &
        "Field not found (STCode):                             ",   &
        STCode(i)
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, ErrMessage )
    ! Reset Errorstatus but return warning back to calling routine at end.
    ErrorStatus = StatusOk
    Field_Not_Found=.TRUE.
    CYCLE
  END IF

  Fields(ifld) % Hdr = UMHdr % Lookup( Fields(ifld) % LookupPos )

  ! want grid info from header for fields

  bzy = Fields(ifld) % Hdr % ZerothLat
  bdy = Fields(ifld) % Hdr % LatInt
  bzx = Fields(ifld) % Hdr % ZerothLon
  bdx = Fields(ifld) % Hdr % LonInt

  ! Allocate space for unpacked field
  ! Need to check if land packed.
  IF (MOD(Fields(ifld) % Hdr % LBPack / 10, 10) == PC_BitMask_CompressType) THEN
    IF (ASSOCIATED(LSMField % RData)) THEN
      Fields(ifld) % Hdr % NumCols = LSMField % Hdr % NumCols
      Fields(ifld) % Hdr % NumRows = LSMField % Hdr % NumRows
    ELSE
      ErrorStatus = StatusWarning

      CALL EReport( RoutineName, ErrorStatus, &
                  "Cannot store land/sea packed field without LSM - skipping.")
      ! Reset Errorstatus but return warning back to calling routine at end.
      ErrorStatus = StatusOk
      Field_Not_Found = .TRUE.
      CYCLE
    END IF
  END IF

  NumCols = Fields(ifld) % Hdr % NumCols
  NumRows = Fields(ifld) % Hdr % NumRows

  ALLOCATE ( Fields(ifld) % RData( NumCols*NumRows,1 ) )

  !----------------------------------------
  ! End of search.  Start retrieving data
  !----------------------------------------

  ! DEPENDS ON: readfld
  CALL ReadFld( UMHdr, LDecode, PPHdrMod, Fields(ifld),  &
                ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END IF

  ! Now we need to reshape the read in data to be the correct shape of the
  ! unpacked field.
  TempArray => Fields(ifld) % RData
  NULLIFY(Fields(ifld) % RData)
  ALLOCATE( Fields(ifld) % RData( NumCols,NumRows ) )
  Fields(ifld) % RData = RESHAPE(source = TempArray,                   &
                                 SHAPE  = (/NumCols,NumRows/))
  DEALLOCATE(TempArray)
  NULLIFY(TempArray)

END DO

WRITE(umMessage,'(2(I7,A))') WFldCount, " stored, out of ",     &
                              NumStore, " requested."
CALL umPrint(umMessage,src=RoutineName)
IF ( DFldCount /= 0) THEN
  WRITE(umMessage,'(A,I7,A)') " Warning: ", DFldCount," duplicate fields found."
  CALL umPrint(umMessage,src=RoutineName)
END IF
IF ( PFldCount /= 0) THEN
  WRITE(umMessage,'(A,I7,A)') " Info: ", PFldCount, " previous fields re-used."
  CALL umPrint(umMessage,src=RoutineName)
END IF

IF (Field_Not_Found .AND. ErrorStatus == StatusOK) THEN
  ! Reset errorstatus to pass warning back that a store didnt succeed.
  ErrorStatus = StatusWarning
END IF

!------------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-------------------------------------------------------------------------------
RETURN
END SUBROUTINE get_anc_flds

END MODULE get_anc_flds_mod
