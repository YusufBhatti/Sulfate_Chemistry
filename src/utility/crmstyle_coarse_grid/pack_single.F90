! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routines for packing and unpacking

SUBROUTINE Pack_Single( PackAcc,     &  ! in
                        PackType,    &  ! in
                        XpndField,   &  ! in
                        CompField,   &  ! inout
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
  LenWord,                &
  PP_Field_type,          &
  PP_Header_type,         &
  LSMField
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY: ereport
USE mask_compression, ONLY: compress_to_mask
USE wgdos_packing_mod, ONLY: wgdos_compress_field
USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage

USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY: PC_No_Packing, PC_WGDOS_Packing,                  &
  PC_LandMask_Compression, PC_BitMask_CompressType

IMPLICIT NONE


! Subroutine Arguments:
REAL, INTENT(IN)    :: PackAcc
CHARACTER(LEN=*), INTENT(IN) :: PackType
TYPE(PP_Field_type), INTENT(IN) :: XpndField

TYPE(PP_Field_type), INTENT(INOUT) :: CompField
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "Pack_Single"

! Local variables:
INTEGER, ALLOCATABLE :: int_data(:)
INTEGER :: fake_stash_code
INTEGER :: NumWords
INTEGER :: Num32BitWds
INTEGER :: NumLandPts

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

fake_stash_code = 0

! Can only pack REAL data
IF (XpndField % Hdr % LBUser1 /= 1) THEN
  WRITE(umMessage,'(A,I7)') "Cannot pack non-REAL data in STCode ", &
                    XpndField % Hdr % STCode
  CALL umPrint(umMessage,src='pack_single')
  errorstatus = StatusWarning
  GO TO 9999
END IF


CompField % Hdr = XpndField % Hdr
IF ( ASSOCIATED( Compfield % RData ) ) THEN
  DEALLOCATE( CompField % RData )
  NULLIFY( CompField % RData )
END IF


IF ( TRIM(PackType) == "WGDOS" ) THEN

  CompField % Hdr % bacc = PackAcc
  ALLOCATE( CompField % RData( CompField % Hdr % lblrec, 1 ) )
  ALLOCATE( int_data ( CompField % Hdr % lblrec ) )
  CALL wgdos_compress_field( XpndField % RData,            &
                             XpndField % Hdr % lblrec,     &
                             int_data,                     &
                             CompField % Hdr % lblrec,     &
                             XpndField % Hdr % NumCols,    &
                             Num32BitWds,                  &
                             INT(PackAcc),                 &
                             XpndField % Hdr % bmdi,       &
                             fake_stash_code, ErrorStatus )
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus,                &
                 'Problem reported by wgdos_compress_field' )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF
  CompField % RData(:,1) = TRANSFER(int_data, CompField % RData)
  NumWords = ( Num32BitWds-1 + (LenWord/32) ) * 32/LenWord
  CompField % Hdr % lblrec = NumWords
  CompField % Hdr % lbpack = XpndField % Hdr % lbpack + PC_WGDOS_Packing
ELSE IF ( TRIM(PackType) == "CRAY32" ) THEN
  CALL umPrint( "CRAY32 Packing not supported - leaving unpacked", &
      src='pack_single')
  ErrorStatus = StatusWarning
ELSE IF ( TRIM(PackType) == "GRIB" ) THEN
  CALL umPrint( "GRIB Packing not supported - leaving unpacked", &
      src='pack_single')
  ErrorStatus = StatusWarning
ELSE IF ( TRIM(PackType) == "RUNLEN" ) THEN
  CALL umPrint( "RUNLEN Packing not supported - leaving unpacked", &
      src='pack_single')
  ErrorStatus = StatusWarning
ELSE IF ( PackType == "LAND" ) THEN
  IF (.NOT. ASSOCIATED(LSMField % Rdata)) THEN
    ErrorStatus = StatusWarning
    CALL EReport( RoutineName, ErrorStatus, &
                  "LSM required to unpack land packed fields." )
    ErrorStatus = StatusWarning
  ELSE IF (CompField % Hdr % LBPack /= PC_No_Packing) THEN
    ErrorStatus = StatusWarning
    CALL EReport( RoutineName, ErrorStatus, &
                  "We can land-pack with unpacked fields" )
    ErrorStatus = StatusWarning
  ELSE
    ALLOCATE( CompField % RData( CompField % Hdr % lblrec, 1 ) )
    CALL compress_to_mask(XpndField % Rdata, CompField % Rdata, &
                          LSMField % Rdata /= 0.0, CompField % Hdr % lblrec, &
                          NumLandPts)
    CompField % Hdr % lblrec = NumLandPts
    CompField % Hdr % LBPack = CompField % Hdr % LBPack +                      &
                               100 * PC_LandMask_Compression +                 &
                               10 * PC_BitMask_CompressType + PC_No_Packing
  END IF
ELSE
  CALL umPrint( "PackType not recognised - leaving unpacked", &
      src='pack_single')
  ErrorStatus = StatusWarning
END IF

9999 CONTINUE

END SUBROUTINE Pack_Single

