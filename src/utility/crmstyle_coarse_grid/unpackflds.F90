! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Routines for packing and unpacking

SUBROUTINE UnPackFlds( PackType,        &  ! in
                       NumCols,         &  ! in
                       NumRows,         &  ! in
                       Field,           &  ! inout
                       ErrorStatus )       ! inout

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
  LSMField
USE mask_compression, ONLY: expand_from_mask
USE wgdos_packing_mod, ONLY: wgdos_expand_field
USE Err_Mod, ONLY:        &
  StatusOK,               &
  StatusWarning
USE ereport_mod, ONLY: ereport

USE errormessagelength_mod, ONLY: errormessagelength

USE packing_codes_mod, ONLY: PC_BitMask_CompressType, PC_LandMask_Compression, &
  PC_SeaMask_Compression, PC_RunLength_Packing, PC_Cray32_Packing,             &
  PC_WGDOS_Packing

IMPLICIT NONE

! Subroutine Arguments:
INTEGER, INTENT(IN) :: PackType
INTEGER, INTENT(IN) :: NumCols
INTEGER, INTENT(IN) :: NumRows

TYPE(PP_Field_type), INTENT(INOUT) :: Field
INTEGER, INTENT(INOUT) :: ErrorStatus

! Local Constants:
CHARACTER(LEN=*), PARAMETER :: RoutineName = "UnPackFlds"

! Local Variables:
INTEGER :: FldSize
INTEGER :: idum                     ! Dummy integer
INTEGER :: NumWords                 ! Size of compressed field
INTEGER :: int_data(Field % Hdr % lblrec) ! Integer array expected by packing
INTEGER :: fake_stash_item          ! Dummy item number used in packing

REAL :: WorkArray(NumCols*NumRows,1) ! array used for un_packing

LOGICAL :: PackMask        ! Do we have some sort of mask
LOGICAL :: LandMask        ! Is it a landmask?
INTEGER :: Num_Land_Points ! Number of land points

! End of header --------------------------------------------------------

IF ( ErrorStatus /= StatusOK ) THEN
  ! Previous error - do not proceed
  GO TO 9999
END IF

NumWords = Field % Hdr % lblrec
FldSize  = NumCols * NumRows
PackMask = MOD(Field % Hdr % LBPack/10,10) == PC_BitMask_CompressType
LandMask = MOD(Field % Hdr % LBPack/100,10) == PC_LandMask_Compression
fake_stash_item = 0

IF ( PackType == PC_WGDOS_Packing ) THEN
  ! WGDOS packing
  int_data = TRANSFER(Field % RData, int_data)
  CALL wgdos_expand_field( WorkArray,  FldSize, int_data,                 &
                           NumWords,   idum,    NumCols,                  &
                           NumRows,    idum,    Field % Hdr % bmdi,       &
                           fake_stash_item,     ErrorStatus)
  ! StashItem not available, so set it to 0
  IF ( ErrorStatus /= StatusOK ) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, &
                  "wgdos_expand_field: Reported errors" )
    ErrorStatus = StatusWarning
    GO TO 9999
  END IF
ELSE IF ( PackType == PC_Cray32_Packing ) THEN
  ! 32 bit packing (in dumps)
  WorkArray = Field % RData
  ! DEPENDS ON: expand32b
  CALL expand32b( numwords , WorkArray,  Field % Hdr % LBuser7 )

ELSE IF ( PackType == PC_RunLength_Packing ) THEN
  ! Run length encoded data
   !CALL RUNLEN_DECODE( WorkArray, FldSize, Field % RData, NumWords,  &
   !                    Field % Hdr % BMDI, ErrorStatus, ErrMessage)
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
                "RUNLEN packing not supported." )
  ErrorStatus = StatusWarning
  GO TO 9999

ELSE
  ErrorStatus = StatusWarning

  CALL EReport( RoutineName, ErrorStatus, &
                "UNPACKFLDS - packing type not supported." )
  ErrorStatus = StatusWarning
  GO TO 9999

END IF

Field % RData = WorkArray

! We now unpack land mask

IF (PackMask) THEN
  ! Check we have LSM field available
  IF (.NOT. ASSOCIATED(LSMField % Rdata)) THEN
    ErrorStatus = StatusWarning

    CALL EReport( RoutineName, ErrorStatus, &
                  "LSM required to unpack land packed fields." )
    ErrorStatus = StatusWarning
  ELSE
    IF (LandMask) THEN
      CALL expand_from_mask(Field % RData, WorkArray, LSMField % RData /= 0.0, &
                            FldSize, Num_Land_Points)
      ! Data unpacked land
      Field % Hdr % lbpack = Field % Hdr % lbpack - 100*PC_LandMask_Compression
    ELSE
      ! For sea-mask packing we can just reverse the land packing mask.
      CALL expand_from_mask(Field % RData, WorkArray, LSMField % RData == 0.0, &
                            FldSize, Num_Land_Points)
      ! Data unpacked sea
      Field % Hdr % lbpack = Field % Hdr % lbpack - 100*PC_SeaMask_Compression
    END IF
    ! Data unpacked mask
    Field % Hdr % lbpack = Field % Hdr % lbpack - 10*PC_BitMask_CompressType
  END IF
END IF

Field % Hdr % lbpack = Field % Hdr % lbpack - PackType ! Data unpacked
Field % Hdr % lblrec = FldSize
Field % Hdr % LBUser1 = 1  ! Data type is now real
Field % Hdr % bacc = -99.0


9999 CONTINUE


END SUBROUTINE UnPackFlds

