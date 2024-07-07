! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Allocates space in the header data-type

MODULE Rcf_AllocHdr_Mod

!  Subroutine Rcf_AllocHdr    - allocates space for the header
!
! Description:
!   Allocates space for all the constituant parts of the header
!   data-type
!
! Method:
!   Space is allocated according to sizes in the fixed-length header.
!   If present, FFLookup is the size of Len2Lookup - ie the actual
!   rather than FieldsFile maximum size.
!
!   Adapted from VAR code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_ALLOCHDR_MOD'

CONTAINS

SUBROUTINE Rcf_AllocHdr( UMhdr, FFLookup )

USE Rcf_UMhead_Mod        ! All of it is used...

USE Rcf_HeadAddress_Mod   ! Supplies various FH_* header constants.
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE
TYPE (UM_header_type), INTENT(INOUT) :: UMhdr        ! Data to Allocate
INTEGER, OPTIONAL                    :: FFLookup     ! Len2Lookup for

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_ALLOCHDR'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
                                                   ! FieldsFile

!-----------------------------------------------------------------------
! Section 1.  Set up lengths of header arrays from Fixed-Length-Header
!-----------------------------------------------------------------------

IF ( UMhdr % FixHd(FH_IntCStart) > 0 ) THEN          ! Integer Consts
  UMhdr % LenIntC = UMhdr % FixHd(FH_IntCSize)
ELSE                                                 ! No Integer Consts
  UMhdr % LenIntC = 0
END IF

IF ( UMhdr % FixHd(FH_RealCStart) > 0 ) THEN         ! Real Consts
  UMhdr % LenRealC = UMhdr % FixHd(FH_RealCSize)
ELSE                                                 ! No Real Consts
  UMhdr % LenRealC = 0
END IF

IF ( UMhdr % FixHd(FH_LevDepCStart) > 0 ) THEN          ! Lev Dep Consts
  UMhdr % Len1LevDepC = UMhdr % FixHd(FH_LevDepCSize1)
  UMhdr % Len2LevDepC = UMhdr % FixHd(FH_LevDepCSize2)
ELSE                                                    ! No Level Dep
  UMhdr % Len1LevDepC = 0                               ! Consts
  UMhdr % Len2LevDepC = 0
END IF

IF ( UMhdr % FixHd(FH_RowDepCStart) > 0 ) THEN          ! Row Dep Consts
  UMhdr % Len1RowDepC = UMhdr % FixHd(FH_RowDepCSize1)
  UMhdr % Len2RowDepC = UMhdr % FixHd(FH_RowDepCSize2)
ELSE                                                    ! No Row Dep
  UMhdr % Len1RowDepC = 0                               ! Consts
  UMhdr % Len2RowDepC = 0
END IF

IF ( UMhdr % FixHd(FH_ColDepCStart) > 0 ) THEN          ! Col Dep Consts
  UMhdr % Len1ColDepC = UMhdr % FixHd(FH_ColDepCSize1)
  UMhdr % Len2ColDepC = UMhdr % FixHd(FH_ColDepCSize2)
ELSE                                                    ! No Column Dep
  UMhdr % Len1ColDepC = 0                               ! Consts
  UMhdr % Len2ColDepC = 0
END IF


IF ( UMhdr % FixHd(FH_FldsOfCStart) > 0 ) THEN          ! Fields of
  UMhdr % Len1FldsOfC = UMhdr % FixHd(FH_FldsOfCSize1)  ! Consts
  UMhdr % Len2FldsOfC = UMhdr % FixHd(FH_FldsOfCSize2)
ELSE                                                    ! No Fields of
  UMhdr % Len1FldsOfC = 0                               ! Consts
  UMhdr % Len2FldsOfC = 0
END IF

IF ( UMhdr % FixHd(FH_ExtraCStart) > 0 ) THEN         ! Extra Consts
  UMhdr % LenExtraC = UMhdr % FixHd(FH_ExtraCSize)
ELSE                                                  ! No Extra Consts
  UMhdr % LenExtraC = 0
END IF

IF ( UMhdr % FixHd(FH_HistStart) > 0 ) THEN         ! Temp History
  UMhdr % LenHistFile = UMhdr % FixHd(FH_HistSize)  ! File present
ELSE
  UMhdr % LenHistFile = 0                           ! No Temp History
END IF

IF ( UMhdr % FixHd(FH_CompFldI1Start) > 0 ) THEN         ! Comp Field
  UMhdr % LenCompFldI1 = UMhdr % FixHd(FH_CompFldI1Size) ! Index 1
ELSE
  UMhdr % LenCompFldI1 = 0                               ! No Comp Field
END IF                                                    ! Index 1

IF ( UMhdr % FixHd(FH_CompFldI2Start) > 0 ) THEN         ! Comp Field
  UMhdr % LenCompFldI2 = UMhdr % FixHd(FH_CompFldI2Size) ! Index 2
ELSE
  UMhdr % LenCompFldI2 = 0                               ! No Comp Field
END IF                                                    ! Index 2

IF ( UMhdr % FixHd(FH_CompFldI3Start) > 0 ) THEN         ! Comp Field
  UMhdr % LenCompFldI3 = UMhdr % FixHd(FH_CompFldI3Size) ! Index 3
ELSE
  UMhdr % LenCompFldI3 = 0                               ! No Comp Field
END IF                                                    ! Index 3

IF ( UMhdr % FixHd(FH_LookupStart) > 0 ) THEN            ! Lookup Table
  UMhdr % Len1Lookup = UMhdr % FixHd(FH_LookupSize1)     ! present.

  IF ( PRESENT( FFLookup ) ) THEN             ! FieldsFile
    UMhdr % FixHd(FH_LookupSize2) = FFLookup  ! Change num of lookup
                                              ! entries for fixhd
  END IF

  UMhdr % Len2Lookup = UMhdr % FixHd(FH_LookupSize2)

ELSE                                                     ! No Lookup
  UMhdr % Len1Lookup = 0                                 ! Table.
  UMhdr % Len2Lookup = 0
END IF

IF ( UMhdr % FixHd(FH_DataStart) > 0 ) THEN              ! Data
  UMhdr % StartData = UMhdr % FixHd(FH_DataStart)        ! present.
  UMhdr % LenData  = UMhdr % FixHd(FH_DataSize)
ELSE                                                     ! No Data.
  UMhdr % StartData = 0
  UMhdr % LenData   = 0
END IF


!-----------------------------------------------------------------------
! Section 2.  Allocate space for each UM header array
!           (Most arrays dimensioned by ...+1 for F77 dynamic
!           allocation in UM.)
!
! Don't check allocation status, as its better for program to
! fail internally
!-----------------------------------------------------------------------

ALLOCATE ( UMhdr % IntC      (UMhdr % LenIntC))
ALLOCATE ( UMhdr % RealC     (UMhdr % LenRealC))

IF ( UMhdr % Len1LevDepC  /= 0 .AND. UMhdr % Len2LevDepC  /= 0 ) THEN
  ALLOCATE ( UMhdr % LevDepC   (UMhdr % Len1LevDepC, &
                                UMhdr % Len2LevDepC ))
ELSE
  ALLOCATE ( UMhdr % LevDepC (1,1)  )
END IF

!Use 2 dimensional Row and Column dependent consts.
IF ( UMhdr % Len1RowDepC  /= 0 .AND. UMhdr % Len2RowDepC  /= 0 ) THEN
  ALLOCATE ( UMhdr % RowDepC   (UMhdr % Len1RowDepC, &
                                UMhdr % Len2RowDepC ))
ELSE
  ALLOCATE ( UMhdr % RowDepC (1,1)  )
END IF

IF ( UMhdr % Len1ColDepC  /= 0 .AND. UMhdr % Len2ColDepC  /= 0 ) THEN
  ALLOCATE ( UMhdr % ColDepC   (UMhdr % Len1ColDepC, &
                                UMhdr % Len2ColDepC ))
ELSE
  ALLOCATE ( UMhdr % ColDepC (1,1)  )
END IF
ALLOCATE ( UMhdr % FldsOfC   (UMhdr % Len1FldsOfC * UMhdr % &
                                                    Len2FldsOfC + 1))
ALLOCATE ( UMhdr % ExtraC    (UMhdr % LenExtraC + 1))
ALLOCATE ( UMhdr % HistFile  (UMhdr % LenHistFile + 1))
ALLOCATE ( UMhdr % CompFldI1 (UMhdr % LenCompFldI1 + 1))
ALLOCATE ( UMhdr % CompFldI2 (UMhdr % LenCompFldI2 + 1))
ALLOCATE ( UMhdr % CompFldI3 (UMhdr % LenCompFldI3 + 1))
ALLOCATE ( UMhdr % Lookup    (UMhdr % Len1Lookup, UMhdr % Len2Lookup))

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
RETURN
END SUBROUTINE Rcf_AllocHdr
END MODULE Rcf_AllocHdr_Mod
