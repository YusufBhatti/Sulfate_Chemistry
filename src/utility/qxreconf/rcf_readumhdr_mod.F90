! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Read in the header from a UM dump.

MODULE Rcf_ReadUMhdr_Mod

!  Subroutine Rcf_ReadUMhdr - read a header from the dump
!
! Description:
!   Read in a model header from a UM dump.
!
! Method:
!   Call READ_FLH (UM routine) to read in Fixed-length Header.
!   Call Rcf_Allochdr to allocate array sizes for whole header, and
!   Call READHEAD (UM routine) to read in whole header.
!
!  Based on VAR code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_READUMHDR_MOD'

CONTAINS

SUBROUTINE Rcf_ReadUMhdr (  UMhdr )

USE Rcf_UMhead_Mod, ONLY: &
    LenFixHd,          &
    UM_Header_type

USE Ereport_Mod, ONLY: &
    Ereport

USE Rcf_AllocHdr_Mod, ONLY: &
    Rcf_AllocHdr

USE io
USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim
IMPLICIT NONE

! Subroutine arguments
TYPE (UM_header_type), INTENT (INOUT) ::  UMhdr ! UM header from dump

! Local constants
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_READUMHDR'

! Local variables
INTEGER             ::  ReturnCode   ! Return code from UM routines
INTEGER             ::  WordAddress  ! Position on file, used in SETPOS

INTEGER             ::  FFLookup
INTEGER             ::  Len_IO
INTEGER, ALLOCATABLE::  Lookup(:)
REAL                ::  a_io
CHARACTER (LEN=errormessagelength)  ::  Cmessage

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------------
!  Section 0.1:  Initialisations
!-----------------------------------------------------------------------

CMessage   = ' '
ReturnCode = 0


!-----------------------------------------------------------------------
!  Section 1.    Allocate & read in Fixed-Length Header from Model
!                dump file
!-----------------------------------------------------------------------

ALLOCATE (UMhdr % FixHd(LenFixHd))

WordAddress = 0

CALL setpos (UMhdr % UnitNum, WordAddress, ReturnCode)  ! Start of dump

IF (ReturnCode /= 0) THEN
  Cmessage = 'SETPOS failure'
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF

! DEPENDS ON: read_flh
CALL read_flh (                   &
  UMhdr % UnitNum,  &   ! in
  UMhdr % FixHd,    &   ! out
  LenFixHd,         &   ! in
  ReturnCode,       &   ! out
  CMessage )            ! out

IF (ReturnCode /= 0 ) THEN
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF


!-----------------------------------------------------------------------
!  Section 2.    Allocate space for the UM header from the Model
!                dump file
!-----------------------------------------------------------------------

! First need to count Lookups if dealing with a FieldsFile.
! Fieldsfiles produced by the UM use a system of reserved headers.
! The user estimates how many lookup headers they will need before
! runtime.  This often leads to empty lookups as not all of the
! reserved headers will be used, the number of used lookups will not
! always be the same as the second dimension of the lookup table (fixed
! length header 152).  For these files the recon must count the number
! of used lookups.  The UM fills empty lookup entries with -99.
! For other programs (such as the Mule utilities) the number of lookups
! will correspond exactly to the second dimension of the lookup table.
IF ( Umhdr % FixHd( 5 ) == 3 ) THEN        ! FieldsFile

  ! Check if Lookups exist first
  FFLookup = 0
  IF ( Umhdr % FixHd( 150 ) > 0 ) THEN

    CALL setpos( UMhdr % UnitNum, UMhdr % FixHd( 150 ), ReturnCode )

    IF ( ReturnCode /= 0 ) THEN
      Cmessage = 'SETPOS failure'
      CALL Ereport( RoutineName, ReturnCode, Cmessage )
    END IF

    ALLOCATE( Lookup( UMhdr % FixHd( 151 ) ) )
    Lookup(:) = 0

    ! Loop over the number of lookups specified by fixed header 152
    DO WHILE ( FFLookup < UMhdr % FixHd( 152 ) )
      CALL buffin( UMhdr % UnitNum, Lookup, UMhdr % FixHd( 151 ), &
              Len_IO, a_IO )
      ! -99 indicates the start of a block of empty lookups
      IF ( Lookup(1) == -99 ) THEN
        EXIT
      END IF
      FFLookup = FFLookup + 1
    END DO

    DEALLOCATE( Lookup )
  END IF

  CALL Rcf_AllocHdr( UMhdr, FFLookup )
ELSE
  CALL Rcf_AllocHdr ( UMhdr )        ! inout
END IF


!-----------------------------------------------------------------------
!  Section 3.    Read in complete UM header from the Model dump file
!-----------------------------------------------------------------------

WordAddress = 0

CALL setpos (UMhdr % UnitNum, WordAddress, ReturnCode)  ! Start of dump

IF (ReturnCode /= 0) THEN
  Cmessage = 'SETPOS failure'
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF

! DEPENDS ON: readhead
CALL readhead (                                                &
  UMhdr % UnitNum,                                             &
  UMhdr % FixHd,     LenFixHd,                                 &
  UMhdr % IntC,      UMhdr % LenIntC,                          &
  UMhdr % RealC,     UMhdr % LenRealC,                         &
  UMhdr % LevDepC,   UMhdr % Len1LevDepC, UMhdr % Len2LevDepC, &
  UMhdr % RowDepC,   UMhdr % Len1RowDepC, UMhdr % Len2RowDepC, &
  UMhdr % ColDepC,   UMhdr % Len1ColDepC, UMhdr % Len2ColDepC, &
  UMhdr % FldsOfC,   UMhdr % Len1FldsOfC, UMhdr % Len2FldsOfC, &
  UMhdr % ExtraC,    UMhdr % LenExtraC,                        &
  UMhdr % HistFile,  UMhdr % LenHistFile,                      &
  UMhdr % CompFldI1, UMhdr % LenCompFldI1,                     &
  UMhdr % CompFldI2, UMhdr % LenCompFldI2,                     &
  UMhdr % CompFldI3, UMhdr % LenCompFldI3,                     &
  UMhdr % Lookup,    UMhdr % Len1Lookup,  UMhdr % Len2Lookup,  &
  UMhdr % LenData,                                             &
  UMhdr % StartData,                                           &
  ReturnCode,        CMessage )

IF (ReturnCode /= 0) THEN
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF

WordAddress = 0

CALL setpos (UMhdr % UnitNum, WordAddress, ReturnCode)  ! Start of dump

IF (ReturnCode /=  0) THEN
  Cmessage = 'SETPOS failure'
  CALL Ereport( RoutineName, ReturnCode, Cmessage )
END IF
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

RETURN
END SUBROUTINE  Rcf_ReadUMhdr

END MODULE Rcf_ReadUMhdr_Mod
