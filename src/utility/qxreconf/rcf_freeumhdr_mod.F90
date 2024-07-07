! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!   Clears space from UM_header_type

MODULE Rcf_FreeUMhdr_Mod
IMPLICIT NONE

!  Subroutine Rcf_FreeUMhdr - frees up used arrays etc.
!
! Description:
!   Allocated pointers to hold header components are released and
!   sizes set to MDI as required.
!
! Method:
!   Deallocate Arrays, Set sizes to IMDI.
!
! Inspired by VAR code.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.


CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_FREEUMHDR_MOD'

CONTAINS

SUBROUTINE Rcf_FreeUMhdr ( UMhdr )

USE Rcf_UMhead_Mod, ONLY: &
    UM_header_type
USE missing_data_mod, ONLY: &
    imdi
USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (UM_header_type), INTENT(INOUT) ::  UMhdr

! Local variables
CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_FREEUMHDR'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------
! 1. Set constants to IMDI
!-----------------------------------------------------------------
UMhdr % LenIntC      = imdi
UMhdr % LenRealC     = imdi
UMhdr % Len1LevDepC  = imdi
UMhdr % Len2LevDepC  = imdi
UMhdr % Len1RowDepC  = imdi
UMhdr % Len2RowDepC  = imdi
UMhdr % Len1ColDepC  = imdi
UMhdr % Len2ColDepC  = imdi
UMhdr % Len1FldsOfC  = imdi
UMhdr % Len2FldsOfC  = imdi
UMhdr % LenExtraC    = imdi
UMhdr % LenHistFile  = imdi
UMhdr % LenCompFldI1 = imdi
UMhdr % LenCompFldI2 = imdi
UMhdr % LenCompFldI3 = imdi
UMhdr % Len1Lookup   = imdi
UMhdr % Len2Lookup   = imdi
UMhdr % StartData    = imdi
UMhdr % LenData      = imdi
UMhdr % NumFlds      = imdi
UMhdr % MaxFldSize   = imdi
UMhdr % UnitNum      = imdi

!--------------------------------------------------------------
! 2. Deallocate header information
!--------------------------------------------------------------

IF (ALLOCATED( UMhdr % FixHd     )) DEALLOCATE( UMhdr % FixHd     )
IF (ALLOCATED( UMhdr % IntC      )) DEALLOCATE( UMhdr % IntC      )
IF (ALLOCATED( UMhdr % CompFldI1 )) DEALLOCATE( UMhdr % CompFldI1 )
IF (ALLOCATED( UMhdr % CompFldI2 )) DEALLOCATE( UMhdr % CompFldI2 )
IF (ALLOCATED( UMhdr % CompFldI3 )) DEALLOCATE( UMhdr % CompFldI3 )
IF (ALLOCATED( UMhdr % Lookup    )) DEALLOCATE( UMhdr % Lookup    )
IF (ALLOCATED( UMhdr % RealC     )) DEALLOCATE( UMhdr % RealC     )
IF (ALLOCATED( UMhdr % LevDepC   )) DEALLOCATE( UMhdr % LevDepC   )
IF (ALLOCATED( UMhdr % RowDepC   )) DEALLOCATE( UMhdr % RowDepC   )
IF (ALLOCATED( UMhdr % ColDepC   )) DEALLOCATE( UMhdr % ColDepC   )
IF (ALLOCATED( UMhdr % FldsOfC   )) DEALLOCATE( UMhdr % FldsOfC   )
IF (ALLOCATED( UMhdr % ExtraC    )) DEALLOCATE( UMhdr % ExtraC    )
IF (ALLOCATED( UMhdr % HistFile  )) DEALLOCATE( UMhdr % HistFile  )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_FreeUMhdr

END MODULE Rcf_FreeUMhdr_Mod
