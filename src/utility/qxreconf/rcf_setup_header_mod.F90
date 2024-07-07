! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  sets up the output dump header

MODULE Rcf_Setup_Header_Mod
IMPLICIT NONE

!  Subroutine Rcf_Setup_Header - sets up the output dump header.
!
! Description:
!   The header for the output dump is constructed.
!
! Method:
!   The header is setup section by section - see UMDP F3 for
!   details.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SETUP_HEADER_MOD'

CONTAINS

SUBROUTINE Rcf_Setup_Header( Hdr_In, Hdr_Out )

USE Rcf_UMhead_Mod, ONLY: &
    Um_header_type,    &
    LenFixHd

USE Rcf_Setup_FixHd_Mod, ONLY: &
    Rcf_Setup_FixHd

USE Rcf_Setup_IntC_Mod, ONLY: &
    Rcf_Setup_IntC

USE Rcf_Setup_RealC_Mod, ONLY: &
    Rcf_Setup_RealC

USE Rcf_Setup_Lookup_Mod, ONLY: &
    Rcf_Setup_Lookup

USE Rcf_Grid_Type_Mod, ONLY: &
    Output_Grid

USE Rcf_AllocHdr_Mod, ONLY: &
    Rcf_AllocHdr

USE Rcf_Setup_LevDepC_Mod, ONLY: &
    Rcf_Setup_LevDepC

USE Rcf_Setup_RowDepC_Mod, ONLY: &
    Rcf_Setup_RowDepC

USE Rcf_Setup_ColDepC_Mod, ONLY: &
    Rcf_Setup_ColDepC

USE missing_data_mod, ONLY: rmdi

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE (Um_Header_Type), INTENT(IN)    :: Hdr_In
TYPE (Um_Header_Type), INTENT(INOUT) :: Hdr_Out

! Local arguments
INTEGER                              :: i          ! looper
INTEGER                              :: Len1_min   ! Min(Len1)
INTEGER                              :: Len2_min   ! Min(Len2)
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SETUP_HEADER'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!--------------------------------------------------------------
! Fixed Header
!--------------------------------------------------------------
ALLOCATE( Hdr_Out % FixHd( LenFixHd ) )
CALL Rcf_Setup_FixHd( Hdr_In, Hdr_Out )
CALL Rcf_AllocHdr( Hdr_Out )

!--------------------------------------------------------------
! Integer Constants
!--------------------------------------------------------------
CALL Rcf_Setup_IntC( Hdr_In, Hdr_Out )

!--------------------------------------------------------------
! Real Constants
!--------------------------------------------------------------
CALL Rcf_Setup_RealC( Hdr_In, Hdr_Out )

!-------------------------------------------------------------
! Below here things could be put into their own routines -
! I've not done this for simplicity, but it should be done if
! there is significant expansion of any of the sections.
!-------------------------------------------------------------

!-------------------------------------------------------------
! Level dependent constants
!-------------------------------------------------------------
! Initialise to RMDI
Hdr_Out % LevDepC(:,:) = rmdi

CALL Rcf_Setup_LevDepC( Hdr_Out, Hdr_In, Output_Grid )

!--------------------------------------------------------------
! row dependent constants
!--------------------------------------------------------------
CALL Rcf_Setup_RowDepC( Hdr_Out, Output_Grid )

!----------------------------------------------------------------
! Copy column dependent constants
!----------------------------------------------------------------
CALL Rcf_Setup_ColDepC( Hdr_Out, Output_Grid )

!-----------------------------------------------------------------
! Copy fields of constants
!-----------------------------------------------------------------
Hdr_Out % FldsOfC( : ) = Hdr_In % FldsOfC( : )

!-----------------------------------------------------------------
! Copy extra constants
!-----------------------------------------------------------------
DO i = 1, Hdr_Out % LenExtraC
  Hdr_Out % ExtraC( i ) = Hdr_In % ExtraC( i )
END DO

!-----------------------------------------------------------------
! Copy History block
!-----------------------------------------------------------------
Hdr_Out % HistFile( : ) = Hdr_In % HistFile( : )

!-----------------------------------------------------------------
! Copy Compressed field indexes - not needed since ocean removed.
!-----------------------------------------------------------------

!--------------------------------------------------------------
! Lookup Tables
!--------------------------------------------------------------
CALL Rcf_Setup_Lookup( Hdr_In, Hdr_Out )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Setup_Header
END MODULE Rcf_Setup_Header_Mod
