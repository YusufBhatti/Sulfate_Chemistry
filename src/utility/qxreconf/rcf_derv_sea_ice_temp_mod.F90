! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Derive sea ice temp from surface temp.

MODULE Rcf_Derv_Sea_Ice_Temp_Mod

!  Subroutine: Rcf_Derv_Sea_Ice_Temp
!
! Description:
!   Copies T* to sea ice temp.
!
! Method:
!   A basic copy with no calculation.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_SEA_ICE_TEMP_MOD'

CONTAINS

SUBROUTINE Rcf_Derv_Sea_Ice_Temp( fields_out, field_count_out,        &
                                  sea_ice_temp, hdr_out )

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE decomp_params, ONLY: &
    decomp_rcf_output

USE um_stashcode_mod, ONLY: &
    stashcode_tstar,           &
    stashcode_prog_sec

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
INTEGER,               INTENT(IN)    :: field_count_out
TYPE( field_type ),    TARGET        :: fields_out( field_count_out )
TYPE( field_type ),    INTENT(INOUT) :: sea_ice_temp
TYPE( um_header_type), INTENT(IN)    :: hdr_out

! Local variables
INTEGER                        :: pos
INTEGER                        :: i
TYPE( field_type ), POINTER    :: t_star

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_DERV_SEA_ICE_TEMP'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!-----------------------------------------------------------------
! Write out our action if appropriate
!-----------------------------------------------------------------
IF ( mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,*) 'Copying T* to Sea Ice Temp'
  CALL umPrint(umMessage,src='rcf_derv_sea_ice_temp_mod')
END IF

!-----------------------------------------------------------------
! Get T* from output dump
!-----------------------------------------------------------------
CALL Rcf_Locate( stashcode_prog_sec, stashcode_tstar,                &
                 fields_out, field_count_out, pos )
t_star => fields_out( pos )
CALL Rcf_Alloc_Field( t_star )
CALL Rcf_Read_Field( t_star, Hdr_Out, decomp_rcf_output )

!------------------------------------------------------------------
! Do the copy...
!------------------------------------------------------------------
sea_ice_temp % DATA( : , 1) = t_star % DATA( : , 1)

!-------------------------------------------------------------------
! clean up
!-------------------------------------------------------------------
CALL Rcf_DeAlloc_Field( t_star )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Derv_Sea_Ice_Temp
END MODULE Rcf_Derv_Sea_Ice_Temp_Mod
