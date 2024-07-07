! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reinitialises sea ice surface temperature

MODULE Rcf_Derv_Ice_Temp_Mod

! Description:
!     This subroutine is the initialises the sea ice temperature
!     on ice catagories when these are called for in stash
!     It sets all ice category pseudo levels to have temperatures
!     equal to the single level sea ice temperature (stash=49)
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_ICE_TEMP_MOD'

CONTAINS

SUBROUTINE Rcf_Derv_Ice_Temp( fields_out, field_count_out, hdr_out,  &
                            ice_temp_cat )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE um_stashcode_mod, ONLY: &
    stashcode_sea_ice_temp,    &
    stashcode_prog_sec

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal

USE UM_ParCore, ONLY: &
    mype

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER       :: fields_out(:)
TYPE( um_header_type), INTENT(IN) :: hdr_out
TYPE( field_type ), INTENT(INOUT) :: ice_temp_cat
INTEGER, INTENT(IN)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       ::  ice_temp


INTEGER                           ::  pos   ! position in array
INTEGER                           ::  i     ! loop index

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_DERV_ICE_TEMP'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Reinitialising ice surface temperature'
  CALL umPrint(umMessage,src='rcf_derv_ice_temp_mod')
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------
! GBM ice temperature; will abort if ice_temp not found

CALL Rcf_Locate( stashcode_prog_sec, stashcode_sea_ice_temp,       &
                 fields_out,field_count_out,pos)

ice_temp => fields_out(pos)
CALL Rcf_Alloc_Field( ice_temp )
CALL Rcf_Read_Field( ice_temp, hdr_out, decomp_rcf_output )


!----------------------------------------------------------------------
! Loop through ice_temp_cat
!----------------------------------------------------------------------

DO i = 1,ice_temp_cat % levels
  ice_temp_cat % DATA(:,i) = ice_temp % DATA(:,1)
END DO

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------
CALL Rcf_Dealloc_Field( ice_temp )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Derv_Ice_Temp
END MODULE Rcf_Derv_Ice_Temp_Mod
