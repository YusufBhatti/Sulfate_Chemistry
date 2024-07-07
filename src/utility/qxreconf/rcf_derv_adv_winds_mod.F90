! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reinitialises advected winds as prognostic winds.

MODULE Rcf_Derv_Adv_Winds_Mod

!  Subroutine Rcf_Derv_Adv_Winds
!
! Description:
!   Derive advected wind components from actual wind components
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

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_DERV_ADV_WINDS_MOD'

CONTAINS

SUBROUTINE Rcf_Derv_Adv_Winds( stash_item, fields_out,            &
                               field_count_out, hdr_out,          &
                               advect_wind )

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE um_stashcode_mod, ONLY: &
    stashcode_u_adv,     stashcode_v_adv,         &
    stashcode_u,         stashcode_v,             &
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
TYPE( field_type ), INTENT(INOUT), TARGET :: advect_wind
INTEGER, INTENT(IN)               :: STASH_Item
INTEGER, INTENT(IN)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       ::  wind

INTEGER                           ::  pos   ! position in array
INTEGER                           ::  i,j,k ! loop index

CHARACTER (LEN=*), PARAMETER  :: RoutineName='RCF_DERV_ADV_WINDS'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  IF ( STASH_Item == stashcode_u_adv ) THEN
    WRITE(umMessage,*) 'Reinitialising Advected U Wind as U Wind'
    CALL umPrint(umMessage,src='rcf_derv_adv_winds_mod')
  ELSE IF ( STASH_Item == stashcode_v_adv ) THEN
    WRITE(umMessage,*) 'Reinitialising Advected V Wind as V Wind'
    CALL umPrint(umMessage,src='rcf_derv_adv_winds_mod')
  END IF
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in where available
!----------------------------------------------------------------------
! Prognostic U or V Wind; will abort if wind not found

IF ( STASH_Item == stashcode_u_adv ) THEN
  CALL Rcf_Locate(stashcode_prog_sec, stashcode_u, &
                  fields_out, field_count_out, pos)
ELSE IF ( STASH_Item == stashcode_v_adv ) THEN
  CALL Rcf_Locate(stashcode_prog_sec, stashcode_v, &
                  fields_out, field_count_out, pos)
END IF

wind => fields_out(pos)
CALL Rcf_Alloc_Field( wind )
CALL Rcf_Read_Field( wind, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Loop through wind levels
!----------------------------------------------------------------------

DO i = 1,wind % levels
  advect_wind % DATA(:,i) = wind % DATA(:,i)
END DO

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------

CALL Rcf_Dealloc_Field( wind )

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Derv_Adv_Winds
END MODULE Rcf_Derv_Adv_Winds_mod
