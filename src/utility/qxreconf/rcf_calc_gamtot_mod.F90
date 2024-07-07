! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Set up Gamtot for large-scale hydrology

MODULE Rcf_Calc_Gamtot_Mod

! Subroutine Rcf_Calc_Gamtot
!
! Description:
!   Gets the total gamma function for topographic index for use
!   as a scale factor in the large-scale hydrology scheme.
!
! Method:
!   Calls routine CALC_GAMMA
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
!   Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Reconfiguration

 
USE timer_mod, ONLY: timer
 
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_CALC_GAMTOT_MOD'

CONTAINS

SUBROUTINE Rcf_Calc_Gamtot( fields_out, field_count_out, &
                            gamtot_pos, hdr_out )

USE calc_fsat_mod,         ONLY: calc_fsat

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE decomp_params, ONLY: &
    decomp_rcf_output

USE Rcf_Field_Type_Mod, ONLY: &
    Field_Type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE um_stashcode_mod, ONLY: &
    stashcode_gamtot,          &
    stashcode_vol_smc_sat,     &
    stashcode_ti_mean,         &
    stashcode_ti_sig,          &
    stashcode_prog_sec


USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal

USE nlstcall_mod, ONLY:   &
    LTimer

USE UM_ParCore, ONLY: &
    mype

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE Ereport_mod, ONLY: &
    Ereport

USE jules_hydrology_mod, ONLY: l_top

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER          :: fields_out(:)
TYPE( um_header_type ), INTENT(IN) :: hdr_out
INTEGER, INTENT(IN)                :: gamtot_pos
INTEGER, INTENT(IN)                :: field_count_out

! Local variables
TYPE( field_type ), POINTER          :: ti_mean
TYPE( field_type ), POINTER          :: ti_sig
TYPE( field_type ), POINTER          :: vol_smc_sat
TYPE( field_type ), POINTER          :: gamtot

REAL                                 :: dummy  &
                                       (fields_out(gamtot_pos) % level_size)


INTEGER                              :: soil_index  &
                                       (fields_out(gamtot_pos) % level_size)
INTEGER                              :: soil_pts
INTEGER                              :: SIZE
INTEGER                              :: pos  ! field position
INTEGER                              :: i    ! looper

INTEGER      :: ErrorStatus

CHARACTER (LEN=errormessagelength)           :: Cmessage
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_CALC_GAMTOT'

LOGICAL :: L_Gamtot = .TRUE.
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
IF (LTimer) CALL Timer( RoutineName, 3)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Setting up Gamtot from topographic variables'
  CALL umPrint(umMessage,src='rcf_calc_gamtot_mod')
END IF

IF (.NOT. l_top) THEN
  ErrorStatus = 10
  WRITE(CMessage, '(A)') 'Field only required in TOPMODEL-based hydrology'
  CALL Ereport ( RoutineName, ErrorStatus, CMessage)
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in:
!----------------------------------------------------------------------

gamtot => fields_out(gamtot_pos)

CALL Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,        &
                 fields_out, field_count_out, pos)
vol_smc_sat => fields_out(pos)
CALL Rcf_Alloc_Field( vol_smc_sat )
CALL Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_ti_mean,            &
                 fields_out, field_count_out, pos)
ti_mean => fields_out(pos)
CALL Rcf_Alloc_Field( ti_mean )
CALL Rcf_Read_Field( ti_mean, hdr_out, decomp_rcf_output )

CALL Rcf_Locate( stashcode_prog_sec, stashcode_ti_sig,             &
                 fields_out, field_count_out, pos)
ti_sig => fields_out(pos)
CALL Rcf_Alloc_Field( ti_sig )
CALL Rcf_Read_Field( ti_sig, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Set up soil index:
!----------------------------------------------------------------------
soil_pts=0
soil_index(:)=0
DO i=1 , vol_smc_sat % level_size
  IF (vol_smc_sat % DATA(i,1) > 0.0) THEN
    soil_pts = soil_pts + 1
    soil_index(soil_pts) = i
  END IF
END DO

!----------------------------------------------------------------------
! Calculate gamtot:
!----------------------------------------------------------------------

gamtot % DATA(:,1) = 0.0

CALL Calc_Fsat( L_gamtot,              &
                soil_pts,              &
                soil_index,            &
                ti_mean % level_size,  &
                ti_mean % DATA(:,1),   &
                ti_sig % DATA(:,1),    &
                dummy,                 &
                dummy,                 &
                gamtot % DATA(:,1),    &
                dummy,                 &
                dummy)

!----------------------------------------------------------------------
! Tidy Up
!----------------------------------------------------------------------

CALL Rcf_Dealloc_Field( vol_smc_sat )
CALL Rcf_Dealloc_Field( ti_mean )
CALL Rcf_Dealloc_Field( ti_sig )
IF (LTimer) CALL Timer( RoutineName, 4)

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Calc_Gamtot

END MODULE Rcf_Calc_Gamtot_Mod
