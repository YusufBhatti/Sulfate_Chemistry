! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Reinitialises frozen/unfrozen soil moisture

MODULE Rcf_Freeze_Soil_Mod

! Description:
!     This subroutine is the interface to FREEZE_SOIL which
!     reinitialises frozen/unfrozen soil moisture fields after
!     soil moisture/temperature has been updated from ancillary.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_FREEZE_SOIL_MOD'

CONTAINS

SUBROUTINE Rcf_Freeze_Soil( fields_out, field_count_out, hdr_out,  &
                            grid_out )

USE freeze_soil_mod, ONLY: freeze_soil

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Grid_Type_Mod, ONLY: &
    Grid_Type

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_Dealloc_Field

USE um_stashcode_mod, ONLY: &
    stashcode_soil_moist,      &
    stashcode_soil_temp,       &
    stashcode_vol_smc_sat,     &
    stashcode_soil_suction,    &
    stashcode_clapp_hb,        &
    stashcode_unfrozen_soil,   &
    stashcode_frozen_soil,     &
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

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

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
TYPE( grid_type ), INTENT(IN)     :: grid_out
INTEGER, INTENT(IN)               :: field_count_out

! Internal variables
TYPE( field_type ), POINTER       ::  soil_moist
TYPE( field_type ), POINTER       ::  soil_temp
TYPE( field_type ), POINTER       ::  vol_smc_sat
TYPE( field_type ), POINTER       ::  soil_suction
TYPE( field_type ), POINTER       ::  clapp_hb
TYPE( field_type ), POINTER       ::  unfrozen_soil
TYPE( field_type ), POINTER       ::  frozen_soil

INTEGER                           ::  pos   ! position in array
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

CHARACTER(LEN=*), PARAMETER :: RoutineName='RCF_FREEZE_SOIL'

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!----------------------------------------------------------------------
! Write out action if appropriate
!----------------------------------------------------------------------
IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal ) THEN
  WRITE(umMessage,*) 'Reinitialising frozen/unfrozen soil moisture'
  CALL umPrint(umMessage,src='rcf_freeze_soil_mod')
END IF

!----------------------------------------------------------------------
! Find required fields in output dump and read them in
!----------------------------------------------------------------------
! Soil moisture
CALL Rcf_Locate( stashcode_prog_sec, stashcode_soil_moist,           &
                 fields_out, field_count_out, pos)
soil_moist => fields_out(pos)
CALL Rcf_Alloc_Field( soil_moist )
CALL Rcf_Read_Field( soil_moist, hdr_out, decomp_rcf_output )

! Soil temperature
CALL Rcf_Locate( stashcode_prog_sec, stashcode_soil_temp,            &
                 fields_out, field_count_out, pos)
soil_temp => fields_out(pos)
CALL Rcf_Alloc_Field( soil_temp )
CALL Rcf_Read_Field( soil_temp, hdr_out, decomp_rcf_output )

! Volume SMC at Saturation
CALL Rcf_Locate( stashcode_prog_sec, stashcode_vol_smc_sat,          &
                 fields_out, field_count_out,pos)
vol_smc_sat => fields_out(pos)
CALL Rcf_Alloc_Field( vol_smc_sat )
CALL Rcf_Read_Field( vol_smc_sat, hdr_out, decomp_rcf_output )

! Saturated Soil Water Suction
CALL Rcf_Locate( stashcode_prog_sec, stashcode_soil_suction,         &
                 fields_out, field_count_out, pos)
soil_suction => fields_out(pos)
CALL Rcf_Alloc_Field( soil_suction )
CALL Rcf_Read_Field( soil_suction , hdr_out, decomp_rcf_output )

! Clapp-Hornberger "B" Coefficient
CALL Rcf_Locate( stashcode_prog_sec, stashcode_clapp_hb,             &
                 fields_out, field_count_out, pos)
clapp_hb => fields_out(pos)
CALL Rcf_Alloc_Field( clapp_hb )
CALL Rcf_Read_Field( clapp_hb , hdr_out, decomp_rcf_output )

! Unfrozen soil moisture
CALL Rcf_Locate( stashcode_prog_sec, stashcode_unfrozen_soil,        &
                 fields_out, field_count_out, pos)
unfrozen_soil => fields_out(pos)
CALL Rcf_Alloc_Field( unfrozen_soil )
CALL Rcf_Read_Field( unfrozen_soil, hdr_out, decomp_rcf_output )

! Frozen soil moisture
CALL Rcf_Locate( stashcode_prog_sec, stashcode_frozen_soil,          &
                 fields_out, field_count_out, pos)
frozen_soil => fields_out(pos)
CALL Rcf_Alloc_Field( frozen_soil )
CALL Rcf_Read_Field( frozen_soil, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Call Freeze_Soil to do the calculation
!----------------------------------------------------------------------
CALL Freeze_Soil( soil_temp % level_size, grid_out % sm_levels,    &
                  clapp_hb % DATA,        grid_out % soil_depths,  &
                  soil_suction % DATA,    soil_moist % DATA,       &
                  soil_temp % DATA,       vol_smc_sat % DATA,      &
                  unfrozen_soil % DATA,   frozen_soil % DATA   )

!----------------------------------------------------------------------
! Write out corrected soil frozen/unfrozen fractions
!----------------------------------------------------------------------
CALL Rcf_Write_Field( unfrozen_soil, hdr_out, decomp_rcf_output )
CALL Rcf_Write_Field( frozen_soil, hdr_out, decomp_rcf_output )

!----------------------------------------------------------------------
! Clear up dynamic memory used
!----------------------------------------------------------------------
CALL Rcf_Dealloc_Field( frozen_soil )
CALL Rcf_Dealloc_Field( unfrozen_soil )
CALL Rcf_Dealloc_Field( clapp_hb )
CALL Rcf_Dealloc_Field( soil_suction )
CALL Rcf_Dealloc_Field( vol_smc_sat )
CALL Rcf_Dealloc_Field( soil_temp )
CALL Rcf_Dealloc_Field( soil_moist )


IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Freeze_Soil
END MODULE Rcf_Freeze_Soil_Mod
