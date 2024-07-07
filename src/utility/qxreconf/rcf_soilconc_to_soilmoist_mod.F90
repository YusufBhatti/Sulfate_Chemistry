! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Covert soil conc to moisture if soil conc was used in the
!  interpolation in the place of soil moisture

MODULE rcf_soilconc_to_soilmoist_mod

!  Subroutine Rcf_Soilconc_to_soilmoist

! Description:
!   Convert soil conc back to soil moisture after interpolation.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_SOILCONC_TO_SOILMOIST_MOD'

CONTAINS

SUBROUTINE rcf_soilconc_to_soilmoist( fields, field_count, &
           grid, decomp, hdr)

USE rcf_field_type_mod, ONLY: &
    field_type

USE rcf_umhead_mod, ONLY: &
    um_header_type

USE rcf_grid_type_mod, ONLY: &
    grid_type

USE rcf_locate_mod, ONLY: &
    rcf_locate

USE rcf_read_field_mod, ONLY: &
    rcf_read_field

USE rcf_write_field_mod, ONLY: &
    rcf_write_field

USE rcf_alloc_field_mod, ONLY: &
    rcf_alloc_field,            &
    rcf_dealloc_field

USE um_stashcode_mod, ONLY: &
    stashcode_soil_moist,      &
    stashcode_prog_sec

USE um_parcore, ONLY: &
    mype

USE umPrintMgr

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Arguments
TYPE( field_type ), POINTER   :: fields(:)
TYPE( grid_type )             :: grid
TYPE( um_header_type )        :: hdr
INTEGER                       :: field_count
INTEGER                       :: decomp

! Local variables
INTEGER                            :: pos_smc
INTEGER                            :: i
INTEGER                            :: j

TYPE( field_type ), POINTER        :: soil_moist

REAL                               :: smc_min
REAL                               :: smc_max
CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SOILCONC_TO_SOILMOIST'
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. printstatus >= prstatus_normal) THEN
  CALL umPrint(' Converting soil-moisture concentration back to SMCL', &
      src='rcf_soilconc_to_soilmoist_mod')
END IF

!-------------------------------------------------------------------
! Locate required parameters in output dump
!-------------------------------------------------------------------

CALL rcf_locate &
    ( stashcode_prog_sec,stashcode_soil_moist,   &
      fields, field_count, pos_smc,zero_ok_arg = .TRUE. )

!--------------------------------------------------------------------
! Reset smc if appropriate
!--------------------------------------------------------------------

IF (pos_smc /= 0 ) THEN

  soil_moist => fields( pos_smc )

  CALL rcf_alloc_field( soil_moist )

  CALL rcf_read_field( soil_moist, hdr, decomp )

  ! calculate smc
  DO j = 1, soil_moist % levels
    soil_moist % DATA(:,j) =   soil_moist % DATA(:,j) &
                             * grid % soil_depths(j)
  END DO

  !---------------------------------------------------------------------
  ! Write out converted field
  !---------------------------------------------------------------------

  CALL rcf_write_field( soil_moist, hdr, decomp )

  CALL rcf_dealloc_field( soil_moist )

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE rcf_soilconc_to_soilmoist
END MODULE rcf_soilconc_to_soilmoist_mod

