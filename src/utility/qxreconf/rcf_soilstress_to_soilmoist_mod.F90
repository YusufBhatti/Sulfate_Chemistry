! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Covert soil stress to moisture if soil stress was used in the
!  interpolation in the place of soil moisture

MODULE Rcf_Soilstress_to_soilmoist_Mod
IMPLICIT NONE

!  Subroutine Rcf_Soilstress_to_soilmoist
!
! Description:
!   Convert soil stress back to soil moisture after interpolation.
!
! Method:
!    1. Read soil stress, wilt, sat and cri from output dump
!    2. convert soil stress to Vol soil:
!         soil_vol=soilstress*(soil_cri-soil_wilt)+soil_wilt
!    3. Ensure SMC_vol is less or equal to soil_sat:
!         soil_vol=min(soil_vol, soil_sat)
!    4. Convert soil_vol to  soilmoist:
!         soilmoist=soil_vol*rhow*delta_z  ! delta_z is layer depth
!    5. Check for soil moisture being between min/max values and for
!         consistency with land usage fields is done later
!         (in Rcf_soil_Moist_Chk)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  ModuleName='RCF_SOILSTRESS_TO_SOILMOIST_MOD'

CONTAINS

SUBROUTINE Rcf_Soilstress_to_soilmoist( fields, field_count, &
           grid, decomp, hdr)

USE Rcf_Field_Type_Mod, ONLY: &
    field_type

USE Rcf_UMhead_Mod, ONLY: &
    um_header_type

USE Rcf_Grid_Type_Mod, ONLY: &
    grid_type

USE Rcf_Locate_Mod, ONLY: &
    Rcf_Locate

USE Rcf_Read_Field_Mod, ONLY: &
    Rcf_Read_Field

USE Rcf_Write_Field_Mod, ONLY: &
    Rcf_Write_Field

USE Rcf_Alloc_Field_Mod, ONLY: &
    Rcf_Alloc_Field,            &
    Rcf_DeAlloc_Field

USE um_stashcode_mod, ONLY: &
    stashcode_vol_smc_wilt,    &
    stashcode_vol_smc_cri,     &
    stashcode_vol_smc_sat,     &
    stashcode_soil_moist,      &
    stashcode_prog_sec

USE UM_ParCore, ONLY: &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,                &
    PrStatus_Normal


USE water_constants_mod, ONLY: rho_water

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
INTEGER                            :: pos_sat
INTEGER                            :: pos_wilt
INTEGER                            :: pos_smc
INTEGER                            :: pos_cri
INTEGER                            :: i
INTEGER                            :: j

TYPE( field_type ), POINTER        :: soil_moist
TYPE( field_type ), POINTER        :: wilt
TYPE( field_type ), POINTER        :: sat
TYPE( field_type ), POINTER        :: cri

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SOILSTRESS_TO_SOILMOIST'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (mype == 0 .AND. PrintStatus >= PrStatus_Normal) THEN
  WRITE(umMessage,'(A)') 'Converting soil stress to soil moisture...'
  CALL umPrint(umMessage,src='rcf_soilstress_to_soilmoist_mod')
END IF

!-------------------------------------------------------------------
! Locate required peramaters in output dump
!-------------------------------------------------------------------

CALL Rcf_Locate &
    ( stashcode_prog_sec,stashcode_vol_smc_sat,  &
     fields, field_count, pos_sat,  zero_ok_arg = .TRUE. )
CALL Rcf_Locate &
    ( stashcode_prog_sec,stashcode_vol_smc_wilt, &
      fields, field_count, pos_wilt,zero_ok_arg = .TRUE. )
CALL Rcf_Locate &
    ( stashcode_prog_sec,stashcode_soil_moist,   &
      fields, field_count, pos_smc,  zero_ok_arg = .TRUE. )
CALL Rcf_Locate &
    ( stashcode_prog_sec,stashcode_vol_smc_cri,  &
      fields, field_count, pos_cri,  zero_ok_arg = .TRUE. )

!--------------------------------------------------------------------
! Reset smc if appropriate
!--------------------------------------------------------------------

IF (pos_sat  /= 0  .AND. &
    pos_wilt /= 0 .AND. &
    pos_cri  /= 0 .AND. &
    pos_smc  /= 0 ) THEN

  soil_moist => fields( pos_smc )
  wilt       => fields( pos_wilt)
  sat        => fields( pos_sat )
  cri        => fields( pos_cri )


  CALL Rcf_Alloc_Field( soil_moist )
  CALL Rcf_Alloc_Field( wilt       )
  CALL Rcf_Alloc_Field( sat        )
  CALL Rcf_Alloc_Field( cri        )

  CALL Rcf_Read_Field( soil_moist, hdr, decomp )
  CALL Rcf_Read_Field( wilt ,      hdr, decomp )
  CALL Rcf_Read_Field( sat  ,      hdr, decomp )
  CALL Rcf_Read_Field( cri  ,      hdr, decomp )


  DO j = 1, soil_moist % levels
    DO i = 1, soil_moist % level_size

      ! convert smc_stress to VOL SMC
      soil_moist % data(i,j) = soil_moist % data(i,j)*                   &
                      (cri % data(i,1)-wilt % data(i,1))+wilt % data(i,1)


      ! calculate smc
      soil_moist % data(i,j) = soil_moist % data(i,j)*                   &
                               Grid % soil_depths(j)*rho_water
    END DO
  END DO


  IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage,'(A)')'Soil Moisture has been converted from soil stress.'
    CALL umPrint(umMessage,src='rcf_soilstress_to_soilmoist_mod')

  END IF


  !---------------------------------------------------------------------
  ! Write out converted field
  !---------------------------------------------------------------------

  CALL Rcf_Write_Field( soil_moist, hdr, decomp )

  CALL Rcf_Dealloc_Field( soil_moist )
  CALL Rcf_Dealloc_Field( wilt       )
  CALL Rcf_Dealloc_Field( sat        )
  CALL Rcf_Dealloc_Field( cri        )

END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Soilstress_to_soilmoist
END MODULE Rcf_Soilstress_to_soilmoist_Mod

