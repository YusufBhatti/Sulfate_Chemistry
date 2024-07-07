! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  Checks that soil moisture is appropriately set for vegetation

MODULE Rcf_Soil_Moist_Chk_Mod

!  Subroutine Rcf_Soil_Moist_Chk
!
! Description:
!   Ensures that soil moisture is consistent with vegetation fields.
!
! Method:
!    Set Soil Moisture Content (SMC) = 0 for land-ice points
!    Ensure SMC is 0.1 * volumetric wilt as minimum at non-land-ice points
!    (0.1 is fairly arbitary but matches the nudging scheme - see
!     SURF documentation)
!    Ensure SMC is saturation as maximum at non-land-ice points
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_SOIL_MOIST_CHK_MOD'

CONTAINS

SUBROUTINE Rcf_Soil_Moist_Chk( fields, field_count, grid, decomp, hdr, &
                                l_changed )

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
    stashcode_vol_smc_sat,     &
    stashcode_soil_moist,      &
    stashcode_prog_sec

USE UM_ParCore, ONLY: &
    nproc,                  &
    mype

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
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
TYPE( field_type ), POINTER        :: fields(:)
TYPE( grid_type ), INTENT(IN)      :: grid
TYPE( um_header_type ), INTENT(IN) :: hdr
INTEGER, INTENT(IN)                :: field_count
INTEGER, INTENT(IN)                :: decomp
LOGICAL, INTENT(OUT)             :: l_changed

! Local variables
INTEGER                            :: pos_sat
INTEGER                            :: pos_wilt
INTEGER                            :: pos_smc
INTEGER                            :: i
INTEGER                            :: j
INTEGER                            :: istat
INTEGER                            :: COUNT

REAL                               :: smc_min
REAL                               :: smc_max

TYPE( field_type ), POINTER        :: soil_moist
TYPE( field_type ), POINTER        :: wilt
TYPE( field_type ), POINTER        :: sat

CHARACTER (LEN=*), PARAMETER :: RoutineName='RCF_SOIL_MOIST_CHK'

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------
! Locate required paramaters in output dump
!-------------------------------------------------------------------

CALL Rcf_Locate ( stashcode_prog_sec, stashcode_vol_smc_sat,  &
                  fields, field_count, pos_sat,  zero_ok_arg = .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_vol_smc_wilt, &
                  fields, field_count, pos_wilt, zero_ok_arg = .TRUE. )
CALL Rcf_Locate ( stashcode_prog_sec, stashcode_soil_moist,   &
                  fields, field_count, pos_smc,  zero_ok_arg = .TRUE. )

!--------------------------------------------------------------------
! Reset smc if appropriate
!--------------------------------------------------------------------

IF (pos_sat  /= 0 .AND. &
    pos_wilt /= 0 .AND. &
    pos_smc  /= 0 ) THEN

  soil_moist => fields( pos_smc  )
  wilt       => fields( pos_wilt )
  sat        => fields( pos_sat  )

  CALL Rcf_Alloc_Field( soil_moist )
  CALL Rcf_Alloc_Field( wilt       )
  CALL Rcf_Alloc_Field( sat        )

  CALL Rcf_Read_Field( soil_moist, hdr, decomp )
  CALL Rcf_Read_Field( wilt,       hdr, decomp )
  CALL Rcf_Read_Field( sat,        hdr, decomp )

  COUNT = 0
  DO i = 1, soil_moist % level_size

    IF ( wilt % DATA(i,1) > 0.0 ) THEN
      ! Non-ice land point
      ! Make sure soil_moist lies between appropriate bounds

      DO j = 1, soil_moist % levels
        ! min is 0.1 * wilt in absolute (not volumetric) terms
        smc_min = 0.1 * wilt % DATA(i,1) * grid % soil_depths(j) &
                      * rho_water

        ! max is saturation in absolute terms
        smc_max = sat % DATA(i,1) * grid % soil_depths(j) * rho_water

        IF ( soil_moist % DATA(i,j) < smc_min) THEN
          ! Too low - set to the minimum value
          soil_moist % DATA(i,j) = smc_min
          COUNT = COUNT + 1

        ELSE IF ( soil_moist % DATA(i,j) > smc_max) THEN
          ! Too high - set to the maximum value
          soil_moist % DATA(i,j) = smc_max
          COUNT = COUNT + 1

        END IF
      END DO

    ELSE IF (wilt % DATA(i,1) == 0.0 ) THEN
      ! land-ice point
      ! soil_moist should be set to 0.0 if it isn't
      DO j = 1, soil_moist % levels
        IF ( soil_moist % DATA(i,j) /= 0.0 ) THEN
          soil_moist % DATA(i,j) = 0.0
          COUNT = COUNT + 1
        END IF
      END DO
    END IF

  END DO

  CALL gc_isum (1, nproc, istat, COUNT)
  IF (PrintStatus >= PrStatus_Normal .AND. mype == 0) THEN
    WRITE(umMessage,*) 'Soil Moisture: No of values reset ', &
                 COUNT
    CALL umPrint(umMessage,src='rcf_soil_moist_chk_mod')
  END IF


  ! Set the changed flag to true if count is positive
  IF (COUNT > 0) THEN
    l_changed = .TRUE.
  ELSE
    l_changed = .FALSE.
  END IF

  !---------------------------------------------------------------------
  ! Write out changed field
  !---------------------------------------------------------------------

  IF (l_changed) THEN
    CALL Rcf_Write_Field( soil_moist, hdr, decomp )
  END IF

  CALL Rcf_Dealloc_Field( soil_moist )
  CALL Rcf_Dealloc_Field( wilt       )
  CALL Rcf_Dealloc_Field( sat        )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN
END SUBROUTINE Rcf_Soil_Moist_Chk
END MODULE Rcf_Soil_Moist_Chk_Mod
