! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Convert volumetric soil moisture to soil moisture amount

MODULE Rcf_Grib_Spcl_SoilM_Mod

! SUBROUTINE Rcf_Grib_Spcl_SoilM
!
! Description: This routine handles the conversion of volumetric
!              soil moisture into soil moisture amounts.
!              At present there is a double check to ensure the data
!              came from ECMWF as not all input data will be volumetric.
!
! Method: Multiply volumetric quantity by density and layer depth.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_SPCL_SOILM_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Spcl_SoilM(Current,FieldData)

USE Rcf_GRIB_Block_Params_Mod, ONLY: &
    Grib_Record,     &
    LenArrayMax,     &
    p_Orig_cntr,     &
    p_Lvl_Desc_1,    &
    p_Lvl_Desc_2

USE rcf_GRIB_lookups_Mod, ONLY: &
    GrbOrigECMWF

USE umPrintMgr, ONLY:      &
    umPrint,                &
    umMessage,              &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

USE UM_ParCore, ONLY: &
    mype

USE EReport_Mod, ONLY:     &
    EReport

USE water_constants_mod, ONLY: rho_water

USE errormessagelength_mod, ONLY: errormessagelength

USE yomhook,   ONLY: & ! DrHook
    lhook,           &
    dr_hook

USE parkind1,  ONLY: &
    jprb,            &
    jpim

IMPLICIT NONE

! Subroutine arguments

!< Array  arguments with intent(InOut):>
TYPE (Grib_Record),POINTER          :: Current
REAL, INTENT(INOUT)                 :: FieldData(LenArrayMax)

! Local constants
CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_SPCL_SOILM'

! Local variables

CHARACTER (LEN=errormessagelength) :: Cmessage    ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport
INTEGER                          :: ilayer        ! Index of soil layer
REAL                             :: layer_depth   ! Depth of soil layer

! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=======================================================================
!  Routine Code Start :
!=======================================================================

! double check it's an ECMWF field
IF (Current % Block_1(p_Orig_cntr) == GrbOrigECMWF) THEN

  IF ( PrintStatus >= PrStatus_Diag  ) THEN
    IF ( mype == 0 ) THEN
      CALL umPrint( "Converting vol. SoilM to SoilM Amount", &
          src='rcf_grib_spcl_soilm_mod')
    END IF
  END IF

  ! Hard wired for time being - need to derive these some how
  ! Note: Although Block_1(p_Lvl_Desc_1/2) have been modified
  !       in rcf_grib_check, the second read of the data
  !       overwrites the existing header information. Hence
  !       the values are as they are in the original grib file.
  SELECT CASE (Current % Block_1(p_Lvl_Desc_1))
  CASE (0)
    layer_depth = 0.10
  CASE (7)
    layer_depth = 0.25
  CASE (28)
    layer_depth = 0.65
  CASE (100)
    layer_depth = 2.00
  END SELECT

  ! Multiply volumetric values by density and layer depth
  FieldData(:) = layer_depth * rho_water * FieldData(:)

ELSE
  WRITE (Cmessage,*) 'Data did not come from ECMWF'
  ErrorStatus = 10
  CALL EReport( RoutineName, ErrorStatus, Cmessage)
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Spcl_SoilM
END MODULE Rcf_Grib_Spcl_SoilM_Mod
