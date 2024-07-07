! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Convert geopotential field to Orography

MODULE Rcf_Grib_Spcl_Orog_Mod

! SUBROUTINE Rcf_Grib_Spcl_Orog
!
! Description: This routine handles the conversion of Geopotential to
!              Orography.
!              At present there is a double check to ensure the data
!              came from ECMWF but if all centers return geopotential
!              instead of Orography then this check can be removed.
!
! Method: Divide every entry in the data field by G, the gravitational
!         constant.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_SPCL_OROG_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Spcl_Orog(Current,FieldData)

USE Rcf_GRIB_Block_Params_Mod, ONLY: &
    Grib_Record,     &
    LenArrayMax,     &
    p_Orig_cntr

USE rcf_GRIB_lookups_Mod, ONLY: &
    GrbOrigECMWF

USE umPrintMgr, ONLY:       &
    umPrint,                &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

USE UM_ParCore, ONLY: &
    mype

USE EReport_Mod, ONLY:     &
    EReport

USE planet_constants_mod, ONLY: g

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
CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_SPCL_OROG'

! Local variables

CHARACTER (LEN=errormessagelength)  :: Cmessage   ! used for EReport
INTEGER                          :: ErrorStatus   ! used for EReport
! Dr Hook
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!=======================================================================
!  Routine Code Start :
!=======================================================================

! double check it's an ECMWF field (do other sources spec geopot ?)
IF (Current % Block_1(p_Orig_cntr) == GrbOrigECMWF) THEN

  IF ( PrintStatus >= PrStatus_Diag  ) THEN
    IF ( mype == 0 ) THEN
      CALL umPrint( "Converting Geopotential to Orography", &
          src='rcf_grib_spcl_orog_mod')
    END IF
  END IF
  FieldData(:) = FieldData(:) / g
ELSE
  WRITE (Cmessage,*) 'Data did not come from ECMWF'
  ErrorStatus = 10
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Spcl_Orog
END MODULE Rcf_Grib_Spcl_Orog_Mod
