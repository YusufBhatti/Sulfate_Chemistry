! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Convert log(pstar) into pstar

MODULE Rcf_Grib_Spcl_LPstar_Mod
IMPLICIT NONE

! SUBROUTINE Rcf_Grib_Spcl_LPstar
!
! Description: This routine handles the conversion of log(pstar) to
!              pstar.
!              At present there is a double check to ensure the data
!              came from ECMWF as log(pstar) is not a standard variable.
!
! Method: Take the exponential of every entry in the data field.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Reconfiguration
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.
!
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='RCF_GRIB_SPCL_LPSTAR_MOD'

CONTAINS
SUBROUTINE Rcf_Grib_Spcl_LPstar(Current,FieldData)

USE Rcf_GRIB_Block_Params_Mod, ONLY: &
    Grib_Record,     &
    LenArrayMax,     &
    p_Orig_cntr,     &
    p_Param_ID,      &
    EID_Surf_Press

USE rcf_GRIB_lookups_Mod, ONLY: &
    GrbOrigECMWF

USE umPrintMgr, ONLY:      &
    umPrint,                &
    PrintStatus,            &
    PrStatus_Normal,        &
    PrStatus_Diag

USE UM_ParCore, ONLY: &
    mype

USE EReport_Mod, ONLY:     &
    EReport

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
CHARACTER (LEN=*), PARAMETER     :: RoutineName='RCF_GRIB_SPCL_LPSTAR'

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

! double check it's an ECMWF field
IF (Current % Block_1(p_Orig_cntr) == GrbOrigECMWF) THEN

  IF ( PrintStatus >= PrStatus_Diag  ) THEN
    IF ( mype == 0 ) THEN
      CALL umPrint( "Converting log(pstar) to pstar", &
          src='rcf_grib_spcl_lpstar_mod')
    END IF
  END IF
  FieldData(:) = EXP( FieldData(:) )
  Current % Block_1(p_Param_ID) = EID_Surf_Press
ELSE
  WRITE (Cmessage,*) 'Data did not come from ECMWF'
  ErrorStatus = 10
  CALL EReport( RoutineName, ErrorStatus, Cmessage )
END IF

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName, zhook_out, zhook_handle)
RETURN

END SUBROUTINE Rcf_Grib_Spcl_LPstar
END MODULE Rcf_Grib_Spcl_LPstar_Mod
